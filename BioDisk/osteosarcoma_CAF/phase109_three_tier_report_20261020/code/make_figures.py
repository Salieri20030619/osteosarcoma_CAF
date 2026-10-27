import pathlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase109_three_tier_report_20261020')
specs=json.loads((P/'resources/FIGURE_SPECS.json').read_text());six=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
titles={'tier1_source17':'Tier 1 | SCPCP000017: Lung 3 vs Bone 12','tier1_source23':'Tier 1 | SCPCP000023: Lung 3 vs Bone 4','tier2_pooled':'Tier 2 | Pooled Lung 23 vs Bone 16','tier2_NCI':'Tier 2 | NCI Lung 17 vs Bone 16','tier3_pooled':'Tier 3 | Pooled contrast minus normal background','tier3_NCI':'Tier 3 | NCI contrast minus normal background'}
colors={'FDR<0.05':'#C43C39','P<0.05_only':'#DA961F','not_significant':'#B8C1CA'};qa=[]
for spec in specs:
 key=spec['key'];d=pd.read_csv(P/'tables'/(key+'_volcano_source.tsv.gz'),sep='\t');d=d[d.P.notna()].copy()
 assert len(d)==spec['valid'];d['y']=-np.log10(d.P.clip(lower=np.nextafter(0.,1.)))
 d['sig']=np.select([d.FDR.lt(.05),d.P.lt(.05)],['FDR<0.05','P<0.05_only'],default='not_significant')
 fig,ax=plt.subplots(figsize=(11,5.6),dpi=160)
 for c in ['not_significant','P<0.05_only','FDR<0.05']:
  z=d[d.sig.eq(c)];ax.scatter(z.effect,z.y,s=6,c=colors[c],alpha=.65 if c=='not_significant' else .8,lw=0,rasterized=True)
 ax.axhline(-np.log10(.05),color='#64727c',ls='--',lw=.8);ax.axvline(0,color='#7d8791',lw=.7)
 span=max(abs(d.effect.min()),abs(d.effect.max()))*1.07;ax.set_xlim(-span,span)
 ax.set_ylim(-.12,max(d.y.max()*1.2+.2,2.2));ax.set_ylabel('-log10(raw P)',fontsize=10)
 ax.set_xlabel('Lung - Bone effect on log2(CP10k + 0.1) scale'+(' (background subtracted)' if spec['stage']==3 else ''),fontsize=10)
 ax.set_title(titles[key],loc='left',fontsize=13,fontweight='bold',pad=12)
 handles=[Line2D([0],[0],marker='o',ls='',color=colors[k],markersize=6,label=v) for k,v in [('FDR<0.05','FDR < 0.05'),('P<0.05_only','P < 0.05 only'),('not_significant','Not significant')]]
 ax.legend(handles=handles,loc='upper left',frameon=True,facecolor='white',framealpha=.94,fontsize=9)
 ax.text(.99,.97,f"Tested: {spec['valid']:,}   P < 0.05: {spec['nominal']:,}   FDR < 0.05: {spec['FDR']:,}",transform=ax.transAxes,ha='right',va='top',fontsize=9)
 z=d[d.gene.isin(six)].copy();assert len(z)==6
 ax.scatter(z.effect,z.y,s=48,facecolors='none',edgecolors='#172B45',linewidths=1.1,zorder=4)
 fig.subplots_adjust(left=.075,right=.98,bottom=.16,top=.88);fig.canvas.draw()
 texts=[];anchors=[]
 for i,r in enumerate(z.sort_values('effect').itertuples()):
  pt=ax.transData.transform((r.effect,r.y));anchors.append(pt)
  pos=ax.transData.inverted().transform(pt+[(-1 if i%2 else 1)*35,17+i%3*12])
  texts.append(ax.text(*pos,r.gene,fontsize=9,ha='center',va='center',zorder=6,bbox=dict(boxstyle='round,pad=0.16',fc='white',ec='none',alpha=.9)))
 limits=ax.get_window_extent(fig.canvas.get_renderer())
 for step in range(450):
  fig.canvas.draw();rend=fig.canvas.get_renderer();boxes=[t.get_window_extent(rend).expanded(1.12,1.25) for t in texts]
  centers=np.array([ax.transData.transform(t.get_position()) for t in texts]);moves=np.zeros((len(texts),2));overlap=False
  for i in range(len(texts)):
   for j in range(i):
    if boxes[i].overlaps(boxes[j]):
     overlap=True;delta=centers[i]-centers[j];delta[1]=delta[1] if abs(delta[1])>=1 else (1 if i%2 else -1)
     shift=np.array([np.sign(delta[0] or 1)*.8,np.sign(delta[1])*2.4]);moves[i]+=shift;moves[j]-=shift
   if any(boxes[i].contains(*pt) for pt in anchors):moves[i,1]+=2;overlap=True
   if boxes[i].y0<limits.y0+7:moves[i,1]+=limits.y0+7-boxes[i].y0
   if boxes[i].y1>limits.y1-30:moves[i,1]-=boxes[i].y1-(limits.y1-30)
   if boxes[i].x0<limits.x0+6:moves[i,0]+=limits.x0+6-boxes[i].x0
   if boxes[i].x1>limits.x1-6:moves[i,0]-=boxes[i].x1-(limits.x1-6)
  for i,t in enumerate(texts):t.set_position(ax.transData.inverted().transform(centers[i]+moves[i]))
  if not overlap and np.max(abs(moves))<.1:break
 fig.canvas.draw();boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in texts]
 collisions=sum(boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i))
 for t,pt in zip(texts,anchors):ax.annotate('',xy=ax.transData.inverted().transform(pt),xytext=t.get_position(),arrowprops=dict(arrowstyle='-',color='#526273',lw=.55),zorder=5)
 ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=9)
 fig.savefig(P/'figures'/(key+'.png'),dpi=220);plt.close(fig)
 qa.append(dict(key=key,points=len(d),labels=6,label_collisions=int(collisions)));print(qa[-1],flush=True)
(P/'resources/PLOT_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
