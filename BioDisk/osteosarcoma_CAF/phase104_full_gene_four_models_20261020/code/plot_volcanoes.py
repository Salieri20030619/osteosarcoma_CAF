import pathlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase104_full_gene_four_models_20261020')
summary=pd.read_csv(P/'tables/four_analysis_summary.tsv',sep='\t').set_index('model')
titles={'pooled_tumor':'A. Pooled unpaired comparison','source_adjusted_tumor':'B. Source-adjusted comparison','pooled_normal_min20':'C. Pooled + normal-background subtraction','source_adjusted_normal_min20':'D. Source-adjusted + normal-background subtraction'}
colors={'FDR<0.05':'#C43C39','P<0.05_only':'#DA961F','not_significant':'#B8C1CA'}
qa=[]
for model,title in titles.items():
 d=pd.read_csv(P/'tables'/(model+'_allgenes.tsv.gz'),sep='\t')
 z=d[d.valid_test].copy();z['y']=-np.log10(z.P.clip(lower=np.nextafter(0.,1.)))
 fig,ax=plt.subplots(figsize=(11,5.6),dpi=160)
 for cls in ['not_significant','P<0.05_only','FDR<0.05']:
  q=z[z.significance.eq(cls)]
  ax.scatter(q.effect,q.y,s=6,c=colors[cls],alpha=.65 if cls=='not_significant' else .8,linewidths=0,rasterized=True)
 ax.axhline(-np.log10(.05),color='#646e75',ls='--',lw=.8);ax.axvline(0,color='#78848d',lw=.7)
 span=max(abs(z.effect.min()),abs(z.effect.max()))*1.07;ax.set_xlim(-span,span)
 ax.set_ylim(-.15,max(z.y.max()*1.13,3.5));ax.set_ylabel('-log10(raw P)',fontsize=10)
 ax.set_xlabel('Lung - Bone effect on log2(CP10k + 0.1) scale'+(' (normal background subtracted)' if 'normal' in model else ''),fontsize=10)
 ax.set_title(title,loc='left',fontsize=13,fontweight='bold',pad=12)
 handles=[Line2D([0],[0],marker='o',ls='',color=colors[k],markersize=6,label=v) for k,v in [('FDR<0.05','FDR < 0.05'),('P<0.05_only','P < 0.05 only'),('not_significant','Not significant')]]
 ax.legend(handles=handles,loc='upper left',frameon=True,facecolor='white',framealpha=.92,fontsize=9)
 s=summary.loc[model]
 ax.text(.99,.97,f"Tested: {int(s.valid_tests):,}   P < 0.05: {int(s.nominal):,}   FDR < 0.05: {int(s.FDR_significant):,}",transform=ax.transAxes,ha='right',va='top',fontsize=9)
 target=z[z.is_original14].copy()
 ax.scatter(target.effect,target.y,s=44,facecolors='none',edgecolors='#172B45',linewidths=1.0,zorder=4)
 fig.subplots_adjust(left=.075,right=.98,bottom=.16,top=.88)
 fig.canvas.draw()
 # Place 14 labels in pixel coordinates, repel text boxes, keep in axes.
 renderer=fig.canvas.get_renderer();ann=[];anchors=[]
 for i,r in enumerate(target.sort_values('y').itertuples()):
  xy=ax.transData.transform((r.effect,r.y));anchors.append(xy)
  dx=35 if i%2==0 else -35;dy=13+(i%3)*8
  pos=ax.transData.inverted().transform(xy+[dx,dy])
  text=ax.text(*pos,r.gene,fontsize=8.5,ha='center',va='center',zorder=6,bbox=dict(boxstyle='round,pad=0.13',fc='white',ec='none',alpha=.88))
  ann.append(text)
 limits=ax.get_window_extent(renderer)
 for iteration in range(350):
  fig.canvas.draw();renderer=fig.canvas.get_renderer()
  boxes=[t.get_window_extent(renderer).expanded(1.10,1.18) for t in ann];moves=np.zeros((len(ann),2));overlap=False
  centers=np.array([ax.transData.transform(t.get_position()) for t in ann])
  for i in range(len(ann)):
   for j in range(i):
    if boxes[i].overlaps(boxes[j]):
     overlap=True
     delta=centers[i]-centers[j]
     if abs(delta[1])<1:delta[1]=1 if i%2 else -1
     shift=np.array([np.sign(delta[0] or 1)*.6,np.sign(delta[1])*2.2])
     moves[i]+=shift;moves[j]-=shift
   for j,anchor in enumerate(anchors):
    if boxes[i].contains(*anchor):moves[i,1]+=2.0;overlap=True
   if boxes[i].y0<limits.y0+8:moves[i,1]+=limits.y0+8-boxes[i].y0
   if boxes[i].y1>limits.y1-28:moves[i,1]-=boxes[i].y1-(limits.y1-28)
   if boxes[i].x0<limits.x0+6:moves[i,0]+=limits.x0+6-boxes[i].x0
   if boxes[i].x1>limits.x1-6:moves[i,0]-=boxes[i].x1-(limits.x1-6)
  for i,t in enumerate(ann):t.set_position(ax.transData.inverted().transform(centers[i]+moves[i]))
  if not overlap and np.max(abs(moves))<.1:break
 fig.canvas.draw();boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in ann]
 collisions=sum(boxes[i].overlaps(boxes[j]) for i in range(len(boxes)) for j in range(i))
 for t,xy in zip(ann,anchors):
  targetxy=ax.transData.inverted().transform(xy)
  ax.annotate('',xy=targetxy,xytext=t.get_position(),arrowprops=dict(arrowstyle='-',color='#526273',lw=.55),zorder=5)
 ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=9)
 fig.savefig(P/'figures'/(model+'_volcano.png'),dpi=220);plt.close(fig)
 qa.append(dict(model=model,points=len(z),labelled_candidates=len(target),label_collisions=int(collisions),plotted_FDR=int(z.FDR.lt(.05).sum()),plotted_P_only=int((z.P.lt(.05)&~z.FDR.lt(.05)).sum())))
 print(qa[-1],flush=True)
(P/'resources/VOLCANO_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
