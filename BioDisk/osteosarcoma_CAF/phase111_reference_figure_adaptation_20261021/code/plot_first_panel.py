from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from PIL import Image
P=Path('C:/Users/31937/CAF_paper_figures_20261021')
O=Path('C:/Users/31937/Downloads/CAF_reference_figures_600DPI_20261021');O.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(P/'Fig1i_six_candidates_dotplot.tsv',sep='\t')
genes=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
sources=['SCPCP000017','SCPCP000023','NCI']
groups=[(s,c) for s in sources for c in sorted(d.loc[d.source.eq(s),'group'].unique())]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7,'savefig.facecolor':'white','svg.fonttype':'none'})
fig=plt.figure(figsize=(8.8,6.8),dpi=140)
ax=fig.add_axes([.22,.29,.46,.54])
cmap=LinearSegmentedColormap.from_list('expression',['#E6F1F4','#F7C1BE','#CC4052','#7A1236'])
norm=Normalize(0,float(d.mean_log1pCP10k.max()))
for y,(src,cl) in enumerate(groups):
 q=d[(d.source==src)&(d.group==cl)].set_index('gene').loc[genes]
 ax.scatter(np.arange(6),np.repeat(y,6),s=q.RNA_detected_percent.to_numpy()*4.5,c=q.mean_log1pCP10k.to_numpy(),cmap=cmap,norm=norm,edgecolors='#747474',linewidths=.25)
 ax.text(6.03,y,str(int(q.donors.iloc[0])),va='center',ha='center',fontsize=9,clip_on=False)
 ax.text(6.95,y,f'{int(q.cells.iloc[0]):,}',va='center',ha='center',fontsize=9,clip_on=False)
ax.text(6.03,-1.1,'Patients',ha='center',fontsize=9,clip_on=False)
ax.text(6.95,-1.1,'Cells',ha='center',fontsize=9,clip_on=False)
ax.set_xlim(-.6,5.6);ax.set_ylim(len(groups)-.5,-.75)
ax.set_xticks(np.arange(6),genes,rotation=38,ha='right',fontsize=10,fontstyle='italic')
short={'SCPCP000017':'017','SCPCP000023':'023','NCI':'NCI'}
ax.set_yticks(np.arange(len(groups)),[short[s]+'  |  '+c for s,c in groups],fontsize=10)
ax.tick_params(length=0,pad=8)
for j in range(6):ax.axvline(j,color='#EDF0F2',lw=.55,zorder=0)
for y in range(len(groups)):ax.axhline(y,color='#F0F1F3',lw=.55,zorder=0)
for y in range(1,len(groups)):
 if groups[y][0]!=groups[y-1][0]:ax.axhline(y-.5,color='#9AA4AD',lw=.8)
for s in ax.spines.values():s.set_visible(False)
fig.text(.07,.94,'Six candidate genes across CAF candidate groups',fontsize=15,fontweight='bold')
fig.text(.07,.895,'Patient-balanced expression  |  sources displayed separately',fontsize=10,color='#53616B')
cax=fig.add_axes([.22,.105,.23,.017]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal');cb.ax.tick_params(labelsize=8,length=2);cb.outline.set_linewidth(.5)
fig.text(.22,.145,'Mean log1p(CP10k)',fontsize=9)
lax=fig.add_axes([.53,.075,.37,.08]);lax.set_xlim(0,4);lax.set_ylim(0,1);lax.axis('off')
for i,v in enumerate([10,25,50,75]):
 lax.scatter(i+.3,.56,s=v*4.5,color='#B2BBC2',edgecolor='#747474',lw=.25)
 lax.text(i+.3,-.06,str(v)+'%',ha='center',fontsize=8)
fig.text(.53,.17,'RNA-detected cells (%)',fontsize=9)
fig.text(.07,.025,'017 / 023: SCPCP000017 / SCPCP000023; NCI: author-annotated CAFs.',fontsize=8,color='#53616B')
name='Fig1i_six_candidate_expression_dotplot'
for ext in ['png','tiff']:
 kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
 fig.savefig(O/(name+'.'+ext),dpi=600,**kw)
fig.savefig(P/(name+'_preview.png'),dpi=130)
plt.close(fig)
checks=[]
for ext in ['png','tiff']:
 f=O/(name+'.'+ext)
 with Image.open(f) as im:
  assert im.size==(5280,4080);dpi=im.info.get('dpi');assert dpi and all(abs(v-600)<.05 for v in dpi)
  im.load();checks.append(dict(file=str(f),pixels=im.size,dpi=[float(v) for v in dpi],sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
(P/'Fig1i_export_QA.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks))
