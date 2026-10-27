from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,to_rgba
from matplotlib.patches import Patch
from PIL import Image
P=Path('C:/Users/31937/CAF_paper_figures_20261021/CK');O=Path('C:/Users/31937/Downloads/CAF_Fig1c_Fig1k_600DPI_20261021')
O.mkdir(parents=True,exist_ok=True);(P/'previews').mkdir(exist_ok=True);(P/'QA').mkdir(exist_ok=True)
cfg=json.loads((P/'PANELS.json').read_text(encoding='utf-8'));c=pd.read_csv(P/'C_dotplot_values.tsv',sep='\t');profiles=pd.read_csv(P/'K_gene_patient_class_profiles.tsv',sep='\t');program=pd.read_csv(P/'K_program_patient_class_profiles.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7,'savefig.facecolor':'white'})
types=['T cell','Macrophage','Monocyte / DC unresolved','pDC-like','Mixed B / plasma','Plasma-like','Mast','Endothelial','Lymphatic endothelium','Fibroblast-like','Mural / pericyte','Osteogenic-like','Chondrocyte-like','Mesenchymal unresolved','Myogenic','Osteoclast-like','Epithelial','Renal epithelial-like','Erythroid','Mixed endothelial / epithelial','Unresolved lineage']
classcolors={'F1':'#477DA8','F2':'#5A9B86','F3':'#9D82BC','F4':'#D27A79','F5':'#C3A353'};sitecolors={'Bone':'#477DA8','Lung':'#D27A79'}
groupshort={'Fibroblast':'Fibro','Mural':'Mural','Osteogenic':'Osteo','T_NK':'T/NK','B_Plasma':'B/Plasma','Myeloid':'Myeloid','Endothelial':'Endo','Epithelial':'Epi','Mast':'Mast','Myogenic':'Myo','Osteoclast':'OC','Chondrocyte':'Chondro','Erythroid':'Ery'}
grouplabel={'collagen_components':'Collagen components','collagen_biosynthesis':'Collagen processing','collagen_crosslinking':'Collagen crosslinking','ECM_degradation':'ECM remodeling','contraction':'Contraction','mineralization':'Mineralization','inflammation':'Inflammation','Candidates':'Six candidates'}
blockcolors=dict(zip(cfg['K'],['#88AADD','#6898BF','#7DA89A','#7AAB79','#BE9C69','#A68BBB','#D58E85','#657985']))
fix_layout='--fix-layout' in sys.argv
manifest=[f for f in json.loads((P/'QA/MANIFEST.json').read_text(encoding='utf-8'))['figures'] if f['kind']=='CAF_programs_heatmap'] if fix_layout else []
def zscore(x,axis):
 mean=np.nanmean(x,axis=axis,keepdims=True);sd=np.nanstd(x,axis=axis,keepdims=True);z=np.divide(x-mean,sd,out=np.zeros_like(x),where=sd>1e-12);z[~np.isfinite(x)]=np.nan;return z
def save(fig,name,kind,src):
 exports=[]
 for ext in ['png','tiff']:
  f=O/(name+'.'+ext);fig.savefig(f,dpi=600,**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
  with Image.open(f) as im:
   im.load();assert im.size==tuple(round(x*600) for x in fig.get_size_inches());dpi=list(map(float,im.info['dpi']));assert max(abs(v-600) for v in dpi)<.05
   exports.append(dict(file=str(f),pixels=im.size,dpi=dpi,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 fig.savefig(P/'previews'/(name+'.png'),dpi=110);plt.close(fig);manifest.append(dict(name=name,kind=kind,source=src,exports=exports));print('SAVED',name,flush=True)
for src in ['SCPCP000017','SCPCP000023']:
 cc=c[c.source.eq(src)];labels=[t for t in types if t in set(cc.cell_type)];genes=[g for v in cfg['C'].values() for g in v]
 means=cc.pivot(index='cell_type',columns='gene',values='mean_log1pCP10k').loc[labels,genes];det=cc.pivot(index='cell_type',columns='gene',values='detection_percent').loc[labels,genes];zs=zscore(means.to_numpy(),axis=0)
 pd.DataFrame(zs,index=labels,columns=genes).to_csv(P/'QA'/f'{src}_C_color_zscores.tsv',sep='\t')
 fig=plt.figure(figsize=(16.8,9.5),dpi=110);ax=fig.add_axes([.235,.28,.71,.51]);norm=Normalize(-2,2);cmap=plt.get_cmap('RdBu_r')
 yy,xx=np.indices(zs.shape);valid=np.isfinite(zs)&np.isfinite(det.to_numpy())
 ax.scatter(xx[valid],yy[valid],s=det.to_numpy()[valid]*2.1,c=zs[valid],cmap=cmap,norm=norm,lw=.2,edgecolor='#7C858B')
 if (~valid).any():ax.scatter(xx[~valid],yy[~valid],marker='x',s=18,color='#8F969D')
 ns=cc.groupby('cell_type').patients.first();ax.set_yticks(range(len(labels)),[f'{t}  (n={ns[t]})' for t in labels],fontsize=9)
 ax.set_xticks(range(len(genes)),genes,rotation=65,ha='right',fontstyle='italic',fontsize=8.5);ax.set_xlim(-.6,len(genes)-.4);ax.set_ylim(len(labels)-.5,-.6);ax.tick_params(length=0,pad=6)
 for y in range(len(labels)):ax.axhline(y,color='#EEF0F2',lw=.5,zorder=0)
 offset=0
 for group,gs in cfg['C'].items():
  ax.text(offset+(len(gs)-1)/2,-1.0,groupshort[group],ha='center',fontsize=8.5,clip_on=False)
  if offset:ax.axvline(offset-.5,color='#BBC3CA',lw=.65)
  offset+=len(gs)
 for spine in ax.spines.values():spine.set_visible(False)
 fig.text(.055,.94,f'{src} | Marker expression across reviewed cell types',fontsize=18,fontweight='bold')
 fig.text(.055,.89,'All-cell scope | 39 frozen-panel markers | patient-balanced expression and RNA detection',fontsize=12,color='#53616B')
 cbax=fig.add_axes([.235,.115,.21,.018]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbax,orientation='horizontal',extend='both');cb.ax.tick_params(labelsize=9)
 fig.text(.235,.155,'Mean expression: gene-wise z-score',fontsize=10)
 lax=fig.add_axes([.56,.075,.29,.09]);lax.set_xlim(-.4,3.8);lax.set_ylim(-.3,1);lax.axis('off')
 for j,pct in enumerate([10,25,50,100]):lax.scatter(j,.55,s=pct*2.1,c='#A2ACB4',edgecolor='#727C83',lw=.2);lax.text(j,-.02,str(pct)+'%',ha='center',fontsize=9)
 fig.text(.56,.175,'RNA-detected cells (%)',fontsize=10)
 fig.text(.055,.035,'n: patients with that cell type. Frozen expression-lineage labels; marker display does not establish malignancy or confirmed CAF identity.',fontsize=10,color='#53616B')
 save(fig,f'Fig1c_{src}_allcell_marker_dotplot','allcell_marker_dotplot',src)
 for mode,table in [('genes',profiles),('programs',program)]:
  if fix_layout and mode=='programs':continue
  q=table[table.source.eq(src)].sort_values(['class','donor']).reset_index(drop=True)
  if mode=='genes':features=[g for v in cfg['K'].values() for g in v];rowlabels=features
  else:features=['raw__'+k for k in cfg['programs']['sets']];rowlabels=[k.removeprefix('raw__').replace('_',' ') for k in features]
  raw=q[features].to_numpy().T;z=zscore(raw,axis=1);assert np.isfinite(raw).all()
  pd.DataFrame(z,index=rowlabels,columns=q['class']+'|'+q.donor).to_csv(P/'QA'/f'{src}_K_{mode}_zscores.tsv',sep='\t')
  fig=plt.figure(figsize=(14.4,10.8 if mode=='genes' else 8.6),dpi=110)
  left=.19 if mode=='genes' else .25;right=.91;width=right-left;bottom=.28;height=.51
  ax=fig.add_axes([left,bottom,width,height]);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#D9DDE0');im=ax.imshow(z,aspect='auto',interpolation='nearest',cmap=cmap,vmin=-2.5,vmax=2.5)
  ax.set_yticks(range(len(features)),rowlabels,fontsize=9,fontstyle='italic' if mode=='genes' else 'normal')
  ax.set_xticks(range(len(q)),[s.split('|')[-1] for s in q.donor],rotation=90,fontsize=7)
  ax.tick_params(length=0,pad=4);ax.spines[['top','right','bottom','left']].set_visible(False)
  for i in range(1,len(q)):
   if q['class'].iloc[i]!=q['class'].iloc[i-1]:ax.axvline(i-.5,color='white',lw=1.3)
  bar=fig.add_axes([left,bottom+height+.012,width,.027]);rgba=np.array([[to_rgba(classcolors[v]) for v in q['class']],[to_rgba(sitecolors[v]) for v in q.site]])
  bar.imshow(rgba,aspect='auto',interpolation='nearest');bar.set_yticks([0,1],['Class','Site'],fontsize=8);bar.set_xticks([]);bar.tick_params(length=0);bar.spines[:].set_visible(False)
  for cl,gg in q.groupby('class',sort=False):
   mid=float(gg.index.to_numpy().mean());bar.text(mid,-1.3,f'{cl} (n={len(gg)})',ha='center',va='bottom',fontsize=9,clip_on=False)
  if mode=='genes':
   rs=fig.add_axes([.095,bottom,.006,height]);rgb=np.array([to_rgba(blockcolors[k]) for k,gs in cfg['K'].items() for _ in gs]).reshape(len(features),1,4);rs.imshow(rgb,aspect='auto');rs.axis('off')
   pos=0
   for group,gs in cfg['K'].items():
    if pos:ax.axhline(pos-.5,color='white',lw=1.0)
    pos+=len(gs)
   fig.legend(handles=[Patch(facecolor=blockcolors[k],label=grouplabel[k]) for k in cfg['K']],loc='lower left',bbox_to_anchor=(.08,.035),ncol=4,frameon=False,fontsize=9)
  cbax=fig.add_axes([.935,bottom,.012,height]);cb=fig.colorbar(im,cax=cbax,extend='both');cb.set_label('Row z-score',fontsize=10);cb.ax.tick_params(labelsize=8)
  title='Functional genes and six candidates' if mode=='genes' else 'Frozen functional expression programs'
  fig.text(.055,.95,f'{src} | {title}',fontsize=17,fontweight='bold')
  fig.text(.055,.91,f'{int(q.cells.sum()):,} candidate cells | {q.donor.nunique()} patients | {len(q)} patient-class profiles',fontsize=11,color='#53616B')
  fig.text(.055,.875,'Columns: patient-class means, ordered by existing F class; site stripe: blue=Bone, red=Lung',fontsize=10,color='#53616B')
  fig.text(.08,.018,'Descriptive RNA profiles; repeated patients across classes are not independent replicates. Colors compare profiles within each row.',fontsize=9,color='#53616B')
  if mode=='programs':fig.text(.08,.085,'Fixed gene-set means of log1p(CP10k); transcriptional scores do not directly measure functional activity.',fontsize=10,color='#53616B')
  save(fig,f'Fig1k_{src}_CAF_{mode}_heatmap',f'CAF_{mode}_heatmap',src)
(P/'QA/MANIFEST.json').write_text(json.dumps(dict(status='export_complete_visual_pending',figures=manifest,total=6),indent=2),encoding='utf-8')
print('ALL_SIX_RENDERED',flush=True)
