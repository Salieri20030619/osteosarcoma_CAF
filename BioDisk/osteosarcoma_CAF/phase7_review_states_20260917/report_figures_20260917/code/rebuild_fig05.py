from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();root=p/'report_figures_20260917';out=root/'figures';out.mkdir(exist_ok=True);rng=np.random.default_rng(20260917)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white'})
cols=['#327DA8','#DE8650','#6BA484','#8F79A9','#C6A255','#5FA6A5']
E=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
def save(fig,name,caption):
 fig.savefig(out/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white');fig.savefig(out/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig);captions[name]=caption;print('FIGURE',name,flush=True)
def tidy(ax):ax.set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2');ax.set_aspect('equal',adjustable='datalim')
def catplot(ax,xy,lab,palette,title,size=3):
 for k,c in palette.items():
  ix=np.asarray(lab)==k
  if ix.any():ax.scatter(xy[ix,0],xy[ix,1],s=size,c=c,label=k,rasterized=True,linewidths=0,alpha=.85)
 tidy(ax);ax.set_title(title,pad=12)
def group(e):
 z=pd.Series('Other / unresolved',index=e.index)
 z.loc[e.lineage_reviewed.isin(['T_NK','B_Plasma','Myeloid','DC','Mast'])]='Immune'
 z.loc[e.lineage_reviewed=='Endo']='Endothelial'
 z.loc[e.lineage_reviewed=='OC']='Osteoclast'
 z.loc[e.final_identity=='Mural_SMC']='Mural'
 z.loc[e.final_identity=='Malignant_OS_candidate']='Malignant candidate'
 z.loc[e.caf_strict]='CAF candidate';return z
pal={'Other / unresolved':'#d5d8dc','Immune':'#6BA484','Endothelial':'#8F79A9','Osteoclast':'#C6A255','Mural':'#5FA6A5','CAF candidate':'#327DA8','Malignant candidate':'#D26D5B'}
captions={}
# CAF state visualizations.
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');old=sc.read_h5ad(p/'analysis/GSE152048_strict_CAF_states.h5ad');patients=a.obs.patient_id.str.replace('GSE152048:','',regex=False).to_numpy();ps=sorted(set(patients));pp=dict(zip(ps,cols));xy=a.obsm['X_umap_harmony_theta2'];lab=a.obs['harmony_theta2_r0.5'].astype(str);z=np.where(lab=='0','H0 contractile candidate',np.where(lab=='2','H2 inflammatory candidate','Other candidate clusters'));sp={'Other candidate clusters':'#d7dbe0','H0 contractile candidate':'#DE8650','H2 inflammatory candidate':'#327DA8'}
corrected=pd.read_csv(p/'remediation_20260917/state_final_sensitivity/corrected_recluster_cells.tsv.gz',sep='\t',dtype={'corrected_state':str}).set_index('cell_id').reindex(a.obs_names);cl=corrected.corrected_state.to_numpy();csp={'0':'#DE8650','1':'#327DA8','2':'#6BA484','3':'#8F79A9','4':'#C6A255','5':'#5FA6A5'}
fig,axs=plt.subplots(2,2,figsize=(12.5,10));catplot(axs[0,0],old.obsm['X_umap'],patients,pp,'A  Before Harmony | patients',5);catplot(axs[0,1],xy,patients,pp,'B  After Harmony | patients',5);catplot(axs[1,0],xy,z,sp,'C  Frozen candidate states',5);catplot(axs[1,1],xy,cl,csp,'D  Corrected labels on fixed Harmony map',5)
for ax in axs.flat:ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=3,frameon=False,fontsize=8,markerscale=2)
fig.suptitle('Strict CAF candidate set | 2,047 cells from 6 patients',fontsize=15);fig.subplots_adjust(hspace=.55,bottom=.1,top=.9,wspace=.3);save(fig,'Fig05_CAF_Harmony_states_UMAP','2,047个严格CAF候选的UMAP。A为未整合坐标；B为Harmony强度2后的坐标，均按患者着色。C突出冻结H0/H2候选状态；D把去污染后重聚类标签投射到同一Harmony坐标，便于核对边界变化，D不是去污染后重新计算的UMAP。候选簇不等于已确认稳定亚型。')
