from pathlib import Path
import sys,os,json,hashlib,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase46_spatial_functional_20261001';O=P/'phase50_F3_spatial_decomposition_20261002'
sys.path.insert(0,str(S/'code'))
from common import feature_map,corr,test_table,bh
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
pr=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text())
stats=pd.read_csv(O/'tables/statistics.tsv',sep='\t');patients=pd.read_csv(O/'tables/patient_effects.tsv',sep='\t')
names=list(pr['signatures']);labels=['COL1A1','LOX','MMP14','LOX + MMP14','COL1A1 + MMP14','COL1A1 + LOX','Full 3-gene reference']
fig,axs=plt.subplots(1,2,figsize=(12,6),sharey=True)
for ax,ep in zip(axs,pr['endpoints']):
 for j,name in enumerate(names):
  v=patients[(patients.signature==name)&(patients.model=='full_adjusted')&(patients.endpoint==ep)].sort_values('patient_id').effect
  ax.scatter(v,np.full(len(v),j),s=16,alpha=.55,color='#527a9b')
  r=stats[(stats.signature==name)&(stats.model=='full_adjusted')&(stats.endpoint==ep)].iloc[0]
  ax.scatter([r.mean_effect],[j],s=60,marker='D',color='#ae492c',zorder=4)
 ax.axvline(0,color='gray',lw=.8);ax.set_yticks(range(7),labels);ax.set_title('Local CD8 RNA' if ep=='CD8_RNA' else 'Adjacent CD8 RNA');ax.set_xlabel('Patient mean partial rank correlation');ax.grid(axis='x',alpha=.2)
axs[0].invert_yaxis();fig.suptitle('F3 signal decomposition | 19 patients | full adjustment');fig.tight_layout();fig.savefig(O/'figures/patient_decomposition.png',dpi=170);fig.savefig(O/'figures/patient_decomposition.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(13,6),sharey=True)
for ax,ep in zip(axs,pr['endpoints']):
 w=stats[stats.endpoint.eq(ep)].pivot(index='signature',columns='model',values='mean_effect').reindex(index=names,columns=list(pr['models']))
 im=ax.imshow(w,cmap='RdBu_r',vmin=-.12,vmax=.12)
 for j,name in enumerate(names):
  for k,model in enumerate(pr['models']):
   r=stats[(stats.signature==name)&(stats.model==model)&(stats.endpoint==ep)].iloc[0]
   ax.text(k,j,f'{r.mean_effect:.3f}',ha='center',va='center',fontsize=9,color='black')
 ax.set_yticks(range(7),labels);ax.set_xticks(range(5),['Raw','Library + XY','Cancer + library + XY','Fibroblast + library + XY','Full'],rotation=40,ha='right');ax.set_title('Local CD8 RNA' if ep=='CD8_RNA' else 'Adjacent CD8 RNA')
fig.subplots_adjust(left=.17,right=.86,bottom=.25,top=.9,wspace=.14);cax=fig.add_axes([.90,.30,.018,.48]);fig.colorbar(im,cax=cax,label='Patient mean rank correlation');fig.savefig(O/'figures/model_dependence.png',dpi=170);fig.savefig(O/'figures/model_dependence.pdf');plt.close(fig)
