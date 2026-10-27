from pathlib import Path
import json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002');B=R/'blind_review'
def bin_sort(x):
 chrom,pos=x.split(':');c=chrom.replace('chr','');return (int(c) if c.isdigit() else 23 if c=='X' else 24 if c=='Y' else 99,int(pos))
for d in sorted(B.glob('D[0-9][0-9]')):
 source_id=d.name;anon=pd.read_csv(d/'cell_structure.tsv.gz',sep='\t');meta=pd.read_csv(d/'clusters.tsv',sep='\t')
 fig,axs=plt.subplots(1,2,figsize=(12,6.1))
 for ax,col,title in zip(axs,['r10','patient_id'],['Primary clusters','Patients']):
  levels=sorted(anon[col].unique());codes=pd.Categorical(anon[col],categories=levels).codes
  ax.scatter(anon.UMAP1,anon.UMAP2,c=codes,s=.5,cmap='turbo',rasterized=True,alpha=.65)
  ax.set_title(title);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
  if col=='r10':
   for label,t in anon.groupby(col):
    ax.text(t.UMAP1.median(),t.UMAP2.median(),label.split('_')[-1],fontsize=6,ha='center',bbox=dict(facecolor='white',alpha=.6,edgecolor='none',pad=.2))
  else:
   handles=[Line2D([0],[0],marker='o',lw=0,markersize=4,color=plt.get_cmap('turbo')(j/max(len(levels)-1,1)),label=label) for j,label in enumerate(levels)]
   fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.75,.015),ncol=min(7,len(levels)),fontsize=6,frameon=False,title='Patient IDs',title_fontsize=7)
 fig.suptitle(source_id+' | expression-space view');fig.tight_layout(rect=[0,.13,1,1])
 fig.savefig(d/'expression_space.png',dpi=150);fig.savefig(d/'expression_space.pdf');plt.close(fig)
 f=d/'copykat_patient_cluster_profiles.tsv.gz'
 if f.exists():
  patient=pd.read_csv(f,sep='\t');bins=sorted([x for x in patient if x.startswith('chr')],key=bin_sort)
  mean=patient.groupby('cluster_id')[bins].mean().reindex(meta.cluster_id)
  finite=np.abs(mean.to_numpy());finite=finite[np.isfinite(finite)]
  scale=max(float(np.quantile(finite,.95)),.05) if len(finite) else .05
  fig,ax=plt.subplots(figsize=(14,max(5,len(meta)*.14)));cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dddddd')
  im=ax.imshow(mean,aspect='auto',cmap=cm,vmin=-scale,vmax=scale,interpolation='none')
  starts=[];labels=[]
  for j,binid in enumerate(bins):
   chrom=binid.split(':')[0]
   if not labels or labels[-1]!=chrom:starts.append(j);labels.append(chrom)
  ends=starts[1:]+[len(bins)];centers=[(a+b-1)/2 for a,b in zip(starts,ends)]
  for x in starts:ax.axvline(x-.5,color='black',lw=.25)
  ax.set_xticks(centers,[x.replace('chr','') for x in labels],fontsize=7)
  ax.set_yticks(range(len(meta)),[x.split('_')[-1] for x in meta.cluster_id],fontsize=6)
  ax.set_xlabel('Chromosome | native coordinates, 5 Mb display bins');ax.set_ylabel('Anonymous primary cluster')
  ax.set_title(source_id+' | patient-mean native CopyKAT signal');fig.colorbar(im,ax=ax,label='Native signal (display scale only)',shrink=.8);fig.tight_layout()
  fig.savefig(d/'copykat_cluster_profiles.png',dpi=150);fig.savefig(d/'copykat_cluster_profiles.pdf');plt.close(fig)
 print('RENDERED',source_id,flush=True)
(R/'BLIND_FIGURES_FINAL_COMPLETE').write_text('SUCCESS')
