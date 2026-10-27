from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase51_matrix_lineage_attribution_20261002')
s=pd.read_csv(R/'tables/lineage_statistics.tsv',sep='\t')
summary=pd.read_csv(R/'tables/lineage_source_summary.tsv',sep='\t')
receivers=['Fibroblast_other','Mural','Endothelial','Myeloid','T_NK','CNV_altered_unresolved','Other_nonfrozen_unresolved','Annotated_osteogenic_nonfrozen','CNV_altered_osteogenic_annotated']
labels=['Other fibroblast annotation','Mural annotation','Endothelial annotation','Myeloid annotation','T / NK annotation','CNV altered, unresolved','Other unresolved','Osteogenic annotation','Osteogenic + CNV altered']
genes=['COL1A1','LOX','MMP14']
fig,axs=plt.subplots(1,2,figsize=(11.8,7.5),sharey=True)
for ax,target,title in zip(axs,['supported__F3','expanded__F4'],['Supported F3','Expanded F4']):
 t=s[s.dataset.eq('SCPCP000017')&s.target.eq(target)]
 w=t.pivot(index='receiver',columns='gene',values='mean_effect').reindex(index=receivers,columns=genes)
 cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#eeeeee')
 im=ax.imshow(w,cmap=cm,norm=matplotlib.colors.TwoSlopeNorm(vmin=-3,vcenter=0,vmax=7),aspect='auto')
 for j,receiver in enumerate(receivers):
  for k,gene in enumerate(genes):
   r=t[t.receiver.eq(receiver)&t.gene.eq(gene)].iloc[0]
   if r.N==0:txt='Not evaluable\nN = 0'
   else:txt=f'{r.mean_effect:+.2f}'+('*' if r.FDR<.05 else '')+f'\nN={int(r.N)} | P={r.P:.3g}' if np.isfinite(r.P) else f'{r.mean_effect:+.2f}\nN={int(r.N)} | no test'
   ax.text(k,j,txt,ha='center',va='center',fontsize=8,color='white' if r.mean_effect>5.2 else 'black')
 ax.set_xticks(range(3),genes);ax.set_yticks(range(9),labels);ax.set_title(title)
fig.subplots_adjust(left=.27,right=.89,top=.87,bottom=.1,wspace=.12);cax=fig.add_axes([.92,.20,.018,.55]);fig.colorbar(im,cax=cax,label='Mean log2 normalized-expression difference')
fig.suptitle('SCPCP000017 | target minus comparator',y=.96)
fig.text(.27,.035,'* BH FDR < 0.05 within 27 hypotheses per target; exact patient sign flip',fontsize=9)
fig.savefig(R/'figures/main_matched_lineage_effects.png',dpi=170);fig.savefig(R/'figures/main_matched_lineage_effects.pdf');plt.close(fig)
groups=['CAF_frozen_all','Fibroblast_other','Mural','Endothelial','Myeloid','T_NK','CNV_altered_unresolved','Other_nonfrozen_unresolved']
legend=['Frozen CAF candidates','Other fibroblast annotation','Mural','Endothelial','Myeloid','T / NK','CNV altered, unresolved','Other unresolved']
colors=['#237f78','#67b9a9','#dca653','#6d99c2','#bc6a71','#947ca8','#666666','#c6c6c6']
fig,ax=plt.subplots(figsize=(10,5.8));bottom=np.zeros(3)
for group,label,color in zip(groups,legend,colors):
 d=summary[summary.dataset.eq('SCPCP000017')&summary.group.eq(group)].set_index('gene').reindex(genes)
 y=d.mean_gene_capture_share.to_numpy()*100
 ax.bar(genes,y,bottom=bottom,color=color,label=label,width=.65)
 for j,v in enumerate(y):
  if v>=4:ax.text(j,bottom[j]+v/2,f'{v:.1f}%',ha='center',va='center',fontsize=10)
 bottom+=y
assert np.allclose(bottom,100)
ax.set_ylabel('Mean captured transcript share (%)');ax.set_ylim(0,100);ax.set_title('SCPCP000017 | 20 source-scoped patients')
ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),frameon=False)
fig.subplots_adjust(left=.09,right=.64,bottom=.15,top=.9)
fig.text(.09,.035,'All called cells; equal site/patient means. Captured RNA, not spatial-cell attribution.',fontsize=9)
fig.savefig(R/'figures/main_RNA_source_shares.png',dpi=170);fig.savefig(R/'figures/main_RNA_source_shares.pdf');plt.close(fig)
# Source-specific contrasts: no pooled inference or independence claim across datasets.
t=s[s.target.eq('supported__F3')&s.receiver.eq('Fibroblast_other')]
datasets=sorted(t.dataset.unique());fig,axs=plt.subplots(1,3,figsize=(12,5),sharey=True)
for ax,gene in zip(axs,genes):
 d=t[t.gene.eq(gene)].set_index('dataset').reindex(datasets)
 for j,(ds,r) in enumerate(d.iterrows()):
  if np.isfinite(r.mean_effect):
   ax.scatter(r.mean_effect,j,color='#237f78' if r.P<.05 else '#777777',s=45)
   ax.annotate(f'N={int(r.N)}',(r.mean_effect,j),xytext=(5,4),textcoords='offset points',fontsize=8)
  else:ax.text(0,j,'N=0',va='center',fontsize=8,color='#aaaaaa')
 ax.axvline(0,color='#bbbbbb',lw=.8);ax.set_title(gene);ax.set_yticks(range(len(datasets)),datasets);ax.set_xlabel('Mean log2 expression difference');ax.grid(axis='x',alpha=.2);ax.margins(x=.3)
axs[0].set_ylim(len(datasets)-.5,-.5);fig.suptitle('Supported F3 vs other fibroblast annotation | sources kept separate')
fig.tight_layout(rect=[0,.065,1,.94]);fig.text(.12,.015,'Green: raw P < 0.05; gray: non-nominal or descriptive; N=0: not evaluable',fontsize=9);fig.savefig(R/'figures/source_specific_F3_effects.png',dpi=170);fig.savefig(R/'figures/source_specific_F3_effects.pdf');plt.close(fig)
(R/'delivery/FIGURE_RENDER_COMPLETE').write_text('SUCCESS')
