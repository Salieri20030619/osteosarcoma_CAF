from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase52_blind_lineage_validation_20261002')
b=pd.read_csv(R/'tables/control_agreement_summary.tsv',sep='\t').sort_values('annotation')
fig,ax=plt.subplots(figsize=(10,6.5));y=np.arange(len(b))
ax.barh(y-.18,b.mean_full_agreement*100,height=.34,label='Full-panel agreement',color='#467994')
ax.barh(y+.18,b.mean_correct_consensus_fraction*100,height=.34,label='All controls: agreeing consensus',color='#6eae99')
ax.set_yticks(y,[f'{r.annotation} (N={r.N_patients})' for r in b.itertuples()]);ax.invert_yaxis();ax.set_xlim(0,110)
for j,r in enumerate(b.itertuples()):
 ax.text(r.mean_full_agreement*100+.8,j-.18,f'{100*r.mean_full_agreement:.1f}%',va='center',fontsize=8)
 ax.text(r.mean_correct_consensus_fraction*100+.8,j+.18,f'{100*r.mean_correct_consensus_fraction:.1f}%',va='center',fontsize=8)
ax.set_xlabel('Patient-mean agreement with existing annotations (%)');ax.set_title('Target-gene-heldout panel benchmark');ax.legend(loc='lower right',fontsize=9,frameon=False)
fig.tight_layout(rect=[0,.055,1,1]);fig.text(.07,.015,'Annotation agreement is not external accuracy. Osteogenic controls: none.',fontsize=9)
fig.savefig(R/'figures/control_annotation_agreement.png',dpi=170);fig.savefig(R/'figures/control_annotation_agreement.pdf');plt.close(fig)
s=pd.read_csv(R/'tables/query_prediction_share_summary.tsv',sep='\t')
s=s[s.scope.eq('query_primary_other')]
cats=['Fibroblast','Osteogenic','Other','Panel_discordant','Unresolved'];labels=['Fibroblast-like','Osteogenic-like','Other panel labels','Panel-discordant','Unresolved']
colors=['#268278','#dcad60','#7b99bc','#ae7685','#c9c9c9'];genes=['COL1A1','LOX','MMP14']
fig,axs=plt.subplots(1,2,figsize=(11,6),sharey=True)
for ax,scheme,title in zip(axs,['full_panel','split_panel_consensus'],['Full panel','Disjoint-panel consensus']):
 d=s[s.scheme.eq(scheme)].copy();d['category']=d.predicted_label.where(d.predicted_label.isin(cats),'Other')
 w=d.groupby(['category','gene']).mean_within_query_gene_share.sum().unstack().reindex(index=cats,columns=genes).fillna(0)
 assert np.allclose(w.sum(0),1)
 bottom=np.zeros(3)
 for category,label,color in zip(cats,labels,colors):
  v=w.loc[category].to_numpy()*100;ax.bar(genes,v,bottom=bottom,color=color,label=label,width=.66)
  for j,a in enumerate(v):
   if a>=5:ax.text(j,bottom[j]+a/2,f'{a:.1f}%',ha='center',va='center',fontsize=9)
  bottom+=v
 ax.set_ylim(0,100);ax.set_title(title)
axs[0].set_ylabel('Mean share of RNA within primary unassigned cells (%)')
handles,labs=axs[1].get_legend_handles_labels();fig.legend(handles,labs,loc='lower center',bbox_to_anchor=(.5,.07),ncol=3,frameon=False,fontsize=9)
fig.suptitle('Primary unassigned subset | 93,283 cells | 20 patients')
fig.subplots_adjust(left=.09,right=.98,top=.85,bottom=.24,wspace=.12)
fig.text(.08,.025,'Patient/site equal weighting; computational labels do not establish CAF or malignant identity.',fontsize=9)
fig.savefig(R/'figures/query_RNA_prediction_sensitivity.png',dpi=170);fig.savefig(R/'figures/query_RNA_prediction_sensitivity.pdf');plt.close(fig)
st=pd.read_csv(R/'tables/target_gene_statistics.tsv',sep='\t')
contrasts=['query_Fibroblast_minus_query_Osteogenic','query_Fibroblast_minus_supported_F3','query_Osteogenic_minus_supported_F3']
clabels=['Fibroblast-like minus osteogenic-like','Fibroblast-like minus supported F3','Osteogenic-like minus supported F3']
fig,axs=plt.subplots(1,2,figsize=(12,5.7),sharey=True)
for ax,scheme,title in zip(axs,['full_panel','split_panel_consensus'],['Full panel','Disjoint-panel consensus']):
 t=st[st.scheme.eq(scheme)];w=t.pivot(index='contrast',columns='gene',values='mean_effect').reindex(index=contrasts,columns=genes)
 im=ax.imshow(w,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-4.5,vcenter=0,vmax=1))
 for j,k in enumerate(contrasts):
  for i,g in enumerate(genes):
   r=t[t.contrast.eq(k)&t.gene.eq(g)].iloc[0]
   ax.text(i,j,f'{r.mean_effect:+.2f}\nN={int(r.N)}; P={r.P:.3g}\nq={r.FDR:.3g}',ha='center',va='center',fontsize=8,color='white' if r.mean_effect<-3 else 'black')
 ax.set_xticks(range(3),genes);ax.set_yticks(range(3),clabels);ax.set_title(title)
fig.subplots_adjust(left=.32,right=.89,top=.81,bottom=.14,wspace=.1);ca=fig.add_axes([.92,.2,.018,.56]);fig.colorbar(im,cax=ca,label='Mean log2 expression difference')
fig.suptitle('Held-out target-gene expression | patient-level comparisons')
fig.text(.31,.045,'9 hypotheses per scheme. N < 10 remains exploratory.',fontsize=9)
fig.savefig(R/'figures/heldout_gene_patient_effects.png',dpi=170);fig.savefig(R/'figures/heldout_gene_patient_effects.pdf');plt.close(fig)
(R/'delivery/FIGURE_RENDER_COMPLETE').write_text('SUCCESS')
