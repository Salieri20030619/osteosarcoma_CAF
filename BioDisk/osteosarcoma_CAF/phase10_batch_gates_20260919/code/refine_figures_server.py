from pathlib import Path
import pandas as pd,numpy as np,matplotlib,json,hashlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase10_batch_gates_20260919')
plt.rcParams.update({'font.size':10})
m=pd.read_csv(O/'tables/batch_model_comparison.tsv',sep='\t');mm=m[~m.analysis.str.startswith('leaveout_')].copy();y=np.arange(len(mm));fig,axs=plt.subplots(1,2,figsize=(15,8))
axs[0].barh(y,mm.same_batch_neighbors,color='#3a7ca5');axs[0].scatter(mm.composition_expected_same_batch,y,color='#c45536',marker='|',s=90,label='Composition expectation');axs[0].set_yticks(y,mm.analysis);axs[0].set_ylim(-.6,len(mm)-.4);axs[0].invert_yaxis();axs[0].set_xlabel('Same-batch neighbor fraction');axs[0].legend(loc='upper center',bbox_to_anchor=(.5,-.1),frameon=False);axs[0].set_title('Batch mixing, interpreted with source composition')
v=mm.ARI_vs_unified_primary.to_numpy().copy();v[mm.analysis.eq('unified_primary').to_numpy()]=1.;axs[1].barh(y,v,color='#577d5d');axs[1].set_yticks(y,[str(s)+' (self)' if s=='unified_primary' else s for s in mm.analysis]);axs[1].set_ylim(-.6,len(mm)-.4);axs[1].invert_yaxis();axs[1].set_xlim(0,1.02)
for i,z in enumerate(v):
 if pd.isna(z):axs[1].text(.025,i,'not computed',va='center',color='#666666',fontsize=8)
axs[1].set_title('Agreement on shared cells');axs[1].set_xlabel('ARI vs unified primary');fig.tight_layout();fig.savefig(O/'figures/03_batch_sensitivity.png',dpi=170);plt.close(fig)
s=pd.read_csv(O/'tables/primary_cluster_stability.tsv',sep='\t');labels=pd.read_csv(O/'tables/unified_primary_labels.tsv.gz',sep='\t');cs=sorted(labels.unified_primary.unique());ct=pd.crosstab(labels.unified_primary,labels.dataset).reindex(cs);ct=ct.div(ct.sum(1),axis=0)
fig,axs=plt.subplots(1,2,figsize=(15,7.5));ct.plot.bar(stacked=True,ax=axs[0],colormap='tab10');axs[0].set_title('Study composition of unified clusters');axs[0].set_ylabel('Fraction of cells');axs[0].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=3,frameon=False)
s.set_index('cluster')[['LOO_min_Jaccard','LOO_median_Jaccard','HVG_PC_min_Jaccard']].plot.bar(ax=axs[1]);axs[1].set_ylim(0,1.03);axs[1].set_title('LOO overlap conditional on >=20 remaining cells');axs[1].set_ylabel('Best-match Jaccard');axs[1].set_xticklabels([f'U{int(r.cluster)}\n{int(r.LOO_evaluable_runs)}/14' for r in s.itertuples()],rotation=0,fontsize=8);axs[1].set_xlabel('Cluster / evaluable leave-out runs');axs[1].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.18),ncol=2,frameon=False);fig.subplots_adjust(bottom=.26,top=.9,wspace=.27);fig.savefig(O/'figures/04_cluster_support.png',dpi=170);plt.close(fig)
m=json.loads((O/'MANIFEST_SHA256.json').read_text(encoding='utf-8'))
for n in ['figures/03_batch_sensitivity.png','figures/04_cluster_support.png']:m[n]=hashlib.sha256((O/n).read_bytes()).hexdigest()
(O/'MANIFEST_SHA256.json').write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')

print('REMOTE_FIGURES_REFINED')
