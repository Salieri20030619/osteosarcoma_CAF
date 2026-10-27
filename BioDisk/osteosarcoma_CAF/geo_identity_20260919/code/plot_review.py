from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_identity_20260919')
s=pd.read_csv(O/'tables/sample_identity_summary.tsv',sep='\t')
t=pd.read_csv(O/'tables/candidate_sensitivity_tiers.tsv',sep='\t')
s=s.merge(t[['sample','count_threshold_sensitive_candidate']],on='sample',validate='one_to_one')
s=s[s.FibLike>0];y=np.arange(len(s));left=np.zeros(len(s))
v=[s.supported_candidate,s.count_threshold_sensitive_candidate,s.raw_strict_candidate-s.supported_candidate-s.count_threshold_sensitive_candidate,s.dual_low_FibLike-s.raw_strict_candidate,s.FibLike-s.dual_low_FibLike]
labels=['Primary supported candidate','Detection-threshold-sensitive candidate','Other corrected-lineage unresolved','Other dual-low cells','CNV elevated / unresolved / unassessed']
colors=['#287D65','#DFA83D','#B36C35','#718DAC','#C8CDD2']
fig,ax=plt.subplots(figsize=(13,6.7))
for val,lab,col in zip(v,labels,colors):
 assert (val>=0).all()
 ax.barh(y,val,left=left,label=lab,color=col);left+=val
assert np.array_equal(left,s.FibLike)
ax.set_yticks(y,[v.split('__')[1]+(' [lung]' if site=='lung' else '') for v,site in zip(s['sample'],s.site)])
ax.invert_yaxis();ax.set_xlabel('Cells (descriptive counts)');ax.set_title('CAF identity review | all original FibLike cells retained')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2,frameon=False,fontsize=9)
fig.tight_layout();fig.savefig(O/'figures/identity_funnel.png',dpi=160,bbox_inches='tight');plt.close(fig)
print('IDENTITY_FIGURE_RENDERED')
