from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
s=pd.read_csv(O/'tables/program_response_statistics.tsv',sep='\t');g=pd.read_csv(O/'tables/fixed_program_gene_response_statistics.tsv',sep='\t')
print('ADJUSTED_NOMINAL',s[(s.model=='adjusted_HC3')&(s.p<.05)].to_string(index=False))
print('GENE_COUNTS',g.groupby('cohort').agg(n=('gene','size'),nominal=('p',lambda x:(x<.05).sum()),fdr=('fdr_gene_family',lambda x:(x<.05).sum())).to_string())
w=g.pivot(index='gene',columns='cohort',values=['p','effect','fdr_gene_family'])
both=w[(w[('p','GSE39055')]<.05)&(w[('p','GSE21257')]<.05)]
print('BOTH_NOMINAL',both.to_string())
a=g[g.cohort=='GSE39055'].set_index('gene');b=g[g.cohort=='GSE21257'].set_index('gene');common=a.join(b,lsuffix='_39055',rsuffix='_21257',how='outer')
common['same_direction']=common.effect_39055*common.effect_21257>0
common['both_nominal']= (common.p_39055<.05)&(common.p_21257<.05)
common.reset_index().to_csv(O/'tables/cross_cohort_gene_evidence.tsv',sep='\t',index=False)
print('COMPLETE')
