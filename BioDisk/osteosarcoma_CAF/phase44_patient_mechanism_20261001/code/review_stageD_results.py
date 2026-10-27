from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');D=P/'phase44_patient_mechanism_20261001'
assert (D/'CORE_COMPLETE').read_text().strip()=='SUCCESS' and (D/'qa_resume.exit_code').read_text().strip()=='0'
g=pd.read_csv(D/'tables/all_nominal_gene_patient_stability.tsv.gz',sep='\t');tf=pd.read_csv(D/'tables/all_TF_signed_patient_statistics.tsv.gz',sep='\t')
cols=['layer','contrast','symbol','N','mean_delta','positive_patients','P','FDR','patient_leaveout_same_direction_fraction','patient_leaveout_mean_min','largest_absolute_patient_effect_fraction']
v=g[g.dataset.eq('SCPCP000017')&g.contrast.isin(['F3_vs_rest','F4_vs_rest'])&g.symbol.isin(['COL1A1','LOX','MMP14','CCN2','LDHA','CA9','RUNX2','ANGPT1','SMAD3'])];print('MAIN_GENE_STABILITY',v[cols].to_string(index=False));v.to_csv(D/'delivery/main_selected_gene_patient_stability.tsv',sep='\t',index=False)
x=tf[tf.dataset.eq('SCPCP000017')&tf.contrast.isin(['F3_vs_rest','F4_vs_rest'])]
print('SELECTED_TF',x[x.regulator.isin(['SMAD3','TEAD1','STAT3','HIF1A','RUNX2','JUN','FOS','RELA'])][['layer','contrast','regulator','N','defined_targets','measured_targets','mean_delta','positive_patients','P','FDR']].to_string(index=False))
nom=x[x.P.lt(.05)].sort_values(['layer','contrast','P','mean_delta'],ascending=[True,True,True,False]);nom.to_csv(D/'delivery/main_all_nominal_TF_response.tsv',sep='\t',index=False)
print('TOP_TF',nom.groupby(['layer','contrast']).head(6)[['layer','contrast','regulator','N','mean_delta','positive_patients','P','FDR']].to_string(index=False))
for root in [D,D/'clinical_reporting_extension']:
 cl=pd.read_csv(root/'tables/all_clinical_fraction_statistics.tsv',sep='\t');adj=pd.read_csv(root/'tables/adjusted_clinical_fraction_statistics.tsv',sep='\t');print('CLINICAL_COUNTS',root.name,'unadjusted_valid',cl.P.notna().sum(),'unadjusted_nominal',(cl.P<.05).sum(),'adjusted_valid',adj.P.notna().sum(),'adjusted_nominal',(adj.P<.05).sum());print('ADJUSTED_NOMINAL',adj[adj.P.lt(.05)].to_string(index=False))
q=pd.read_csv(D/'tables/independent_stageD_FDR_QA.tsv',sep='\t');print('INDEPENDENT_R_FDR',len(q),'families',int(q.rows.sum()),'hypotheses','maxerr',q.max_error.max())
print('CURRENT_LR',(P/'phase45_LR_expression_20261001/STATUS.json').read_text());print('CURRENT_C2',(P/'phase43_full_transcriptome_discovery_20261001/STATUS_C2.json').read_text())