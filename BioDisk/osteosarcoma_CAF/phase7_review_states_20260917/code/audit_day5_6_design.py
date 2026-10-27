"""SYNTAX_CHECK: python -m py_compile code/audit_day5_6_design.py
RUN_COMMAND: python code/audit_day5_6_design.py
EXPECTED_SUCCESS_MARKER: DAY5_6_DESIGN_AUDIT_COMPLETE
"""
from pathlib import Path
import json
import pandas as pd,numpy as np,anndata as ad
from scipy.stats import spearmanr
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917');out=p/'analysis'
meta=pd.read_csv(out/'patient_metadata.tsv',sep='\t').fillna('not_reported');assert meta.patient_id.is_unique
meta['site_group']=np.where(meta.site=='Lung','Lung','Bone');meta['treatment_group']=np.where(meta.treatment_raw.str.contains('GT',regex=False),'GT_based','other_recorded_regimen')
for factor in ['treatment_group','pathological_type_raw','disease_stage_raw']:pd.crosstab(meta.site_group,meta[factor]).to_csv(out/('site_by_'+factor+'.tsv'),sep='\t')
models=[]
for factors in [['site_group'],['site_group','treatment_group'],['site_group','treatment_group','pathological_type_raw']]:
 X=pd.get_dummies(meta[factors],drop_first=True,dtype=float);X.insert(0,'intercept',1.);rank=int(np.linalg.matrix_rank(X.to_numpy()));groups=meta.site_group.value_counts()
 models.append({'factors':'+'.join(factors),'patients':len(meta),'columns':X.shape[1],'rank':rank,'full_rank':rank==X.shape[1],'lung_patients':int(groups.get('Lung',0)),'bone_patients':int(groups.get('Bone',0)),'clinical_inference_allowed':False,'reason':'fewer than two lung patients; matrix rank cannot create biological replication' if groups.get('Lung',0)<2 else 'exploratory cohort; no inferential clinical model fitted in this analysis'})
pd.DataFrame(models).to_csv(out/'clinical_design_identifiability.tsv',sep='\t',index=False)
a=ad.read_h5ad(out/'GSE152048_strict_CAF_states.h5ad');assert a.obs.patient_id.nunique()==len(meta)
vars=['n_counts','n_genes','pct_mt','program_Stress_core','program_Cell_cycle'];rows=[]
for j in range(min(5,a.obsm['X_pca'].shape[1])):
 for variable in vars:
  value=float(spearmanr(a.obsm['X_pca'][:,j],a.obs[variable]).statistic)
  rows.append({'component':'PC'+str(j+1),'variable':variable,'spearman_rho':value,'unit':'cell-level descriptive association','P_value':'not_computed_no_cell_pseudoreplication'})
pd.DataFrame(rows).to_csv(out/'PCA_quality_stress_cycle_associations.tsv',sep='\t',index=False)
assert int(pd.read_csv(out/'patient_by_state_cells.tsv',sep='\t',index_col=0).to_numpy().sum())==a.n_obs
ext=pd.read_csv(out/'extended_patient_by_state_cells.tsv',sep='\t',index_col=0);summary=json.loads((out/'analysis_summary.json').read_text());assert int(ext.to_numpy().sum())==summary['extended_cells']
# Test-family membership and missing-test behavior must remain explicit.
genes=pd.read_csv(out/'state_gene_effects_all.tsv.gz',sep='\t');valid=genes.test_status.str.startswith('exploratory')
assert genes.loc[valid,'P_value'].between(0,1).all() and genes.loc[valid,'FDR_BH'].between(0,1).all();assert genes.loc[~valid,'P_value'].isna().all()
meta.to_csv(out/'patient_design_metadata.tsv',sep='\t',index=False)
(out/'DESIGN_AUDIT_COMPLETE').write_text('DAY5_6_DESIGN_AUDIT_COMPLETE\n');print('DAY5_6_DESIGN_AUDIT_COMPLETE',a.n_obs,len(meta),flush=True)
