"""SYNTAX_CHECK: python -m py_compile code/mycaf_effects.py
RUN_COMMAND: python -u code/mycaf_effects.py
EXPECTED_SUCCESS_MARKER: MYCAF_EFFECT_REVIEW_COMPLETE
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json
from itertools import product
import numpy as np,pandas as pd,anndata as ad
from statsmodels.stats.multitest import multipletests
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917/analysis')
info=json.loads((p/'mycaf_subdivision_status.json').read_text())
if not info['executed']:
 (p/'MYCAF_EFFECT_COMPLETE').write_text('MYCAF_EFFECT_REVIEW_COMPLETE conditional subdivision not supported\n');print('MYCAF_EFFECT_REVIEW_COMPLETE',info);raise SystemExit(0)
a=ad.read_h5ad(p/'mycaf_subclusters.h5ad');x=a.layers['counts'];rows=[];support=[];programs=[c for c in a.obs if c.startswith('program_')]
for state in sorted(a.obs.mycaf_cluster.unique()):
 mask=(a.obs.mycaf_cluster==state).values;patient_counts=a.obs.loc[mask].patient_id.value_counts();deltas=[]
 for patient in a.obs.patient_id.unique():
  within=(a.obs.patient_id==patient).values;yes=mask&within;no=~mask&within
  if yes.sum()<20 or no.sum()<20:continue
  y=np.asarray(x[yes].sum(0)).ravel();n=np.asarray(x[no].sum(0)).ravel();deltas.append(np.log2(y/y.sum()*1e6+1)-np.log2(n/n.sum()*1e6+1))
 target=np.asarray(x[mask].sum(0)).ravel();rest=np.asarray(x[~mask].sum(0)).ravel();pooled=np.log2(target/target.sum()*1e6+1)-np.log2(rest/rest.sum()*1e6+1) if rest.sum()>0 else np.full(a.n_vars,np.nan)
 arr=np.asarray(deltas);m=pd.DataFrame({'mycaf_state':state,'ensembl_id':a.var_names,'symbol':a.var.symbol.values,'paired_patients':len(deltas),'median_log2CPM_difference':np.median(arr,axis=0) if len(deltas) else pooled,'target_pct_detected':np.asarray((x[mask]>0).mean(0)).ravel()*100,'patient_positive_fraction':(arr>0).mean(0) if len(deltas) else np.full(a.n_vars,np.nan)})
 if len(deltas)>=3:
  assert len(deltas)<=12;signs=np.array(list(product([-1.,1.],repeat=len(deltas))));pv=(np.abs(signs@arr)>=np.abs(arr.sum(0))[None,:]-1e-10).mean(0);m['P_value']=pv;m['FDR_BH']=multipletests(pv,method='fdr_bh')[1];m['test_status']='exploratory_patient_paired_sign_flip'
 else:m['P_value']=np.nan;m['FDR_BH']=np.nan;m['test_status']='not_tested_insufficient_paired_patients'
 m['FDR_family']='myCAF_'+state+'_vs_other_myCAF_all_genes';rows.append(m)
 support.append({'mycaf_state':state,'cells':int(mask.sum()),'patients_at_least_10_cells':int((patient_counts>=10).sum()),'max_patient_fraction':float(patient_counts.max()/patient_counts.sum()),'cross_patient_support':bool((patient_counts>=10).sum()>=3),**{c:float(a.obs.loc[mask,c].mean()) for c in programs}})
full=pd.concat(rows,ignore_index=True);full.to_csv(p/'mycaf_gene_effects_all.tsv.gz',sep='\t',index=False)
full[full.target_pct_detected>=20].sort_values(['mycaf_state','median_log2CPM_difference'],ascending=[True,False]).groupby('mycaf_state').head(30).to_csv(p/'mycaf_top30_effect_markers.tsv',sep='\t',index=False)
pd.DataFrame(support).to_csv(p/'mycaf_subcluster_support_programs.tsv',sep='\t',index=False)
(p/'MYCAF_EFFECT_COMPLETE').write_text('MYCAF_EFFECT_REVIEW_COMPLETE\n');print('MYCAF_EFFECT_REVIEW_COMPLETE',len(support),flush=True)
