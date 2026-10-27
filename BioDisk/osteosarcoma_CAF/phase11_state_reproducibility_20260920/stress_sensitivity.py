"""Targeted sensitivity triggered by Phase11 patient-consistent stress signal.
SYNTAX_CHECK: python -m py_compile stress_sensitivity.py
RUN_COMMAND: python -u stress_sensitivity.py
EXPECTED_SUCCESS_MARKER: STRESS_SENSITIVITY_COMPLETE
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,hashlib
import numpy as np,pandas as pd,scanpy as sc,harmonypy as hm
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'); B=P/'phase10_batch_gates_20260919'; O=P/'phase11_state_reproducibility_20260920'
R=P/'phase7_review_states_20260917/remediation_20260917/resources'
st=pd.read_csv(R/'stress_core_author.csv');cc=pd.read_csv(R/'seurat_cc_updated2019.tsv',sep='\t')
assert len(st)==512 and len(cc)==97
a=sc.read_h5ad(B/'objects/unified_strict_integrated.h5ad');SEED=20260919
stress=a.var_names.isin(st.ensembl_gene_id);cycle=a.var.symbol.isin(cc.gene).to_numpy()
protocol={'trigger':'U0 Stress_core is directionally elevated in four eligible non-GSE250015 donor records; targeted diagnostic only',
 'stress_requested':512,'stress_present_shared':int(stress.sum()),'cycle_requested':97,'cycle_present_shared':int(cycle.sum()),
 'seed':SEED,'fixed':'same 3976 cells, log1p(CP10k), HVG2000 PC30 Harmony study_technology theta2 dynamic lambda, k15 Leiden0.5',
 'changed':'exclude stress, cycle, or their union before reselecting HVGs; no score regression or cell removal',
 'interpretation':'source of boundary sensitivity, not proof stress is technical or biological',
 'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [R/'stress_core_author.csv',R/'seurat_cc_updated2019.tsv']}}
(O/'resources/stress_sensitivity_protocol.json').write_text(json.dumps(protocol,indent=2))
rows=[]
for mode,exclude in [('exclude_stress',stress),('exclude_cycle',cycle),('exclude_stress_cycle',stress|cycle)]:
    mask=~exclude&~a.var.symbol.str.startswith(('MT-','RPL','RPS')).to_numpy();b=a[:,mask].copy();sc.pp.filter_genes(b,min_cells=3)
    sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor='seurat',batch_key='study_technology');b=b[:,b.var.highly_variable].copy()
    assert not b.var_names.isin(a.var_names[exclude]).any()
    (O/'resources'/(mode+'_HVG.json')).write_text(json.dumps(b.var_names.tolist()))
    sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=SEED)
    h=hm.run_harmony(b.obsm['X_pca'],b.obs,'study_technology',theta=2,lamb=None,random_state=SEED,ncores=4,max_iter_harmony=30,verbose=False)
    z=np.asarray(h.Z_corr);z=z if z.shape==b.obsm['X_pca'].shape else z.T
    assert np.isfinite(z).all();b.obsm['diagnostic']=z
    sc.pp.neighbors(b,use_rep='diagnostic',n_neighbors=15,random_state=SEED);sc.tl.leiden(b,resolution=.5,key_added=mode,random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
    b.obs[['dataset','sample','patient_id','donor_record','donor_known','technology','site',mode]].to_csv(O/'tables'/(mode+'_labels.tsv.gz'),sep='\t',index_label='cell_id')
    rows.append(dict(model=mode,cells=b.n_obs,excluded_genes=int(exclude.sum()),HVG=b.n_vars,clusters=b.obs[mode].nunique(),ARI=adjusted_rand_score(a.obs.unified_cluster.astype(str),b.obs[mode].astype(str))))
    print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(O/'tables/stress_sensitivity_metrics.tsv',sep='\t',index=False)
(O/'STRESS_SENSITIVITY_COMPLETE').write_text('STRESS_SENSITIVITY_COMPLETE\n')
