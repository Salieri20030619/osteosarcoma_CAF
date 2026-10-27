"""SYNTAX_CHECK: python -m py_compile code/validate_stage.py
RUN_COMMAND: python code/validate_stage.py
EXPECTED_SUCCESS_MARKER: PHASE5_OBJECT_AUDIT_PASS
"""
from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np,anndata as ad
from PIL import Image
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase5_nonpdx_qc_20260916')
manifest=pd.read_csv(r/'manifest.tsv',sep='\t').fillna('');summary=pd.read_csv(r/'results/sample_summary.tsv',sep='\t')
assert len(summary)==len(manifest)==17
assert set(summary['sample'])==set(manifest['sample'])
assert (r/'results/SUCCESS.txt').read_text().strip()=='PHASE5_QC_CLUSTERING_COMPLETE'
audit=[]
for row in manifest.to_dict('records'):
 p=r/'results'/row['sample'];a=ad.read_h5ad(p/'annotated.h5ad');x=a.layers['counts'];qc=pd.read_csv(p/'all_cells_qc.tsv.gz',sep='\t');s=json.loads((p/'complete.json').read_text())
 assert len(qc)==int(row['input_cells'])==s['input_cells']
 assert int(qc.analysis_keep.sum())==a.n_obs==s['retained']
 assert set(qc.loc[qc.analysis_keep,'cell_id'])==set(a.obs_names)
 assert a.obs_names.is_unique and a.var_names.is_unique and 'symbol' in a.var
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
 assert np.allclose(np.asarray(x.sum(1)).ravel(),a.obs.n_counts)
 for field in ['sample','patient_id','site','group']:
  expected=row['group'] if field=='group' else row[field]
  assert set(a.obs[field].astype(str))=={str(expected)},(row['sample'],field)
 rng=np.random.default_rng(270231);idx=rng.choice(a.n_obs,min(a.n_obs,256),replace=False)
 b=x[idx].copy();b=b.multiply((1e4/np.asarray(b.sum(1)).ravel())[:,None]).tocsr();b.data=np.log1p(b.data)
 diff=(a.X[idx]-b).tocsr();assert not len(diff.data) or np.max(np.abs(diff.data))<1e-5
 assert np.isfinite(a.obsm['X_umap']).all()
 mk=pd.read_csv(p/'cluster_markers.tsv',sep='\t');assert np.isfinite(mk[['scores','pvals','pvals_adj']].to_numpy()).all()
 image_checks={}
 for name in ['QC.png','UMAP.png','Markers.png']:
  with Image.open(p/name) as im:image_checks[name]=list(im.size);im.verify()
 audit.append({'sample':row['sample'],'cells':a.n_obs,'genes':a.n_vars,'qc_accounting':'PASS','counts':'PASS','clinical_fields':'PASS','normalization_256_cells':'PASS','marker_tables':'PASS','image_files':image_checks,'visual_layout_review':'not_performed_local_transfer_approval_rejected'})
 del a,x,qc,b,diff
(r/'review/object_audit.json').write_text(json.dumps(audit,indent=2))
(r/'review/OBJECT_AUDIT_COMPLETE').write_text('PHASE5_OBJECT_AUDIT_PASS\n')
print('PHASE5_OBJECT_AUDIT_PASS',len(audit),sum(x['cells'] for x in audit),flush=True)
