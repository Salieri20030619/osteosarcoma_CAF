from pathlib import Path
import sys,json,gzip,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927');gsm=sys.argv[1]
manifest=json.loads((O/'tables/GSE237070_sra_download_manifest.json').read_text());rows=[x for x in manifest if x['gsm']==gsm];assert rows
d=O/'objects/GSE237070_reconstructed'/gsm;f=d/'STAR/Solo.out/Gene/filtered'
def path(name):
 p=f/name
 if p.exists():return p
 p=f/(name+'.gz');assert p.exists();return p
x=mmread(str(path('matrix.mtx'))).tocsr();features=pd.read_csv(path('features.tsv'),sep='\t',header=None);barcodes=pd.read_csv(path('barcodes.tsv'),sep='\t',header=None)[0].astype(str)
assert x.shape==(len(features),len(barcodes));assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all();assert barcodes.is_unique
obs=pd.DataFrame({'sample':gsm,'dataset':'GSE237070','patient':'GSE237070_patient_'+rows[0]['patient'],'tissue':rows[0]['tissue'],'cell_call_source':'reconstructed_STARsolo_EmptyDrops_CR','project_QC_filter_applied':False},index=gsm+'_'+barcodes.to_numpy())
obs['n_counts']=np.asarray(x.sum(0)).ravel();obs['n_genes']=x.getnnz(0)
var=pd.DataFrame({'symbol':features[1].astype(str).to_numpy()},index=features[0].astype(str).to_numpy());a=ad.AnnData(x.T.tocsr(),obs=obs,var=var);a.write_h5ad(d/'full_native_called_counts.h5ad',compression='gzip')
panel=json.loads((O/'tables/GSE237070_VALIDATION_SCOPE.json').read_text())['markers']
coverage=[]
for gene in panel:
 ix=np.where(a.var.symbol.to_numpy()==gene)[0];coverage.append({'gene':gene,'feature_count':len(ix),'total_UMI':int(a.X[:,ix].sum()) if len(ix) else None})
result={'state':'PASS','sample':gsm,'patient':rows[0]['patient'],'tissue':rows[0]['tissue'],'source_runs':[r['run'] for r in rows],'cells':a.n_obs,'features':a.n_vars,'total_UMI':int(a.X.sum()),'source_MTX_integer_nonnegative_finite':True,'all_native_called_cells_retained':True,'additional_project_QC_thresholds':False,'markers':coverage,'CAF_biological_validation':'pending annotation and patient-level comparison; this is matrix integrity QA'}
(d/'MATRIX_QA.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
