
from pathlib import Path
import json,gzip,datetime,traceback
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy.sparse import csc_matrix
O=Path(__file__).resolve().parents[1];(O/'objects/GSE293065').mkdir(parents=True,exist_ok=True)
m=pd.read_csv(O/'tables/GSE293065_sample_patient_manifest.tsv',sep='\t').fillna('').set_index('sample')
manifest=pd.read_csv(O/'tables/matrix_download_manifest.tsv',sep='\t')
outcomes=json.loads((O/'delivery/MATRIX_DOWNLOAD_OUTCOMES.json').read_text());good={x['path'] for x in outcomes if x['status']=='complete'}
ready=[];failures=[]
for sample,r in m.iterrows():
 try:
  files=manifest[manifest['sample'].eq(sample)].set_index('role');assert len(files)==4
  assert set(files.path)<=good,'At least one source asset failed download/integrity checks'
  with h5py.File(files.loc['matrix','path'],'r') as f:
   z=f['matrix'];shape=tuple(z['shape'][:]);bc=z['barcodes'].asstr()[:];ids=z['features/id'].asstr()[:];symbols=z['features/name'].asstr()[:]
   counts=csc_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=shape).T.tocsr()
  assert len(set(ids))==len(ids) and len(set(bc))==len(bc)
  pos=pd.read_csv(files.loc['spatial_coordinate','path'],header=None)
  if str(pos.iloc[0,0]).lower() in ['barcode','barcodes']:pos=pos.iloc[1:].copy()
  assert pos.shape[1]==6;pos.columns=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres'];pos=pos.set_index('barcode');assert pos.index.is_unique
  obs=pos.reindex(bc);assert not obs.isna().any().any()
  for col in obs.columns:obs[col]=pd.to_numeric(obs[col])
  for k,v in r.items():obs[k]=str(v)
  obs['sample']=sample;obs['source_barcode']=bc;obs.index=pd.Index([sample+'__'+x for x in bc],name='spot_id')
  a=ad.AnnData(X=counts,obs=obs,var=pd.DataFrame({'gene_symbol':symbols},index=pd.Index(ids,name='gene_id')))
  a.obsm['spatial']=obs[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()
  with gzip.open(files.loc['spatial_scale','path'],'rt') as f:a.uns['spatial_scalefactors']=json.load(f)
  a.uns['provenance']={'dataset':'GSE293065','matrix_path':files.loc['matrix','path'],'lowres_image_path':files.loc['low_resolution_image','path'],'input':'author provided filtered_feature_bc_matrix; no additional filtering','modality':'Visium FFPE; mixed-cell spots','analysis_status':'prepared_only_no_hypothesis_test'}
  dest=O/'objects/GSE293065'/f'{sample}_counts_spatial.h5ad';a.write_h5ad(dest,compression='gzip')
  check=ad.read_h5ad(dest,backed='r');assert check.shape==a.shape and list(check.obs_names)==list(a.obs_names);check.file.close()
  ready.append(dict(sample=sample,patient_id=r['patient_id'],gsm=r['gsm'],site_group=r['site_group'],spots=a.n_obs,features=a.n_vars,object_path=str(dest),coordinate_alignment_exact=True,no_extra_filter=True,candidate_genes=';'.join(x for x in ['PDGFRB','MCAM','GPC3','THY1','FAP','JAG1','CSPG4','GPRC5A'] if x in set(symbols))))
 except Exception:failures.append(dict(sample=sample,error=traceback.format_exc()))
pd.DataFrame(ready).to_csv(O/'tables/GSE293065_ready_objects.tsv',sep='\t',index=False)
qa={'at':datetime.datetime.now().astimezone().isoformat(),'prepared_samples':len(ready),'expected_samples':26,'patients_in_ready_objects':len(set(x['patient_id'] for x in ready)),'spots_in_ready_objects':sum(x['spots'] for x in ready),'failed_samples':failures,'validation_tests_performed':False,'spot_counts_not_biological_N':True,'no_extra_QC_gate':True}
(O/'delivery/SPATIAL_PREPARATION_QA.json').write_text(json.dumps(qa,indent=2));(O/'preparation.exit_code').write_text('0' if len(ready)==26 else '1');(O/'PREPARATION_FINISHED').write_text(qa['at']);print(json.dumps(qa))
