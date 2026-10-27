from pathlib import Path
import anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
files=[P/'phase3_identity_20260914/identity/objects/caf_strict.h5ad',P/'phase7_review_states_20260917/identity/objects/GSE152048_caf_strict.h5ad']+sorted((P/'geo_identity_20260919/objects').glob('*.h5ad'))
shared=None
for p in files:
 a=ad.read_h5ad(p,backed='r');symbols=a.var.symbol.astype(str)
 shared=set(symbols) if shared is None else shared&set(symbols)
 fields=[x for x in ['sample','sample_id','short','patient_id','site','dataset','technology','technology_reported','patient_mapping_status','patient_independence_status'] if x in a.obs]
 print(json.dumps({'file':str(p),'shape':a.shape,'layers':list(a.layers),'gene_examples':a.var_names[:3].tolist(),'symbol_duplicates':int(symbols.duplicated().sum()),'mapping':a.obs[fields].drop_duplicates().to_dict('records')},ensure_ascii=False))
 a.file.close()
print('COMMON_SYMBOLS',len(shared))
