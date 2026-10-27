import pathlib,pandas as pd,anndata as ad,json,hashlib,numpy as np
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase99_six_patient_admission_check_20261020';a=root/'phase93_expanded_admission_20261016'
patients=pd.read_csv(a/'tables/admitted_patients.tsv',sep='\t',encoding='utf-8',dtype=str)
manifest=json.loads((a/'objects/COUNT_COLLECTION.json').read_text(encoding='utf-8'))
reg=pd.read_csv(a/'tables/admitted_expanded_cells.tsv.gz',sep='\t',encoding='utf-8')
print('registry columns',reg.columns.tolist())
mapping={lib:r.patient_id for r in patients.itertuples() for lib in r.libraries.split(';')}
genes=json.loads((root/'phase96_all14_literature_20261016/project_registry.json').read_text(encoding='utf-8'))
rows=[];libs=[]
for z in manifest['added_objects']:
 f=pathlib.Path(z['path']);assert hashlib.sha256(f.read_bytes()).hexdigest()==z['sha256']
 o=ad.read_h5ad(f);assert o.n_obs==z['cells'];print(z['sample'],o.shape,'var',o.var.columns.tolist(),'first',o.var_names[:2].tolist())
 denom=float(o.X.sum());assert denom>0
 for gene in genes:
  ids=o.var_names.astype(str).str.split('.').str[0];ix=np.flatnonzero(ids==gene['feature_ids'])
  assert len(ix)<=1
  raw=float(o.X[:,ix].sum()) if len(ix) else np.nan
  rows.append(dict(patient_id=mapping[z['sample']],sample=z['sample'],gene=gene['gene'],feature_id=gene['feature_ids'],cells=o.n_obs,library_UMIs=denom,raw_count=raw,measured=bool(len(ix)),log2_CP10k_plus_0_1=float(np.log2(raw/denom*1e4+.1)) if len(ix) else np.nan,site='Unknown',inferential_comparison_eligible=False))
 libs.append(dict(sample=z['sample'],patient_id=mapping[z['sample']],cells=o.n_obs,source='NCC_Zheng2024',sequenced_lesion_site='Unknown',site_evidence='',baseline_overlap_status='global_patient_independence_not_verified'))
pd.DataFrame(rows).to_csv(p/'tables/new_library_14gene_descriptive_expression.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
pd.DataFrame(libs).to_csv(p/'tables/sample_site_information_needed.tsv',sep='\t',encoding='utf-8',index=False)
print('verified libraries',len(libs),'patients',len(set(mapping.values())),'cells',sum(z['cells'] for z in libs))
print('measured genes per library',pd.DataFrame(rows).groupby('sample').measured.sum().to_dict())
