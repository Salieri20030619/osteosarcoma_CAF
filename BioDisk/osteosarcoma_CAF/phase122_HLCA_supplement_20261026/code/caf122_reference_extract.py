from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import cellxgene_census,tiledbsoma,pandas as pd,numpy as np,anndata as ad,json,h5py,hashlib
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';I=R/'phase119_Bone_Lung_functional_subgroups_20261024'
ctx=tiledbsoma.SOMATileDBContext(tiledb_config={'vfs.s3.region':'us-west-2','vfs.s3.no_sign_request':'true','sm.compute_concurrency_level':'2','sm.io_concurrency_level':'2','vfs.s3.max_parallel_ops':'1'})
dataset='9f222629-9e39-47d0-b83f-e08d610c7479';types=['CL:4028006','CL:4028004','CL:2000093','CL:0000057','CL:0009089','CL:0000186','CL:0000192','CL:0019019','CL:0000669']
filter_obs=f"dataset_id == '{dataset}' and cell_type_ontology_term_id in {types} and disease in ['normal','pulmonary fibrosis']"
programs0=json.loads((I/'resources/PROGRAMS.json').read_text());names=['collagen_biosynthesis','collagen_crosslinking','fibrillar_I_III_V','mCAF_matrix','F2C3_collagen_support','myCAF_customer','Osteogenic_customer','mineralization','inflammation'];programs={k:programs0[k] for k in names};programs['HLCA_IPF_reference']=['CCL2','COL1A1','CTHRC1','MMP19','SERPINE1','HIF1A'];programs['collagen_processing_nonCOL']=['SERPINH1','P4HA1','P4HA2','P4HA3','PLOD1','PLOD2','PLOD3','PCOLCE','PCOLCE2','BMP1','LOX','LOXL1','LOXL2'];programs['FNDC1_paper']=['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3']
targetgenes=set(pd.read_csv(I/'tables/genes.tsv',sep='\t').symbol);desired=sorted(set(sum(programs.values(),[]))&targetgenes);print('FEATURES',len(desired),flush=True)
with cellxgene_census.open_soma(census_version='2025-11-08',context=ctx) as census:
 exp=census['census_data']['homo_sapiens'];obs=exp.obs.read(value_filter=filter_obs,column_names=['soma_joinid','dataset_id','donor_id','disease','cell_type','cell_type_ontology_term_id','raw_sum','is_primary_data','assay','tissue']).concat().to_pandas();print('OBS',len(obs),obs.groupby(['cell_type','disease']).size().to_dict(),flush=True);obs.to_csv(O/'tables/reference_census_selected_obs.tsv.gz',sep='\t',index=False)
 var=exp.ms['RNA'].var.read().concat().to_pandas();print('VAR_COLUMNS',var.columns.tolist(),flush=True);var=var[var.feature_name.isin(desired)].copy();var=var[~var.feature_name.duplicated(keep=False)];available=set(var.feature_name);coverage={k:dict(defined=len(v),measured=len(set(v)&available),missing=sorted(set(v)-available)) for k,v in programs.items()};(O/'resources/REFERENCE_PROGRAM_COVERAGE.json').write_text(json.dumps(coverage,indent=2),encoding='utf-8')
 # Limit only per-donor very large strata; preserve at least20 cells and all smaller strata.
 sel=obs.groupby(['donor_id','disease','cell_type'],group_keys=False,observed=True).apply(lambda g:g.sample(n=min(500,len(g)),random_state=20261026),include_groups=True).sort_values('soma_joinid').reset_index(drop=True)
 print('SELECTED',len(sel),sel.donor_id.nunique(),flush=True)
 cache=O/'objects/reference_program_counts.h5ad'
 if not cache.exists():
  adata=cellxgene_census.get_anndata(census,organism='Homo sapiens',obs_coords=sel.soma_joinid.to_numpy(),var_coords=var.soma_joinid.to_numpy(),X_name='raw',obs_column_names=['soma_joinid','donor_id','disease','cell_type','raw_sum','dataset_id','is_primary_data'],var_column_names=['soma_joinid','feature_id','feature_name']);adata.write_h5ad(cache,compression='gzip')
 else:adata=ad.read_h5ad(cache)
 assert adata.n_obs==len(sel) and set(adata.obs.soma_joinid)==set(sel.soma_joinid);assert np.isfinite(adata.X.data).all() and (adata.X.data>=0).all() and np.equal(adata.X.data,np.floor(adata.X.data)).all();assert np.all(adata.obs.raw_sum.to_numpy()>0)
 den=adata.obs.raw_sum.to_numpy();norm=adata.X.astype(np.float64).multiply((10000/den)[:,None]).tocsr();norm.data=np.log1p(norm.data);features=pd.Index(adata.var.feature_name);sc=adata.obs.copy()
 for name,gs in programs.items():
  got=[g for g in gs if g in features];sc[name]=np.asarray(norm[:,features.get_indexer(got)].mean(1)).ravel() if len(got)/len(gs)>=.8 else np.nan
 sc.to_csv(O/'tables/reference_cell_program_scores.tsv.gz',sep='\t',index=False)
 # Use original HLCA donor-level condition metadata to separate IPF from other pulmonary fibrosis.
 hp=R/'phase25_complete_nonPDX_assets_20260921/references/HLCA/HLCA_full_v1.1_emb.h5ad'
 with h5py.File(hp) as h:
  def col(k):
   g=h['obs'][k]
   if isinstance(g,h5py.Group):
    a=g['categories'].asstr()[:];c=g['codes'][:];out=np.asarray(a,dtype=object)[np.maximum(c,0)];out[c<0]='Unknown';return out
   return g.asstr()[:]
  metadata=pd.DataFrame({k:col(k) for k in ['donor_id','study','dataset','lung_condition']}).drop_duplicates()
 metadata.to_csv(O/'tables/reference_donor_condition_dictionary.tsv',sep='\t',index=False)
 matches=metadata.groupby('donor_id').agg(studies=('study',lambda v:';'.join(sorted(set(v)))),conditions=('lung_condition',lambda v:';'.join(sorted(set(v)))),datasets=('dataset',lambda v:';'.join(sorted(set(v)))))
 sc=sc.merge(matches,left_on='donor_id',right_index=True,how='left',validate='many_to_one');sc['matched_condition']=sc.conditions.notna();sc.to_csv(O/'tables/reference_cell_program_scores.tsv.gz',sep='\t',index=False)
 (O/'resources/SUPPLEMENT_PROGRAMS.json').write_text(json.dumps(programs,indent=2),encoding='utf-8')
 record=dict(status='PASS_REFERENCE_EXTRACTION',census_version='2025-11-08',dataset_id=dataset,available_stromal_cells=len(obs),extracted_cells=len(sc),selected_donor_records=int(sc.donor_id.nunique()),original_condition_metadata_match_fraction=float(sc.matched_condition.mean()),per_donor_type_cap=500,seed=20261026,primary_data_filter='Not applied; only single integrated HLCA dataset used; no core/extension duplicate concatenation',reference_counts_sha256=hashlib.sha256(cache.read_bytes()).hexdigest());(O/'REFERENCE_EXPRESSION_COMPLETE.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record),flush=True)
