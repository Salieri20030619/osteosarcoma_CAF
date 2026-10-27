from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import pandas as pd,numpy as np,anndata as ad,json,shutil,importlib,datetime
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';I=R/'phase119_Bone_Lung_functional_subgroups_20261024'
m=pd.read_csv(O/'tables/frozen_registry.tsv.gz',sep='\t');model=pd.read_csv(O/'reference/HLCA_reference_model_gene_order_ids_and_symbols.csv',index_col=0);refids=pd.Index(model.index.astype(str));full=ad.read_h5ad(R/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',backed='r');current_symbols=full.var.symbol.reindex(refids);full.file.close();objs=json.loads((I/'resources/OBJECTS.json').read_text());parts=[];order=[];cov=[]
for item in objs:
 a=ad.read_h5ad(item['path'],backed='r');src=item['source'];ix=m.index[~m.doublet&m.source.eq(src)&m.cell_id.isin(a.obs_names)].to_numpy()
 if not len(ix):a.file.close();continue
 mp={v:i for i,v in enumerate(a.var_names)};symbols=a.var.symbol.astype(str);symok=~symbols.duplicated(keep=False);symmap={v:i for i,v in enumerate(symbols) if symok.iloc[i]};cols=[];methods=[]
 for gid,old,new in zip(refids,model.gene_symbol,current_symbols):
  if gid in mp:cols.append(mp[gid]);methods.append('stable_ensembl')
  elif isinstance(new,str) and new in symmap:cols.append(symmap[new]);methods.append('ensembl_to_current_symbol')
  elif old in symmap:cols.append(symmap[old]);methods.append('published_symbol')
  else:cols.append(-1);methods.append('absent_in_available_input')
 cols=np.array(cols);ok=cols>=0;raw=a.X[a.obs_names.get_indexer(m.loc[ix,'cell_id']),:].tocsr();assert np.allclose(np.asarray(raw.sum(1)).ravel(),m.loc[ix,'library_counts']);sub=raw[:,cols[ok]].tocoo();out=sparse.csr_matrix((sub.data,(sub.row,np.flatnonzero(ok)[sub.col])),shape=(len(ix),2000));parts.append(out);order.extend(ix)
 cov.append(dict(source=src,sample=';'.join(m.loc[ix,'sample'].unique()),cells=len(ix),present=int(ok.sum()),absent=int((~ok).sum()),mapping_methods=pd.Series(methods).value_counts().to_dict(),missing_ids=refids[~ok].tolist()))
 pd.DataFrame(dict(gene_id=refids,legacy_symbol=model.gene_symbol,current_symbol=current_symbols,column=cols,method=methods)).to_csv(O/'tables'/('model_feature_map_'+src+'_'+str(len(cov))+'.tsv'),sep='\t',index=False)
 a.file.close()
indices=m.index[~m.doublet];assert len(order)==len(set(order))==6423;X=sparse.vstack(parts).tocsr()[pd.Index(order).get_indexer(indices)];obs=m.loc[indices].copy().set_index('cell_id',drop=False);obs['dataset']='CAF122_'+obs.source.astype(str);obs['scanvi_label']='unlabeled';q=ad.AnnData(X.astype(np.float32),obs=obs,var=model.copy());p=O/'objects/query_model_counts.h5ad';old=O/'resources/query_model_counts_symbol_only.h5ad'
if p.exists() and not old.exists():shutil.copyfile(p,old)
q.write_h5ad(p,compression='gzip');(O/'resources/FEATURE_ID_RECONCILIATION.json').write_text(json.dumps(dict(status='PASS',coverage=cov,reference_gene_order_exact=q.var_names.equals(refids)),indent=2),encoding='utf-8');print('GENE_FIX',[(a['source'],a['cells'],a['present']) for a in cov],flush=True)
import cellxgene_census
print('CENSUS_OPEN',flush=True)
with cellxgene_census.open_soma(census_version='stable') as census:
 ds=census['census_info']['datasets'].read().concat().to_pandas();keep=ds[ds.astype(str).apply(lambda col:col.str.contains('6f6d381a|9f222629|066943a2|Human Lung Cell Atlas',case=False,regex=True)).any(axis=1)];keep.to_csv(O/'tables/census_HLCA_datasets.tsv',sep='\t',index=False);print('CENSUS_DATASETS',keep.to_dict('records'),flush=True);print('CENSUS_SCHEMA',census['census_data']['homo_sapiens'].obs.schema,flush=True)
print('COMPLETE_GENE_AND_CENSUS',flush=True)
