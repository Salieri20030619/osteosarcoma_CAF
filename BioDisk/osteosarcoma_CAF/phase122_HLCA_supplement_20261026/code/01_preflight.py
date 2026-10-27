from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import pandas as pd,numpy as np,json,hashlib,importlib.metadata as md,zipfile,requests
from scipy import sparse
import anndata as ad
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=R/'phase119_Bone_Lung_functional_subgroups_20261024';O=R/'phase122_HLCA_supplement_20261026'
for sub in ['code','logs','objects','tables','resources','reference','delivery']:(O/sub).mkdir(exist_ok=True,parents=True)
versions={}
for name in ['numpy','pandas','anndata','scanpy','scvi-tools','scarches','torch','cellxgene-census','tiledbsoma','h5py','requests']:
 try:versions[name]=md.version(name)
 except md.PackageNotFoundError:versions[name]=None
print('VERSIONS',json.dumps(versions),flush=True)
(O/'resources/SOFTWARE_PREFLIGHT.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
for url in ['https://zenodo.org/records/7599104/files/HLCA_full_v1.1_emb.h5ad?download=1','https://api.cellxgene.cziscience.com/curation/v1/collections']:
 try:
  r=requests.get(url,stream=True,timeout=(8,15));print('HTTP',url,r.status_code,r.headers.get('Content-Length'),flush=True);r.close()
 except Exception as e:print('HTTP_ERROR',url,str(e),flush=True)
m=pd.read_csv(I/'tables/cohort_registry.tsv.gz',sep='\t');labels=pd.read_csv(I/'tables/hierarchical_cell_labels.tsv.gz',sep='\t');assert m.cell_id.equals(labels.cell_id);m['parent']=labels.parent;m['leaf']=labels.leaf;fit=~m.doublet
genes=pd.read_csv(I/'tables/genes.tsv',sep='\t').symbol.to_numpy();model=pd.read_csv(O/'reference/HLCA_reference_model_gene_order_ids_and_symbols.csv',index_col=0);print('MODEL',model.shape,model.columns.tolist(),model.head(2).to_dict(),flush=True)
can=pd.read_csv(R/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');shared=set(genes);modelgenes=model.gene_symbol.tolist();assert len(modelgenes)==2000
parts=[];order=[];coverage=[]
objects=json.loads((I/'resources/OBJECTS.json').read_text())
for item in objects:
 path=Path(item['path']);src=item['source'];a=ad.read_h5ad(path,backed='r');ix=m.index[fit&m.source.eq(src)&m.cell_id.isin(a.obs_names)].to_numpy()
 if not len(ix):a.file.close();continue
 if src!='NCI':mp={r.symbol:a.var_names.get_indexer([r.feature])[0] for r in can[can.dataset.eq(src)].itertuples()}
 else:
  sy=a.var.symbol.astype(str);ok=~sy.duplicated(keep=False);mp={s:j for j,s in enumerate(sy) if ok.iloc[j]}
 ids=a.obs_names.get_indexer(m.loc[ix,'cell_id']);raw=a.X[ids,:].tocsr();den=np.asarray(raw.sum(1)).ravel();assert np.allclose(den,m.loc[ix,'library_counts']);assert np.max(abs(raw.data-np.round(raw.data)),initial=0)<1e-5
 valid=[g for g in modelgenes if mp.get(g,-1)>=0];missing=sorted(set(modelgenes)-set(valid));assert len(valid)>0
 cols=[mp[g] for g in valid];outcols=[modelgenes.index(g) for g in valid];coo=raw[:,cols].tocoo();out=sparse.csr_matrix((coo.data,(coo.row,np.array(outcols)[coo.col])),shape=(len(ix),2000));parts.append(out);order.extend(ix);coverage.append(dict(source=src,path=str(path),cells=len(ix),measured_model_genes=len(valid),missing_genes=missing));a.file.close();print('COUNTS',src,len(ix),len(valid),flush=True)
assert len(order)==int(fit.sum())==6423 and len(set(order))==6423
allix=m.index[fit].to_numpy();X=sparse.vstack(parts).tocsr()[pd.Index(order).get_indexer(allix)];obs=m.loc[allix].copy().set_index('cell_id',drop=False);obs['dataset']=obs.source.astype(str)+'__'+obs['sample'].astype(str);obs['batch']=obs['dataset'];var=model.copy();var.index=var.index.astype(str);query=ad.AnnData(X=X.astype(np.float32),obs=obs,var=var);query.write_h5ad(O/'objects/query_model_counts.h5ad',compression='gzip')
m.to_csv(O/'tables/frozen_registry.tsv.gz',sep='\t',index=False)
with zipfile.ZipFile(O/'reference/HLCA_reference_model.zip') as z:
 assert z.testzip() is None
 for name in z.namelist():
  if name.startswith('HLCA_reference_model/') and not name.endswith('/') and '.DS_Store' not in name:z.extract(name,O/'reference')
r=dict(status='PASS_RAW_QUERY_PREFLIGHT',main_cells=6423,model_genes=2000,coverage=coverage,missing_padding='Zeros for model architecture only; missing coverage explicitly recorded and not treated as measured zeros in expression inference',input_plan_sha256=hashlib.sha256((O/'PLAN_ZH.txt').read_bytes()).hexdigest());(O/'QUERY_PREFLIGHT.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True)
