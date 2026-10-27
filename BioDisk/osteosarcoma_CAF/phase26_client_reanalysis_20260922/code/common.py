from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='4'
import json,datetime,hashlib,gzip,io,tarfile,subprocess,shutil,gc
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.io import mmread,mmwrite
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase26_client_reanalysis_20260922'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
R='/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def write_json(p,x):
 p=Path(p);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));tmp.replace(p)
def now():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def truth(x):return pd.Series(x).fillna(False).astype(str).str.lower().isin(['true','1','1.0']).to_numpy()
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def clean_frame(df):
 for c in df:
  if df[c].dtype==object or isinstance(df[c].dtype,pd.StringDtype):df[c]=df[c].fillna('not_reported').astype(str).astype(object)
 df.index=pd.Index(df.index.astype(str).to_numpy(dtype=object));return df

def valid_tar(p):
 try:
  with tarfile.open(p,'r:gz') as t:
   for m in t.getmembers():
    if m.isfile():
     with t.extractfile(m) as f:
      for b in iter(lambda:f.read(4*1024*1024),b''):pass
  return True
 except Exception:return False

def obtain_tar(row):
 name=row['source_url'].rsplit('/',1)[1]
 for p in [P/'phase23_loose_first_20260921/raw'/name,P/'phase25_complete_nonPDX_assets_20260921/recovered_GEO'/name,O/'raw'/name]:
  if p.exists() and valid_tar(p):return p
 dest=O/'raw'/name;part=dest.with_name(dest.name+'.part')
 # Truncated historical archives remain untouched. New files get a separate versioned destination.
 base=['curl','--fail','--location','--silent','--show-error','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--retry','3','--retry-delay','10','--continue-at','-','--output',str(part)]
 with (O/'logs'/('source_download_'+row['sample']+'.log')).open('a') as log:
  for proxy in ['http://127.0.0.1:7897','']:
   cmd=base+(['--proxy',proxy] if proxy else ['--proxy','','--noproxy','*'])+[row['source_url']]
   rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
   if rc==0 and valid_tar(part):part.replace(dest);return dest
 raise RuntimeError('Original GEO archive not yet complete; no substitution by historical QC subset')

def load_tar(p):
 with tarfile.open(p) as t:
  members=t.getmembers()
  def extract(s):
   matches=[m for m in members if m.isfile() and m.name.endswith(s)];assert len(matches)==1,(s,len(matches))
   return gzip.decompress(t.extractfile(matches[0]).read())
  features=[s.split('\t') for s in extract('features.tsv.gz').decode().splitlines()];barcodes=extract('barcodes.tsv.gz').decode().splitlines();x=mmread(io.BytesIO(extract('matrix.mtx.gz'))).tocsr()
 keep=np.array([len(f)<3 or f[2]=='Gene Expression' for f in features]);features=[f for f,k in zip(features,keep) if k]
 return ad.AnnData(x[keep].T.tocsr(),obs=pd.DataFrame(index=barcodes),var=pd.DataFrame({'symbol':[f[1] for f in features]},index=[f[0] for f in features]))

def prepare(row):
 sm=row['sample'];d=O/'inputs'/sm;d.mkdir(parents=True,exist_ok=True)
 if (d/'COMPLETE').exists():return json.loads((d/'input_QA.json').read_text())
 if row['kind']=='ScPCA':
  a=ad.read_h5ad(row['counts_path']);a.obs_names=[sm+'_'+str(x) for x in a.obs_names];a.var['symbol']=a.var.gene_symbol.astype(str)
  assert set(a.obs.is_xenograft.astype(str).str.lower())=={'false'} and set(a.obs.is_cell_line.astype(str).str.lower())=={'false'}
  assert set(a.obs.sample_id.astype(str))=={sm}
  proc=ad.read_h5ad(row['processed_path'],backed='r');po=proc.obs.copy();po.index=[sm+'_'+str(x) for x in po.index]
  for c in ['consensus_celltype_annotation','singler_celltype_annotation','cellassign_celltype_annotation','scimilarity_celltype_annotation','is_infercnv_reference','infercnv_total_cnv']:
   if c in po:a.obs['provider_'+c]=po[c].reindex(a.obs_names)
  ref=truth(a.obs.get('provider_is_infercnv_reference',pd.Series(False,index=a.obs_names)))
  infstatus=str(proc.uns.get('infercnv_status','missing'));t=proc.uns.get('infercnv_table')
  if isinstance(t,pd.DataFrame):t.to_csv(d/'provider_infercnv_table.tsv.gz',sep='\t',index=False)
  write_json(d/'provider_evidence.json',{'infercnv_status':infstatus,'provider_reference_cells':int(ref.sum()),'source_object':row['processed_path'],'annotation_not_CAF_identity_truth':True})
  proc.file.close();reference_source='provider_is_infercnv_reference';source=row['counts_path']
 else:
  if row['kind']=='cached_full_counts':a=ad.read_h5ad(row['counts_path']);source=row['counts_path']
  else:
   p=obtain_tar(row);a=load_tar(p);a.obs_names=[sm+'_'+str(x) for x in a.obs_names];source=str(p)
   old=read(P/'phase22_no_UMI_sensitivity_20260921/legacy_UMI_audit.tsv').set_index('sample');q=read(old.loc[sm,'QC_table']).set_index('cell_id');assert set(a.obs_names)==set(q.index);a.obs=q.loc[a.obs_names].copy()
  hist=read(P/'phase17_all_cohort_expanded_20260921/tables/all_cell_layer_evidence.tsv.gz').set_index('cell_id');lineage=hist.original_lineage.reindex(a.obs_names).fillna('unannotated')
  a.obs['historical_broad_lineage']=lineage
  # Biological normal lineages only; no marker-count, UMI, reference-count or heldout gate.
  ref=lineage.isin(['T_NK','B_Plasma','Myeloid','DC','Endo','Mast']).to_numpy();reference_source='historical_broad_lineage_annotation_normal_immune_or_endothelial; not independent validation';infstatus='new_native_HMM_required'
 x=a.X.tocsr().astype(np.float32);x.eliminate_zeros();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all();assert a.obs_names.is_unique
 counts=np.asarray(x.sum(1)).ravel()
 for c in ['n_counts','sum']:
  if c in a.obs:assert np.allclose(counts,pd.to_numeric(a.obs[c]).to_numpy()),(sm,'count mismatch',c)
 if 'ensembl_id' in a.var:ids=a.var.ensembl_id.astype(str).to_numpy()
 else:ids=np.array([str(z).split('.')[0] for z in a.var_names])
 if len(set(ids))!=len(ids):
  codes,unique=pd.factorize(ids,sort=False);mapping=sparse.csr_matrix((np.ones(len(ids)),(np.arange(len(ids)),codes)),shape=(len(ids),len(unique)));x=(x@mapping).tocsr();sy=pd.Series(a.var.symbol.astype(str).to_numpy()).groupby(codes).agg(lambda z:';'.join(dict.fromkeys(z))).reindex(range(len(unique))).to_numpy();ids=np.asarray(unique)
 else:sy=a.var.symbol.astype(str).to_numpy()
 assert np.allclose(np.asarray(x.sum(1)).ravel(),counts)
 obs=a.obs.copy();obs['sample']=sm;obs['dataset']=row['dataset'];obs['n_counts']=counts;obs['n_genes']=x.getnnz(1);obs['reference']=ref;obs['tool_id']=[sm+'_'+str(i).zfill(7) for i in range(len(obs))];obs['reference_source']=reference_source
 b=ad.AnnData(x,obs=clean_frame(obs),var=clean_frame(pd.DataFrame({'symbol':sy},index=ids)));b.write_h5ad(d/'full_called_counts.h5ad',compression='gzip')
 cells=pd.DataFrame({'cell_id':b.obs_names,'tool_id':b.obs.tool_id.to_numpy(),'n_counts':counts,'reference':ref});cells.to_csv(d/'cells.tsv',sep='\t',index=False)
 pd.DataFrame({'ensembl':ids,'symbol':sy}).to_csv(d/'genes.tsv',sep='\t',index=False)
 with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
 pd.DataFrame({'id':cells.tool_id,'group':np.where(ref,'reference','query')}).to_csv(d/'annotations.tsv',sep='\t',index=False,header=False)
 qa={'sample':sm,'dataset':row['dataset'],'cells':len(b),'genes':b.n_vars,'reference_cells':int(ref.sum()),'reference_source':reference_source,'provider_infercnv_status':infstatus,'source':source,'source_sha256':digest(source),'counts_sha256':digest(d/'counts.mtx.gz'),'cell_filter_applied_by_project':False,'all_provider_called_cells_retained':True,'software_native_filters_apply_in_R_only':True,'nonPDX_evidence':row['nonPDX_evidence'],'created_at':now()};write_json(d/'input_QA.json',qa);(d/'COMPLETE').write_text('FULL_INPUT_COUNTS_VERIFIED\n');return qa
