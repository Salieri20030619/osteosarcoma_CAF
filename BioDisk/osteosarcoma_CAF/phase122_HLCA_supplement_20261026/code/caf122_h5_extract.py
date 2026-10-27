from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import fsspec,h5py,pandas as pd,numpy as np,json,time,anndata as ad
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';I=R/'phase119_Bone_Lung_functional_subgroups_20261024'
url='https://datasets.cellxgene.cziscience.com/5f863718-02a6-46cb-83e7-52be988b915b.h5ad'
programs=json.loads((I/'resources/PROGRAMS.json').read_text());programs={k:programs[k] for k in ['collagen_biosynthesis','collagen_crosslinking','fibrillar_I_III_V','mCAF_matrix','F2C3_collagen_support','myCAF_customer','Osteogenic_customer','mineralization','inflammation']}
programs.update(HLCA_IPF_reference=['CCL2','COL1A1','CTHRC1','MMP19','SERPINE1','HIF1A'],collagen_processing_nonCOL=['SERPINH1','P4HA1','P4HA2','P4HA3','PLOD1','PLOD2','PLOD3','PCOLCE','PCOLCE2','BMP1','LOX','LOXL1','LOXL2'],FNDC1_paper=['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3'])
desired=set(sum(programs.values(),[]))&set(pd.read_csv(I/'tables/genes.tsv',sep='\t').symbol)
def col(g):
 if isinstance(g,h5py.Group):
  cats=g['categories'].asstr()[:];codes=g['codes'][:];a=np.asarray(cats,dtype=object)[np.maximum(codes,0)];a[codes<0]='Unknown';return a
 return g.asstr()[:] if g.dtype.kind in 'OSU' else g[:]
local_reference=O/'reference/HLCA_full_expression.h5ad'
with (local_reference.open('rb') if local_reference.exists() else fsspec.open(url,'rb',block_size=1048576,cache_type='blockcache')) as f:
 with h5py.File(f,'r') as h:
  obs=pd.DataFrame({k:col(h['obs'][k]) for k in ['donor_id','study','dataset','lung_condition','ann_finest_level','total_counts']});obs['row']=np.arange(len(obs))
  types=['Adventitial fibroblasts','Alveolar fibroblasts','Myofibroblasts','Peribronchial fibroblasts','Subpleural fibroblasts','Pericytes','Smooth muscle','Smooth muscle FAM83D+','SM activated stress response']
  obs=obs[obs.ann_finest_level.isin(types)&obs.lung_condition.isin(['Healthy','Healthy (tumor adjacent)','IPF'])]
  obs.to_csv(O/'tables/reference_http_all_stromal_obs.tsv.gz',sep='\t',index=False)
  sel=obs.groupby(['study','donor_id','lung_condition','ann_finest_level'],group_keys=False).apply(lambda g:g.sample(min(100,len(g)),random_state=20261026),include_groups=True).sort_values('row').reset_index(drop=True)
  sel.to_csv(O/'tables/reference_http_selected_obs.tsv.gz',sep='\t',index=False)
  print('SELECTED',len(sel),sel.groupby('lung_condition').size().to_dict(),flush=True)
  raw=h['raw/X'];print('RAW',dict(raw.attrs),list(h['raw/var']),flush=True)
  names=pd.Index(col(h['raw/var/feature_name']) if 'feature_name' in h['raw/var'] else col(h['var/feature_name']));assert len(names)==raw.attrs['shape'][1]
  genes=sorted(desired&set(names[~names.duplicated(keep=False)]));name_map={v:i for i,v in enumerate(names) if v in genes};inds=np.array([name_map[g] for g in genes]);lookup=np.full(len(names),-1,dtype=int);lookup[inds]=np.arange(len(genes));indptr=raw['indptr'][:]
  sums=np.zeros(len(sel));matrix=np.zeros((len(sel),len(genes)),dtype=np.float32);cache=O/'objects/reference_http_partial.npz';done=0
  if cache.exists():
   z=np.load(cache);assert z['matrix'].shape==matrix.shape;matrix=z['matrix'];sums=z['sums'];done=int(z['done'])
  for j in range(done,len(sel)):
   row=int(sel.row.iloc[j]);lo,hi=map(int,indptr[row:row+2]);v=raw['data'][lo:hi];ci=raw['indices'][lo:hi]
   assert np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all(),'Raw matrix not nonnegative integer counts'
   sums[j]=v.sum();mapped=lookup[ci];keep=mapped>=0;np.add.at(matrix[j],mapped[keep],v[keep])
   if (j+1)%100==0:
    np.savez_compressed(cache,matrix=matrix,sums=sums,done=j+1);print('ROWS',j+1,len(sel),flush=True)
  assert (sums>0).all();sel['raw_sum']=sums
  a=ad.AnnData(sparse.csr_matrix(matrix),obs=sel.copy(),var=pd.DataFrame(index=genes));a.write_h5ad(O/'objects/reference_http_program_counts.h5ad',compression='gzip')
  norm=np.log1p(matrix*10000/sums[:,None]);coverage={}
  for name,gs in programs.items():
   got=[g for g in gs if g in genes];coverage[name]=dict(defined=len(gs),measured=len(got),missing=sorted(set(gs)-set(got)));sel[name]=norm[:,pd.Index(genes).get_indexer(got)].mean(1) if len(got)/len(gs)>=.8 else np.nan
  sel.to_csv(O/'tables/reference_http_cell_program_scores.tsv.gz',sep='\t',index=False)
  (O/'resources/SUPPLEMENT_PROGRAMS.json').write_text(json.dumps(programs,indent=2),encoding='utf-8');(O/'resources/REFERENCE_HTTP_COVERAGE.json').write_text(json.dumps(coverage,indent=2),encoding='utf-8')
  (O/'REFERENCE_HTTP_COMPLETE.json').write_text(json.dumps(dict(status='PASS',url=url,cells=len(sel),genes=len(genes),cap_per_donor_type_condition=100,seed=20261026,normalization='log1p counts per10000 using full raw library',raw_integer_check=True),indent=2),encoding='utf-8');print('COMPLETE',flush=True)
