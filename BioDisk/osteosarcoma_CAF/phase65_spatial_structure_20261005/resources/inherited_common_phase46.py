from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,json,hashlib,datetime
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase46_spatial_functional_20261001'
def dump(path,x):path.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
h=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='\t',keep_default_na=False);h=h[h.status.eq('Approved')];approved=set(h.symbol);aliases={};eid={}
for r in h.itertuples(index=False):
 if r.ensembl_gene_id:eid[r.ensembl_gene_id]=r.symbol
 for z in str(r.alias_symbol+'|'+r.prev_symbol).split('|'):
  if z:aliases.setdefault(z,set()).add(r.symbol)
def canon(z):
 z=str(z)
 return z if z in approved else next(iter(aliases[z])) if z in aliases and len(aliases[z])==1 else None
def feature_map(a):
 col='gene_symbol' if 'gene_symbol' in a.var else 'symbol';rows=[]
 for j,(f,s) in enumerate(zip(a.var_names,a.var[col].astype(str))):
  f0=str(f).split('.')[0];g=eid.get(f0) or canon(s)
  if g:rows.append((g,0 if f0 in eid else 1,str(f),j))
 return {g:j for g,pr,f,j in sorted(rows,reverse=True)}
def bh(v):
 v=np.asarray(v,float);ok=np.isfinite(v);x=np.where(ok,v,1);order=np.argsort(x);out=np.empty(len(v));out[order]=np.minimum(1,np.minimum.accumulate((x[order]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]);out[~ok]=np.nan;return out
def test_table(pat,keys,value,family):
 rows=[];vectors=[]
 for key,t in pat.groupby(keys,dropna=False,sort=True):
  key=key if isinstance(key,tuple) else (key,);d=t[value].to_numpy(float);d=d[np.isfinite(d)];n=len(d)
  row=dict(zip(keys,key));row.update(N=n,mean_effect=np.mean(d) if n else np.nan,median_effect=np.median(d) if n else np.nan,positive_patients=int((d>0).sum()),negative_patients=int((d<0).sum()),P=np.nan,test='unassessable' if n<2 or np.all(d==0) else 'exact_patient_signflip')
  rows.append(row);vectors.append(d)
 for n in sorted({len(d) for d in vectors if 2<=len(d)<=20 and not np.all(d==0)}):
  ids=[i for i,d in enumerate(vectors) if len(d)==n and not np.all(d==0)];v=np.column_stack([vectors[i] for i in ids]);obs=np.abs(v.mean(0));hits=np.zeros(len(ids),np.int64)
  for start in range(0,2**n,16384):
   x=np.arange(start,min(start+16384,2**n),dtype=np.uint64);signs=(((x[:,None]>>np.arange(n,dtype=np.uint64))&1)*2).astype(float)-1
   hits+=(np.abs(signs@v/n)>=obs-1e-12).sum(0)
  for i,p in zip(ids,hits/2**n):rows[i]['P']=p
 for i,d in enumerate(vectors):
  if len(d)>20:
   rng=np.random.default_rng(20261001);s=rng.choice([-1.,1.],(65536,len(d)));rows[i]['P']=(1+(np.abs(s@d/len(d))>=abs(d.mean())-1e-12).sum())/65537;rows[i]['test']='MC_patient_signflip_65536'
 out=pd.DataFrame(rows);out['family']=out.apply(family,axis=1);out['FDR']=out.groupby('family').P.transform(bh);out['nominal']=out.P.lt(.05);return out

def corr(x,y,c=None):
 x=np.asarray(x,float);y=np.asarray(y,float);c=np.empty((len(x),0)) if c is None else np.asarray(c,float);ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(c).all(1)
 if ok.sum()<3:return np.nan
 xr=rankdata(x[ok]);yr=rankdata(y[ok])
 if np.ptp(xr)==0 or np.ptp(yr)==0:return np.nan
 if c.shape[1]:
  m=np.column_stack([np.ones(ok.sum())]+[rankdata(z) for z in c[ok].T]);xr-=m@np.linalg.lstsq(m,xr,rcond=None)[0];yr-=m@np.linalg.lstsq(m,yr,rcond=None)[0]
  if ok.sum()<=np.linalg.matrix_rank(m)+1 or np.linalg.norm(xr)<1e-9 or np.linalg.norm(yr)<1e-9:return np.nan
 return float(np.corrcoef(xr,yr)[0,1])
