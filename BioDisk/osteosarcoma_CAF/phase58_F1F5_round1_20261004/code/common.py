from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,itertools,warnings
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';R=P/'phase49_audit_revision_20261002';O=P/'phase58_F1F5_round1_20261004'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def mean_valid(x,axis=0):
 finite=np.isfinite(x);n=finite.sum(axis);out=np.full(np.shape(n),np.nan,dtype=float);np.divide(np.where(finite,x,0).sum(axis),n,out=out,where=n>0);return out
cache={}
def tests(d):
 p=np.full(d.shape[1],np.nan);N=np.isfinite(d).sum(0);method=np.full(d.shape[1],'not_tested_N_lt2',dtype=object)
 if not len(d):return p,N,method
 patterns,inv=np.unique(np.isfinite(d).T,axis=0,return_inverse=True)
 for k,pat in enumerate(patterns):
  n=int(pat.sum());ix=np.flatnonzero(inv==k)
  if n<2:continue
  z=d[np.ix_(pat,ix)];nz=np.any(z!=0,axis=0);method[ix[~nz]]='not_tested_all_zero';ix=ix[nz];z=z[:,nz]
  if not len(ix):continue
  if n not in cache:cache[n]=np.array(list(itertools.product([-1.,1.],repeat=n))) if n<=15 else np.random.default_rng(20261001).choice([-1.,1.],size=(19999,n))
  signs=cache[n];obs=np.abs(z.mean(0))
  for start in range(0,len(ix),128):
   null=np.abs(signs@z[:,start:start+128]/n);hits=(null>=obs[start:start+128][None,:]-1e-12).sum(0);p[ix[start:start+128]]=hits/len(signs) if n<=15 else (hits+1)/(len(signs)+1)
  method[ix]='exact_patient_signflip' if n<=15 else 'patient_signflip_MC19999_plus_observed'
 return p,N,method
def fdr(df):
 df=df.copy();df['FDR']=np.nan
 for fam,g in df.groupby('family',sort=False):
  q=multipletests(g.P.fillna(1).to_numpy(),method='fdr_bh')[1];df.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);return df