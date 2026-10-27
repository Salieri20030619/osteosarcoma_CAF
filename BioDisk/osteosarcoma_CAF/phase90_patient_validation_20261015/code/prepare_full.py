from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']: os.environ[k]='4'
import json,gzip,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite,mmread
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase90_patient_validation_20261015'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def dump(p,v): p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
panels=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text(encoding='utf-8'))
manifest=pd.read_csv(P/'phase87_NCC_preprocessing_20261014/manifest.tsv',sep='\t').fillna('')
manifest=manifest.assign(priority=manifest['sample'].ne('NCC__BS7').astype(int)).sort_values(['priority','sample'])
manifest.to_csv(O/'resources/manifest.tsv',sep='\t',index=False,encoding='utf-8')
checks=[]
for r in manifest.to_dict('records'):
 sm=r['sample'];d=O/'inputs'/sm;d.mkdir(exist_ok=True)
 src=P/'phase87_NCC_preprocessing_20261014/inputs'/sm/'counts.h5ad'
 a=ad.read_h5ad(src);x=a.X.tocsr();den=np.asarray(x.sum(1)).ravel()
 assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all()
 assert a.obs_names.is_unique and a.var_names.is_unique and (den>0).all()
 symbols=a.var.symbol.astype(str);target=sorted(set(sum(panels.values(),[])))
 gm={}
 for j,sy in enumerate(symbols):
  for g in sy.split(';'):
   if g in target:gm.setdefault(g,[]).append(j)
 raw=np.zeros((len(a),len(target)),dtype=np.float32)
 for j,g in enumerate(target):
  if g in gm:raw[:,j]=np.asarray(x[:,gm[g]].sum(1)).ravel()
 log=np.log1p(raw/np.maximum(den,1)[:,None]*1e4)
 scores={k:log[:,[target.index(g) for g in gs if g in gm]].mean(1) for k,gs in panels.items()}
 v=np.column_stack(list(scores.values()));mx=v.max(1)
 label=np.asarray(list(panels),dtype=object)[v.argmax(1)]
 label[(mx<=0)|(np.sum(v==mx[:,None],axis=1)!=1)]='Unresolved'
 cells=pd.DataFrame({'cell_id':sm+'_'+a.obs_names.astype(str),'barcode':a.obs_names.astype(str),'tool_id':[sm+'_full_c'+str(i).zfill(6) for i in range(len(a))], 'sample':sm,'patient_id':r['patient_id'],'n_counts':den,'n_genes':np.asarray((x>0).sum(1)).ravel(),'provisional_lineage':label})
 old=pd.read_csv(P/'phase88_expansion_identity_20261014/inputs'/sm/'cells.tsv',sep='\t').set_index('cell_id')
 cells['phase88_retained']=cells.cell_id.isin(old.index)
 cells['phase88_cnv_role']=cells.cell_id.map(old.cnv_role).fillna('not_in_phase88')
 cells['reference']=cells.phase88_cnv_role.str.startswith('ref_')
 cells['heldout']=cells.phase88_cnv_role.str.startswith('holdout_')
 for k,z in scores.items():cells['phase27_'+k]=z
 cells.to_csv(d/'cells.tsv',sep='\t',index=False,encoding='utf-8')
 pd.DataFrame({'ensembl':a.var_names.astype(str),'symbol':symbols.to_numpy()}).to_csv(d/'genes.tsv',sep='\t',index=False,encoding='utf-8')
 with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
 back=mmread(d/'counts.mtx.gz').T.tocsr();assert back.shape==x.shape and (back!=x).nnz==0
 pd.DataFrame({'cell':cells.tool_id,'group':np.where(cells.reference,'reference','query')}).to_csv(d/'annotations.tsv',sep='\t',index=False,header=False,encoding='utf-8')
 overlap=cells.set_index('cell_id').loc[old.index]
 assert np.array_equal(overlap.n_counts.to_numpy(),old.n_counts.to_numpy())
 entry={'sample':sm,'source':str(src),'source_sha256':sha(src),'full_cells':len(a),'phase88_retained':len(old),'additional_full_cells':len(a)-len(old),'references':int(cells.reference.sum()),'heldout':int(cells.heldout.sum()),'matrix_roundtrip':True,'overlap_count_totals_equal':True}
 checks.append(entry);dump(d/'INPUT_COMPLETE.json',entry);dump(O/'resources/FULL_INPUT_ACCEPTANCE.json',checks)
 print('PREPARED',sm,len(a),flush=True)
 assert (len(a)-len(old))>=0
dump(O/'PREPARATION_COMPLETE.json',{'input_cells':sum(r['full_cells'] for r in checks),'samples':len(checks),'new_formal_patients':0})
