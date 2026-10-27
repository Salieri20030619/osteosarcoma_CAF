from pathlib import Path
import os,json
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import pandas as pd,numpy as np,anndata as ad
from scipy.io import mmread
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918');Q=P.parent/'geo_preprocess_20260918'
s=pd.read_csv(P/'review/input_summary.tsv',sep='\t');assert len(s)==16
panels=json.loads((P/'sources/author_code_marker_panels.json').read_text());rows=[];checks=[]
for r in s.to_dict('records'):
 sm=r['sample'];a=ad.read_h5ad(Q/'accepted'/sm/'analysis_ready.h5ad');x=a.layers['counts'].tocsr()
 if r['cnv_eligible']:
  d=P/'inputs'/sm;o=pd.read_csv(d/'cells.tsv',sep='\t');g=pd.read_csv(d/'genes.tsv',sep='\t');m=mmread(str(d/'counts.mtx.gz')).T.tocsr()
  assert m.shape==x.shape and (m!=x).nnz==0 and o.cell_id.tolist()==a.obs_names.tolist() and g.source_gene_id.tolist()==a.var_names.tolist()
  assert o.tool_id.is_unique and g.ensembl_id.is_unique and np.array_equal(np.asarray(m.sum(1)).ravel(),o.n_counts.to_numpy())
  refs=o.cnv_role.str.startswith('ref_');holds=o.cnv_role.str.startswith('holdout_');assert refs.sum()>=30 and holds.sum()>=20
  assert (o.loc[refs|holds,'pctUMI_competing_mes']<=.2).all()
  checks.append(dict(sample=sm,export_exactly_equals_accepted_counts=True,reference_cells=int(refs.sum()),holdout_cells=int(holds.sum()),target_cells=int((o.cnv_role=='target').sum())))
 if r['dataset']=='GSE318840':
  x=x.tocsc();syms=a.var.symbol.astype(str).to_numpy();obs=a.obs
  for panel,genes in panels.items():
   det=sum((np.asarray(x[:,np.flatnonzero(syms==gene)].sum(1)).ravel()>0).astype(int) for gene in genes)
   for cl in obs.cluster.astype(str).unique():
    mask=(obs.cluster.astype(str)==cl).to_numpy();rows.append(dict(sample=sm,cluster=cl,cells=int(mask.sum()),current_label=str(obs.loc[mask,'lineage_reviewed_candidate'].iloc[0]),author_code_panel=panel,pct_two_detected=float((det[mask]>=2).mean()*100)))
 print('EXPORT_VALIDATED',sm,flush=True)
pd.DataFrame(rows).to_csv(P/'review/FLEX_author_code_marker_support.tsv',sep='\t',index=False)
(P/'review/export_integrity.json').write_text(json.dumps(checks,indent=2));(P/'EXPORT_VALIDATED').write_text('GEO_CNV_EXPORT_EXACT_COUNTS_VERIFIED\n')
print('CNV_EXPORT_VALIDATION_COMPLETE',len(checks),flush=True)
