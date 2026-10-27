from pathlib import Path
import gzip,json
import pandas as pd,numpy as np,scanpy as sc
from scipy.io import mmread
from scipy import sparse
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();root=p/'remediation_20260917';out=root/'decontx_comparison';out.mkdir(exist_ok=True)
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');lab=a.obs['harmony_theta2_r0.5'].astype(str);panels=json.loads((p/'analysis/gene_panels.json').read_text());e=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t',usecols=['cell_id','lineage_reviewed']).set_index('cell_id')
rows=[];strict=[];summary=[]
for src in sorted((root/'decontx_inputs').iterdir()):
 sample=src.name;d=root/'decontx'/sample
 if not (d/'COMPLETE').exists():summary.append({'sample':sample,'status':'failed_or_incomplete'});continue
 cells=pd.read_csv(src/'cells.tsv',sep='\t');genes=pd.read_csv(src/'genes.tsv',sep='\t');ct=pd.read_csv(d/'cell_contamination.tsv',sep='\t');ids=pd.Index(cells.cell_id);assert ids.equals(pd.Index(ct.cell_id))
 with gzip.open(src/'counts.mtx.gz','rb') as f:raw=sparse.csr_matrix(mmread(f).T)
 with gzip.open(d/'corrected_counts.mtx.gz','rb') as f:corrected=sparse.csr_matrix(mmread(f).T)
 assert raw.shape==corrected.shape and raw.shape==(len(cells),len(genes));assert np.isfinite(corrected.data).all() and (corrected.data>=0).all()
 summary.append({'sample':sample,'status':'complete','cells':len(cells),'median_contamination':ct.contamination.median(),'q90_contamination':ct.contamination.quantile(.9)})
 lineage=e.reindex(ids).lineage_reviewed.to_numpy();isstrict=ids.isin(a.obs_names);states=np.array([str(lab.get(x,'not_strict')) for x in ids])
 for version,mat in [('original',raw),('corrected',corrected)]:
  totals=np.asarray(mat.sum(1)).ravel();log=mat.multiply((1e4/np.maximum(totals,1))[:,None]).tocsr();log.data=np.log1p(log.data)
  for lin in sorted(set(lineage)):
   ix=lineage==lin
   for gene in ['COL1A1','SPP1','IBSP','DCN','LUM','CD74','HLA-DRA','CIITA','CD80','CD86']:
    gj=np.where(genes.symbol==gene)[0]
    if not len(gj):continue
    v=np.asarray(mat[ix][:,gj].sum(1)).ravel();rows.append({'sample':sample,'lineage':lin,'version':version,'gene':gene,'cells':int(ix.sum()),'detected_gt0_pct':float((v>0).mean()*100),'count_ge1_pct':float((v>=1).mean()*100),'mean_count':float(v.mean())})
  for st in sorted(set(states[isstrict])):
   ix=states==st
   for name,gs in panels.items():
    gj=np.where(genes.symbol.isin(gs))[0]
    if not len(gj):continue
    v=np.asarray(log[:,gj].mean(1)).ravel();strict.append({'sample':sample,'state':st,'version':version,'program':name,'cells':int(ix.sum()),'covered_genes':len(gj),'mean_log_program':float(v[ix].mean()),'median_log_program':float(np.median(v[ix])),'median_contamination':ct.loc[ix,'contamination'].median()})
 if isstrict.any():
  order=a.obs_names[a.obs_names.isin(ids)];idx=ids.get_indexer(order);geneidx=pd.Index(genes.gene_id).get_indexer(a.var_names);assert (geneidx>=0).all()
  b=sc.AnnData(corrected[idx][:,geneidx],obs=a.obs.loc[order].copy(),var=a.var.copy());b.write_h5ad(out/(sample+'_strict_corrected.h5ad'),compression='gzip')
pd.DataFrame(summary).to_csv(out/'sample_contamination_summary.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(out/'lineage_gene_before_after.tsv',sep='\t',index=False);pd.DataFrame(strict).to_csv(out/'strict_state_programs_before_after.tsv',sep='\t',index=False)
(out/'COMPLETE').write_text('DECONTX_COMPARISON_COMPLETE\n');print(pd.DataFrame(summary).to_string(index=False))
