from pathlib import Path
import scanpy as sc,pandas as pd,numpy as np,gzip
from scipy.io import mmwrite
p=Path.cwd();out=p/'remediation_20260917/decontx_inputs';out.mkdir(exist_ok=True)
e=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t',usecols=['cell_id','sample','lineage_reviewed']).set_index('cell_id')
source=p.parent/'phase5_nonpdx_qc_20260916/results'
for sample in e['sample'].drop_duplicates():
 d=out/sample;d.mkdir(exist_ok=True)
 if (d/'COMPLETE').exists():continue
 a=sc.read_h5ad(source/sample/'annotated.h5ad');z=e.reindex(a.obs_names);assert z.lineage_reviewed.notna().all()
 with gzip.open(d/'counts.mtx.gz','wb') as f:mmwrite(f,a.layers['counts'].T)
 pd.DataFrame({'cell_id':a.obs_names,'cluster':z.lineage_reviewed.to_numpy()}).to_csv(d/'cells.tsv',sep='\t',index=False)
 pd.DataFrame({'gene_id':a.var_names,'symbol':a.var.symbol.to_numpy()}).to_csv(d/'genes.tsv',sep='\t',index=False)
 (d/'COMPLETE').write_text(str(a.shape));print('EXPORT',sample,a.shape,flush=True)
