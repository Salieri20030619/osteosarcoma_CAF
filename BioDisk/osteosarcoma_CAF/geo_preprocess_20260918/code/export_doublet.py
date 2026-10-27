from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import pandas as pd,anndata as ad,gzip,json
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();out=p/'doublet_review';out.mkdir(exist_ok=True)
s=pd.read_csv(p/'review/geo_processing_summary.tsv',sep='\t');s=s[s.doublet_status!='automatic_threshold_pass'];assert len(s)==8
s[['sample']].to_csv(out/'manifest.tsv',sep='\t',index=False)
for sm in s['sample']:
 d=out/sm;d.mkdir(exist_ok=True);a=ad.read_h5ad(p/'results'/sm/'annotated.h5ad')
 with gzip.open(d/'counts.mtx.gz','wb') as f:mmwrite(f,a.layers['counts'].T,field='integer')
 (d/'genes.tsv').write_text('\n'.join(a.var_names)+'\n');(d/'cells.tsv').write_text('\n'.join(a.obs_names)+'\n')
 print('DOUBLETS_INPUT_READY',sm,a.n_obs,flush=True)
(out/'protocol.json').write_text(json.dumps({'selection':'eight samples whose Scrublet auto threshold failed existing rate gate','method':'scDblFinder','mode':'random artificial doublets; cluster-free','dbr':0.06,'dbr_sd':1,'prior_interpretation':'wide uncertainty removes expected-number penalty; threshold driven by classification error, not forcing 6 percent removal','seed':270231,'scope':'per-sample diagnostic and alternate exclusion; no cross-sample pooling','source':'https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html'},indent=2))
