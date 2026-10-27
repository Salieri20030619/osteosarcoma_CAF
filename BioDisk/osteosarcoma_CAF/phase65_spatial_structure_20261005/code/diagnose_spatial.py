from pathlib import Path
import pandas as pd,numpy as np,anndata as ad
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
d=pd.read_csv(B/'phase46_spatial_functional_20261001/tables/all_spot_RNA_and_components.tsv.gz',sep='\t',usecols=['sample','library_UMI','spot_id'])
sm='MR4';t=d[d['sample']==sm];a=ad.read_h5ad(B/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad')
lib=np.asarray(a.X.sum(1)).ravel();print('sm',sm,'n',len(t),'zero',sum(lib==0),'same',np.array_equal(lib,t.library_UMI),'obs',a.obs_names.tolist()==t.spot_id.tolist())
print('zeroids',t.loc[lib==0,'spot_id'].tolist());print('diff',np.max(np.abs(lib-t.library_UMI.to_numpy())))
