from pathlib import Path
import json, anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for root in [P/'phase15_scpca_download_20260920',P/'phase15_scpca_download_20260920/additional23']:
    inv=json.loads((root/'verified_files.json').read_text())
    for r in inv:
        if '_filtered_rna.h5ad' not in r['path']: continue
        a=ad.read_h5ad(r['path'],backed='r')
        print(json.dumps({'path':r['path'],'shape':a.shape,'var':a.var.head(2).to_dict(),'obs_columns':list(a.obs.columns),'layers':list(a.layers),'uns_keys':list(a.uns)},default=str),flush=True)
        a.file.close()
