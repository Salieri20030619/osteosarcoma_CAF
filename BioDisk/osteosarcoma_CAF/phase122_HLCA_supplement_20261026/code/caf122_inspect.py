from pathlib import Path
import json,anndata as ad,pandas as pd,h5py,requests
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';I=R/'phase119_Bone_Lung_functional_subgroups_20261024'
print((O/'resources/SOFTWARE_PREFLIGHT.json').read_text())
print('HTTP_LOG', '\n'.join(x for x in (O/'logs/01_preflight.log').read_text().splitlines() if x.startswith(('HTTP','VERSIONS','MODEL'))))
print('COVERAGE',[(a['source'],a['cells'],a['measured_model_genes']) for a in json.loads((O/'QUERY_PREFLIGHT.json').read_text())['coverage']])
objs=json.loads((I/'resources/OBJECTS.json').read_text())
for it in [objs[0],next(x for x in objs if x['source']=='NCI')]:
 a=ad.read_h5ad(it['path'],backed='r');print('RAW',it['source'],a.shape,a.var.head(3).to_dict(),a.var_names[:5].tolist());a.file.close()
can=pd.read_csv(R/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');print('CAN',can.head(3).to_dict())
with h5py.File(R/'phase25_complete_nonPDX_assets_20260921/references/HLCA/HLCA_full_v1.1_emb.h5ad') as h:
 print('HLCA',h['X'].shape,list(h['obs'].keys()))
 for k in ['ann_finest_level','lung_condition','disease','dataset']:
  if k in h['obs']:
   g=h['obs'][k];print('CATEGORIES',k,g['categories'].asstr()[:].tolist() if isinstance(g,h5py.Group) else 'stringdataset')
