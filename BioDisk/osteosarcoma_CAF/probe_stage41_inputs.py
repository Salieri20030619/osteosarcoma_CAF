from pathlib import Path
import anndata as ad, pandas as pd, json, importlib.util
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for sample in ['SCPCS000522','SCPCS000246','GSE152048__BC10']:
 f=P/'phase26_client_reanalysis_20260922/inputs'/sample/'full_called_counts.h5ad'
 a=ad.read_h5ad(f,backed='r')
 print('MATRIX',sample,a.shape,'X',type(a.X).__name__,'LAYERS',list(a.layers.keys()),'OBS',list(a.obs.columns),'VAR',list(a.var.columns))
 print('IDS',list(a.obs_names[:3]),list(a.var_names[:3]))
 print('VAR_HEAD',a.var.head(3).to_string())
 a.file.close()
 print('EVIDENCE',pd.read_csv(P/'phase26_client_reanalysis_20260922/tables'/f'{sample}_identity_evidence.tsv.gz',sep='\t',nrows=2).to_string(index=False))
for pkg in ['gseapy','decoupler','scanpy','leidenalg','igraph','statsmodels','scipy','requests']:
 print('PACKAGE',pkg,bool(importlib.util.find_spec(pkg)))
print('MANIFEST',pd.read_csv(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',sep='\t').groupby('dataset').size().to_dict())