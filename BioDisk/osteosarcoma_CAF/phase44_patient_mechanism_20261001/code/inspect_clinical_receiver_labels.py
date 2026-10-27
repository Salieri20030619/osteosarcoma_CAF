from pathlib import Path
import pandas as pd,anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');m=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
print('CLINICAL',m.groupby(['dataset','clinical_group']).size().to_string());print('SITES',m.groupby(['dataset','site']).size().to_string())
for sm in ['SCPCS000246','SCPCS000522','GSE152048__BC10']:
 a=ad.read_h5ad(P/'phase26_client_reanalysis_20260922/inputs'/sm/'full_called_counts.h5ad',backed='r')
 for col in ['provider_consensus_celltype_annotation','historical_broad_lineage']:
  if col in a.obs:print('LINEAGE_ANNOTATIONS',sm,col,a.obs[col].astype(str).value_counts().to_dict())
 a.file.close()