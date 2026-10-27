from pathlib import Path
import pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');m=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
for row in m[m.site.isin(['Paraspinal','Calcaneus'])].itertuples(index=False):
 print('SITE_RECORD',row.sample,row.dataset,row.patient_id,row.site,row.clinical_group)
 a=ad.read_h5ad(P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad',backed='r')
 for c in ['tissue_location','tissue_ontology_term_id','diagnosis','subdiagnosis','provider_consensus_celltype_annotation']:
  if c in a.obs:print(c,list(a.obs[c].astype(str).unique())[:6])
 a.file.close()