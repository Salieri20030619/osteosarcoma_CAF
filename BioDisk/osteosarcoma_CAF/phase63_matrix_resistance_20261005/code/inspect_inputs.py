from pathlib import Path
import json,re,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005'
libs=json.loads((P/'phase58_F1F5_round1_20261004/resources/libraries_canonical.json').read_text())
for lib in ['reactome','go_bp','matrisome']:
 for term,genes in libs[lib].items():
  if re.search('COLLAGEN|CROSS.?LINK|EXTRACELLULAR_MATRIX_ORGANIZATION|EXTRACELLULAR_MATRIX_ASSEMBLY|DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX|ELASTIC_FIBRE|ELASTIC_FIBER|LYSYL_OXID|PEPTIDYL_LYSINE_OXID',term):print('SET',lib,term,len(genes))
paths=[P/'phase33_client_report_20260926/tables/final_clinical_sample_manifest.tsv',P/'phase36_external_validation_20260927/resources/sample_patient_manifest_corrected.tsv',P/'phase48_clinical_programs_delivery_20261001/delivery/client/CONDITIONAL_ENDPOINTS.tsv',P/'phase62_F1F5_internal_20261005/tables/F1F5_internal_cell_labels.tsv.gz']
for path in paths:
 f=pd.read_csv(path,sep='\t',low_memory=False);print('TABLE',str(path),'rows',len(f));print('COLUMNS',list(f.columns))
 for col in f.columns:
  if re.search('huvos|necros|response|resistan|treatment|clinical_group|survival',col,re.I):print('VALUES',col,f[col].value_counts(dropna=False).head(12).to_dict())
 if 'CONDITIONAL' in path.name:print(f.to_string(index=False))
for path in [P/'phase46_spatial_functional_20261001/tables/spatial_patient_statistics.tsv',P/'phase50_F3_spatial_decomposition_20261002/tables/statistics.tsv']:
 f=pd.read_csv(path,sep='\t');print('PREVIOUS_SPATIAL',path.name)
 if 'program' in f:f=f[f.program.isin(['NABA_CORE_MATRISOME','F3_SELECTED_REMODELING'])]
 else:f=f[(f.model=='full_adjusted')&f.signature.isin(['LOX','MMP14','full_reference'])]
 print(f.to_string(index=False))
