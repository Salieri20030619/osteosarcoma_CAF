from pathlib import Path
import json,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for path in ['phase46_spatial_functional_20261001/code/common.py','phase63_matrix_resistance_20261005/resources/PROGRAMS.json']:
 p=B/path
 if p.suffix=='.py':print(p.name,p.read_text()[:15000])
 else:print('program_sets',list(json.loads(p.read_text())['sets']))
for rel in ['phase63_matrix_resistance_20261005/tables/posthoc_collagen_structure_patient_effects.tsv','phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz','phase46_spatial_functional_20261001/tables/all_spot_RNA_and_components.tsv.gz']:
 d=pd.read_csv(B/rel,sep='\t',nrows=2);print(rel,list(d.columns));print(d.iloc[:1,:12].to_string(index=False))
