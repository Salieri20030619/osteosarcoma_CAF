from pathlib import Path
import pandas as pd,numpy as np,json
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase52_blind_lineage_validation_20261002')
panels=json.loads((r/'resources/LINEAGE_PANELS_FROZEN.json').read_text());names=list(panels)
c=pd.read_csv(r/'tables/cell_lineage_scores_and_predictions.tsv.gz',sep='\t')
co=c[c.control_annotation.isin(names)];rows=[]
for scheme,col in [('full_panel','full_label'),('split_panel_consensus','consensus_label')]:
 for (patient,site),t in co.groupby(['patient_id','site']):
  v=t[t[col].eq('Fibroblast')]
  rows.append(dict(scheme=scheme,patient_id=patient,site=site,predicted_Fibroblast_controls=len(v),annotation_Fibroblast_fraction=float(v.control_annotation.eq('Fibroblast').mean()) if len(v) else np.nan))
p=pd.DataFrame(rows);p.to_csv(r/'tables/control_Fibroblast_prediction_site_composition.tsv',sep='\t',index=False)
pp=p.groupby(['scheme','patient_id']).annotation_Fibroblast_fraction.mean().reset_index()
pp.to_csv(r/'tables/control_Fibroblast_prediction_patient_composition.tsv',sep='\t',index=False)
