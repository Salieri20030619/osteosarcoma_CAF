from pathlib import Path
import json,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for name in ['phase34_reference_rerun_20260927','phase32_lung_sensitivity_20260926','phase54_raw_review_20261004','phase62_F1F5_internal_20261005','phase49_audit_revision_20261002']:
 p=B/name;print('DIR',name)
 for sub in ['', 'QA','resources']:
  d=p/sub
  if d.exists():print(sub,[x.name for x in d.iterdir() if x.is_file()][:40])
for p in B.glob('phase34*/**/*STATUS*.json'):
 print('STATUS',str(p),p.read_text()[:4500])
d=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False)
for c in ['doublet_assessed_any','doublet_positive_any','CNV_evidence_current','lineage_confidence','decontx_contamination']:
 print('FIELD',c,d[c].value_counts(dropna=False).to_dict() if c!='decontx_contamination' else d[c].describe().to_dict())
print('DOUBLET_BY_SOURCE',d.groupby('dataset')[['doublet_assessed_any','doublet_positive_any']].agg(['count','sum']).to_string())
