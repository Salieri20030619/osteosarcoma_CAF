from pathlib import Path
import pandas as pd,json
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for rel in ['phase54_raw_review_20261004/tables/doublet_provenance.tsv','phase54_raw_review_20261004/code/doublet_metadata.py']:
 p=B/rel;print('FILE',rel)
 print(p.read_text()[:12000])
p=B/'phase63_matrix_resistance_20261005/code/mechanism.py'
lines=p.read_text().splitlines()
for i,l in enumerate(lines):
 if any(x in l for x in ['corrected','decont','phase31','read_h5','raw_path']):print('MECH',i+1,'\n'.join(lines[max(0,i-1):i+2]))
