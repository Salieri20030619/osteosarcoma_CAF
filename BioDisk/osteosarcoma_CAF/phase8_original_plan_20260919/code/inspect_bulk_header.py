from pathlib import Path
import gzip,csv,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase8_original_plan_20260919')
p=O/'bulk/GSE21257_series_matrix.txt.gz'
if not p.exists():p=p.with_suffix(p.suffix+'.part')
with gzip.open(p,'rt') as f:
 for line in f:
  if line.startswith('!series_matrix_table_begin'):break
  if line.startswith(('!Sample_characteristics','!Sample_data_processing','!Sample_title','!Sample_platform','!Series_summary','!Series_overall_design')):
   a=next(csv.reader([line],delimiter='\t'))
   print(a[0],json.dumps(list(dict.fromkeys(a[1:])),ensure_ascii=False))
print('HEADER_READ_ONLY_COMPLETE')
