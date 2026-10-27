from pathlib import Path
import gzip
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase8_original_plan_20260919')
p=O/'bulk/GPL10295_family.soft.gz.part'
if not p.exists():p=O/'bulk/GPL10295_family.soft.gz'
with gzip.open(p,'rt') as f:
 n=0
 for line in f:
  if line.startswith('!Platform_data_row_count'):print(line.strip())
  if line.startswith('!platform_table_end'):
   print('PLATFORM_TABLE_END_FOUND')
   for j in range(5):print(f.readline().strip())
   break
  n+=1
 print('LINES_TO_END',n)
