from pathlib import Path
import gzip,csv,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase8_original_plan_20260919')
p=O/'bulk/GPL10295_family.soft.gz'
if not p.exists():p=p.with_suffix(p.suffix+'.part')
with gzip.open(p,'rt') as f:
 hit=False
 for line in f:
  if line.startswith('!platform_table_begin'):hit=True;continue
  if hit:
   print(line[:6000].rstrip())
   if not line.startswith('ID'):break
print('PLATFORM_HEADER_READ')
