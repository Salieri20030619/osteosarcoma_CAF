from pathlib import Path
import json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for sub in ['phase66_doublet_repair_DSC1_20261005','phase63_matrix_resistance_20261005']:
 d=P/sub
 print('OBJECTS',sub,[(str(f.relative_to(d)),f.stat().st_size) for f in (d/'objects').glob('*')])
 for fn in ['resources/PROGRAMS.json','code/mechanism.py','code/current_metadata.py']:
  f=d/fn
  if f.exists():
   if f.suffix=='.json':
    v=json.loads(f.read_text());print('PROGRAM_KEYS',sub,list(v));print('SETS',json.dumps({k:len(s) for k,s in v.get('sets',v).items() if isinstance(s,list)}))
   else:
    lines=f.read_text().splitlines();matches=[i for i,s in enumerate(lines) if any(t in s for t in ['np.load','np.save','corrected','PROGRAMS','def load_current_cells'])]
    print('CODE',sub,fn,'\n'.join(f'{i+1}:{lines[i]}' for i in matches[:35]))
print('META_LOADER',(P/'phase66_doublet_repair_DSC1_20261005/code/current_metadata.py').read_text()[:7000])
