from pathlib import Path
import os,signal,json,shutil
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');target=str(O/'code/caf122_gene_and_census.py');found=[]
for d in Path('/proc').iterdir():
 if not d.name.isdigit():continue
 try:args=(d/'cmdline').read_bytes().split(b'\0')
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if target.encode() in args:found.append(int(d.name))
assert len(found)<=1,found
for pid in found:os.kill(pid,signal.SIGTERM)
shutil.copyfile(target,O/'resources/gene_reconciliation_attempt1.py')
(O/'resources/GENE_RECONCILIATION_RUNTIME_FIX.json').write_text(json.dumps(dict(stopped_own_pids=found,reason='Compute duplicate symbol mask once instead of per feature; same mapping definition; avoids quadratic repeated computation',prior_query_preserved=True),indent=2),encoding='utf-8');print(found)
