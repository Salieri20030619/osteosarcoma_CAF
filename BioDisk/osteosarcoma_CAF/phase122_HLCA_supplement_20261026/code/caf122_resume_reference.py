from pathlib import Path
import os,signal,json,subprocess,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026')
target=str(O/'code/caf122_census_probe.py');found=[]
for d in Path('/proc').iterdir():
 if not d.name.isdigit():continue
 try:args=(d/'cmdline').read_bytes().split(b'\0')
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if target.encode() in args:found.append(int(d.name))
assert len(found)<=1,found
for pid in found:os.kill(pid,signal.SIGTERM)
p=O/'code/caf122_reference_extract.py';s=p.read_text();s=s.replace("'CL:0019019']","'CL:0019019','CL:0000669']");p.write_text(s,encoding='utf-8')
log=(O/'logs/reference_launcher.log').open('a')
child=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python','-u',str(O/'code/caf122_run_stage.py'),'reference_extract','caf122_reference_extract.py'],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
rec=dict(stopped_own_probe=found,reason='Remote all-column probe exceeded25min; retry selected columns and ontology IDs',runner_pid=child.pid,time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'resources/REFERENCE_READ_RETRY.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec))
