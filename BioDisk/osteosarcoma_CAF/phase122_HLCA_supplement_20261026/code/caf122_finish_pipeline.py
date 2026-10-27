from pathlib import Path
import os,signal,subprocess,json,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python'
targets=[str(O/'code'/s).encode() for s in ['caf122_h5_extract.py','caf122_reference_extract.py']];stopped=[]
for d in Path('/proc').iterdir():
 if not d.name.isdigit():continue
 try:args=(d/'cmdline').read_bytes().split(b'\0')
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if any(t in args for t in targets):os.kill(int(d.name),signal.SIGTERM);stopped.append(int(d.name))
(O/'resources/REFERENCE_LOCAL_READ_SWITCH.json').write_text(json.dumps(dict(stopped_own_readers=stopped,reason='Random remote HDF5 reads only600rows/225sec; preserve checkpoints and fetch4 verified byte ranges in parallel, then read locally',reference_sampling_cap=100),indent=2),encoding='utf-8')
for stage,script,marker in [('reference_download','caf122_reference_download.py','resources/REFERENCE_FULL_DOWNLOAD_COMPLETE.json'),('h5_extract','caf122_h5_extract.py','REFERENCE_HTTP_COMPLETE.json'),('supplement_analysis','caf122_analyze_supplement.py','ANALYSIS_COMPLETE.json')]:
 if (O/marker).exists():continue
 rc=subprocess.call([py,'-u',str(O/'code/caf122_run_stage.py'),stage,script]);assert rc==0,(stage,rc);assert (O/marker).exists()
print('PIPELINE_COMPLETE',flush=True)
