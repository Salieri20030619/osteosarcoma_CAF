from pathlib import Path
import subprocess,json,concurrent.futures,shlex,time
ROOT=Path(__file__).resolve().parent
REMOTE='/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/CAF_access_recheck_20261013'
kw={'creationflags':subprocess.CREATE_NO_WINDOW,'capture_output':True,'text':True}
code="from pathlib import Path;import json;d=Path("+repr(REMOTE+'/downloads/NCC')+");print(json.dumps({p.name:p.stat().st_size for p in d.iterdir() if p.is_file()}))"
r=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=25','sgs10','python3 -c '+shlex.quote(code)],check=True,**kw)
sizes=json.loads(r.stdout);manifest=json.loads((ROOT/'NCC_DOWNLOAD_COMPLETE.json').read_text(encoding='utf8'))
pending=[f for f in manifest['files'] if sizes.get(f['file'])!=f['bytes']]
print('PENDING',[(f['file'],f['bytes']) for f in pending],flush=True)
def copy(f):
 start=time.time();r=subprocess.run(['scp','-q','-o','BatchMode=yes','-o','ConnectTimeout=25',str(ROOT/'downloads/NCC'/f['file']),'sgs10:'+REMOTE+'/downloads/NCC/'],**kw)
 row={'file':f['file'],'exit_code':r.returncode,'seconds':time.time()-start,'stderr':r.stderr};print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(copy,pending))
(ROOT/'PARALLEL_TRANSFER_RESULT.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
assert all(r['exit_code']==0 for r in rows)
print('TRANSFER_COMPLETE',flush=True)
