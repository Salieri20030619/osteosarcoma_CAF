from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,gzip,hashlib,datetime,threading
p=Path.cwd();root=p/'download_expansion_20260917';rows=json.loads((p/'candidate_sources_20260917/GEO_candidate_sample_manifest.json').read_text());jobs=[]
for row in rows:
 if row['dataset'] in ['GSE318840','GSE162454','GSE250015']:
  for url in row['supplementary_files']:jobs.append({'dataset':row['dataset'],'sample':row['sample'],'url':url.replace('ftp://','https://')})
assert len({x['sample'] for x in jobs})==16
(root/'geo_download_manifest.json').write_text(json.dumps(jobs,indent=2));lock=threading.Lock()
def run(j):
 dest=root/'files'/j['dataset']/j['url'].split('/')[-1];dest.parent.mkdir(exist_ok=True);part=Path(str(dest)+'.part');result=dict(j)
 if dest.exists():result['status']='already_downloaded'
 else:
  cmd=['curl','-4','-f','-L','--retry','3','--retry-delay','20','--connect-timeout','30','--max-time','7200','--speed-limit','100','--speed-time','180','-C','-','-o',str(part),j['url']]
  with (root/'logs'/(dest.name+'.log')).open('a') as log:rc=subprocess.run(cmd,stdout=log,stderr=log).returncode
  if rc:result.update(status='download_failed',exit_code=rc,partial_bytes=part.stat().st_size if part.exists() else 0);return result
  try:
   if part.name.endswith('.gz.part'):
    with gzip.open(part,'rb') as f:
     while f.read(1024*1024):pass
   part.replace(dest);result['status']='complete'
  except Exception as e:result.update(status='integrity_failed',error=str(e));return result
 h=hashlib.sha256()
 with dest.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 result.update(bytes=dest.stat().st_size,sha256=h.hexdigest());return result
results=[]
with ThreadPoolExecutor(2) as ex:
 for result in ex.map(run,jobs):
  results.append(result);(root/'geo_download_progress.json').write_text(json.dumps(results,indent=2));print(result['dataset'],result['sample'],result['status'],flush=True)
state={'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'total_files':len(jobs),'completed_files':sum(x['status'] in ['complete','already_downloaded'] for x in results),'results':results,'next_action':'STOP_AFTER_DOWNLOAD; no QC CNV integration or report changes authorized tonight'}
(root/'geo_download_completion.json').write_text(json.dumps(state,indent=2));print('GEO_DOWNLOAD_ATTEMPT_FINISHED',state['completed_files'],len(jobs),flush=True)
