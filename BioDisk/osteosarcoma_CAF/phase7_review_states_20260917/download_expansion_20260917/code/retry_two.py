import pathlib,subprocess,json,datetime,os,gzip,hashlib,concurrent.futures,sys,fcntl
r=pathlib.Path(__file__).resolve().parents[1]
lock=(r/'logs/retry_two.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
os.environ.update(http_proxy='http://127.0.0.1:7897',https_proxy='http://127.0.0.1:7897',HTTP_PROXY='http://127.0.0.1:7897',HTTPS_PROXY='http://127.0.0.1:7897')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(path,data):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2));tmp.replace(path)
old=json.loads((r/'geo_download_completion.json').read_text())
jobs=[dict(x) for x in old['results'] if x['sample'] in ['GSM5155199','GSM5155200'] and x['status']=='download_failed']
assert len(jobs)==2
def run(j):
 dest=r/'files'/j['dataset']/j['url'].split('/')[-1];part=pathlib.Path(str(dest)+'.part')
 result={k:j[k] for k in ['dataset','sample','url']}
 cmd=['curl','-4','--http1.1','-f','-L','-sS','--retry','5','--retry-all-errors','--retry-delay','10','--connect-timeout','30','--max-time','3600','--speed-limit','1024','--speed-time','120','-C','-','-o',str(part),j['url']]
 print('START',j['sample'],now(),flush=True)
 with (r/'logs'/(j['sample']+'_retry_two.curl.log')).open('a') as log:
  rc=subprocess.run(cmd,stdout=log,stderr=log).returncode
 result['curl_exit_code']=rc
 if rc:
  result.update(status='download_failed',exit_code=rc,partial_bytes=part.stat().st_size if part.exists() else 0)
 else:
  try:
   with gzip.open(part,'rb') as f:
    while f.read(4*1024*1024):pass
   h=hashlib.sha256()
   with part.open('rb') as f:
    for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
   result.update(status='complete',bytes=part.stat().st_size,sha256=h.hexdigest(),gzip_integrity='PASS')
   part.replace(dest)
  except Exception as e:result.update(status='integrity_failed',error=str(e))
 save(r/'logs'/(j['sample']+'_retry_two.result.json'),result)
 print('FINISHED',j['sample'],result['status'],now(),flush=True)
 return result
try:
 with concurrent.futures.ThreadPoolExecutor(2) as pool:results=list(pool.map(run,jobs))
 mapping={x['url']:x for x in results}
 old['results']=[mapping.get(x['url'],x) for x in old['results']]
 old['finished_utc']=now();old['completed_files']=sum(x['status'] in ['complete','already_downloaded'] for x in old['results'])
 save(r/'geo_download_completion.json',old);save(r/'geo_download_progress.json',old['results'])
 ok=old['completed_files']==old['total_files']==44
 completion={'process_exit_code':0 if ok else 2,'all_files_complete':ok,'total_files':old['total_files'],'completed_files':old['completed_files'],'finished_utc':now(),'retried_samples':[x['sample'] for x in results]}
 save(r/'logs/retry_two_completion.json',completion)
 save(r/'logs/geo_proxy_completion.json',completion)
 s=json.loads((r/'STATUS.json').read_text());s['state']='PAUSED_AFTER_DOWNLOAD_ATTEMPTS';s['GEO_completion']=completion;s['finished_utc']=now()
 save(r/'STATUS.json',s)
 p=r.parent.parent/'PROJECT_STATE.json';z=json.loads(p.read_text());z['current_status']=s['state'];z['download_expansion']=s;save(p,z)
 code=0 if ok else 2
 print('RETRY_TWO_FINISHED',json.dumps(completion),flush=True)
except BaseException as e:
 code=1;print('RETRY_TWO_ERROR',repr(e),flush=True)
(r/'logs/retry_two.exit_code').write_text(str(code)+'\n')
sys.exit(code)
