
from pathlib import Path
import os,json,datetime,selectors,fcntl,subprocess,hashlib,traceback,time
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase34_reference_rerun_20260927';S=P/'phase26_client_reanalysis_20260922/resource_recovery_SCPCS000525_20260925'
lock=(O/'supervision.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now():return datetime.datetime.now().astimezone().isoformat()
def read(f):
 try:return json.loads(f.read_text())
 except Exception:return None
def proc(pid):
 try:
  t=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split();return dict(pid=pid,state=t[0],ppid=int(t[1]),ticks=int(t[11])+int(t[12]),start_ticks=int(t[19]),rss_pages=int(t[21]))
 except (FileNotFoundError,ProcessLookupError,PermissionError):return None
def tree(pid):
 rows={}
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   r=proc(int(p.name))
   if r:rows[r['pid']]=r
 ids={pid}
 while True:
  found={k for k,r in rows.items() if r['ppid'] in ids}-ids
  if not found:break
  ids|=found
 vals=[rows[k] for k in ids if k in rows];return dict(processes=vals,cpu_ticks_total=sum(x['ticks'] for x in vals),rss_pages_total=sum(x['rss_pages'] for x in vals))
def tail(p):
 if not p.exists():return None
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-700));return dict(bytes=p.stat().st_size,mtime=p.stat().st_mtime,tail=f.read().decode(errors='replace'))
sel=selectors.DefaultSelector();watched={};audits={};done={};previous={};quiet={}
roots={'rerun':int((O/'controller.pid').read_text()),'recovery':int((S/'runner.pid').read_text())}
def watch(key,pid):
 if key in watched:return
 try:
  fd=os.pidfd_open(pid);sel.register(fd,selectors.EVENT_READ,(key,pid));watched[key]=fd
 except ProcessLookupError:watched[key]=None
def audit(kind):
 if kind in audits or kind in done:return
 log=(O/'logs'/f'supervision_{kind}_audit.log').open('x')
 p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/supervision_acceptance.py'),kind],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT);log.close();audits[kind]=p;watch('audit:'+kind,p.pid)
def snapshot(event):
 status=read(O/'STATUS.json');recovery=read(S/'STATUS.json');ts={k:tree(pid) for k,pid in roots.items() if k not in done}
 logs={f.name:tail(f) for f in (O/'logs').glob('SCPCS*.log')};logs['S525']=tail(S/'logs/copykat_SCPCS000525.log')
 notes=[]
 for k,t in ts.items():
  marker=(t['cpu_ticks_total'],tuple((n,v['bytes']) for n,v in logs.items() if v and (n=='S525')==(k=='recovery')))
  if event=='scheduled_snapshot':quiet[k]=quiet.get(k,0)+1 if previous.get(k)==marker else 0
  previous[k]=marker
  if quiet.get(k,0)>=3:notes.append(k+': CPU ticks and log sizes unchanged across 3 monitoring intervals; review required, no automatic kill/restart.')
 item=dict(at=now(),event=event,monitor_pid=os.getpid(),interval_seconds=900,native_exit_events=True,rerun=status,recovery=recovery,process_trees=ts,logs=logs,review_notes=notes,audit_outcomes=done)
 f=O/'SUPERVISION_STATUS.json';tmp=f.with_suffix('.tmp');tmp.write_text(json.dumps(item,indent=2));tmp.replace(f)
 with (O/'logs/supervision_history.jsonl').open('a') as f:f.write(json.dumps(item)+'\n')
for k,pid in roots.items():watch(k,pid)
state=read(O/'STATUS.json') or {}
for name,pid in state.get('active',{}).items():watch('native:'+name,pid)
snapshot('supervision_started')
for k in roots:
 if watched.get(k) is None:audit(k)
while len(done)<2:
 events=sel.select(timeout=900)
 if not events:snapshot('scheduled_snapshot');continue
 for key,_ in events:
  name,pid=key.data;sel.unregister(key.fd);os.close(key.fd)
  if name in roots:audit(name)
  elif name.startswith('audit:'):
   kind=name.split(':',1)[1];rc=audits[kind].wait();done[kind]=dict(exit_code=rc,acceptance=read(O/'delivery'/f'SUPERVISION_{kind.upper()}_ACCEPTANCE.json'))
  snapshot('process_exit:'+name)
(O/'SUPERVISION_FINISHED').write_text(now());snapshot('all_watched_jobs_audited')
