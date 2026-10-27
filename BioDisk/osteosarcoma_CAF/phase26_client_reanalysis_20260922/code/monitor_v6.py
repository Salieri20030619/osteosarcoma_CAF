from pathlib import Path
import os,json,datetime,select,ctypes
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase26_client_reanalysis_20260922')
M=O/'monitoring'
pid=int((O/'runner.pid').read_text())
(M/'monitor.pid').write_text(str(os.getpid())+'\n')
def now():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def write(p,x):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2));t.replace(p)
def resource():
 processes={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   f=(p/'stat').read_text().rsplit(')',1)[1].split()
   if (p/'status').stat().st_uid!=os.getuid():continue
   processes[int(p.name)]={'pid':int(p.name),'ppid':int(f[1]),'state':f[0],'rss_GiB':int(f[21])*os.sysconf('SC_PAGE_SIZE')/2**30,'cpu_seconds':(int(f[11])+int(f[12]))/os.sysconf('SC_CLK_TCK'),'command':(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')}
  except (OSError,ValueError,IndexError):continue
 roots={pid}
 h=json.loads((O/'RESERVE50_HANDOVER.json').read_text())
 for c in h['controllers']:
  old=processes.get(c['pid'])
  if old and c['script'] in old['command']:roots.add(c['pid'])
 selected=set(roots)
 while True:
  added={k for k,v in processes.items() if v['ppid'] in selected}-selected
  if not added:break
  selected.update(added)
 rows=[processes[k] for k in sorted(selected) if k in processes]
 for row in rows:
  try:
   lines=(Path('/proc')/str(row['pid'])/'smaps_rollup').read_text().splitlines()
   row['pss_GiB']=sum(int(l.split()[1]) for l in lines if l.startswith('Pss:'))/2**20
  except OSError:row['pss_GiB']=None
 mem={l.split(':')[0]:int(l.split()[1]) for l in Path('/proc/meminfo').read_text().splitlines()}
 return {'available_memory_GiB':mem['MemAvailable']/2**20,'total_memory_GiB':mem['MemTotal']/2**20,'process_RSS_sum_GiB_shared_pages_may_repeat':sum(r['rss_GiB'] for r in rows),'process_PSS_sum_GiB':sum(r['pss_GiB'] or 0 for r in rows),'load_average':list(os.getloadavg()),'processes':rows}
libc=ctypes.CDLL(None,use_errno=True)
fd=libc.syscall(434,pid,0)
if fd<0:raise OSError(ctypes.get_errno(),'pidfd_open')
while True:
 alive=(Path('/proc')/str(pid)/'cmdline').exists()
 state=json.loads((O/'STATUS.json').read_text())
 rec={'recorded_at':now(),'runner_pid':pid,'runner_present':alive,'interval_seconds':900,'state':state,'resources':resource(),'exit_code':(O/'exit_code').read_text().strip() if (O/'exit_code').exists() else None,'stage_finished_marker':(O/'CNV_STAGE_FINISHED').exists(),'stage_QA':json.loads((O/'delivery/CNV_STAGE_QA.json').read_text()) if (O/'delivery/CNV_STAGE_QA.json').exists() else None}
 rec['memory_admission']=json.loads((O/'MEMORY_GUARD_STATE.json').read_text()) if (O/'MEMORY_GUARD_STATE.json').exists() else None
 write(M/'LATEST.json',rec)
 with (M/'progress.jsonl').open('a') as f:f.write(json.dumps(rec,ensure_ascii=False)+'\n')
 if not alive or (O/'RUN_FINISHED').exists():
  write(M/'FINAL.json',rec);break
 select.select([fd],[],[],900)
os.close(fd)
