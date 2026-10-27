from pathlib import Path
import os,json,time,datetime,fcntl,traceback
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase25_complete_nonPDX_assets_20260921');M=O/'monitoring';M.mkdir(exist_ok=True)
lock=(M/'monitor.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
(M/'monitor.pid').write_text(str(os.getpid())+'\n')
INTERVAL=900
TZ=datetime.timezone(datetime.timedelta(hours=8))
def read_json(p,default):
 for attempt in range(3):
  try:return json.loads(p.read_text()) if p.exists() else default
  except (json.JSONDecodeError,OSError):time.sleep(1)
 return default

def save(name,data):
 p=M/name;t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(data,ensure_ascii=False,indent=2));t.replace(p)

def active(pid,needle):
 try:return needle in (Path('/proc')/str(pid)/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
 except OSError:return False

def snapshot():
 s=read_json(O/'STATUS.json',{});files=read_json(O/'verified_files.json',[]);inv=read_json(O/'matrix_inventory.json',[]);plan=read_json(O/'selected_file_manifest.json',[])
 try:pid=int((O/'runner.pid').read_text().strip())
 except (ValueError,OSError):pid=None
 rc=None
 if (O/'exit_code').exists():
  try:rc=int((O/'exit_code').read_text().strip())
  except (ValueError,OSError):pass
 kinds={}
 for r in inv:
  kinds.setdefault(r['sample'],set())
  for k in ['filtered','processed','unfiltered']:
   if r['file'].endswith('_'+k+'_rna.h5ad'):kinds[r['sample']].add(k)
 total=sum(r['uncompressed_bytes'] for r in plan);done=sum(r['bytes'] for r in files)
 alive=active(pid,'run_official_download.py') if pid else False
 d={'recorded_at':datetime.datetime.now(TZ).isoformat(),'interval_seconds':INTERVAL,'runner_pid':pid,'runner_active':alive,'download_state':s.get('status','not_recorded'),'verified_complete_files':len(files),'verified_complete_bytes':done,'planned_bytes':total,'verified_percent':round(done/total*100,2) if total else None,'verified_matrix_objects':len(inv),'new_samples_with_all_three_objects':sum(len(v)==3 for v in kinds.values()),'new_samples_requested':55,'current_file':s.get('current_file'),'exit_code':rc,'run_finished_marker':(O/'RUN_FINISHED').exists(),'all67_marker':(O/'ALL67_ASSETS_VERIFIED').exists(),'bone_reference_complete':(O/'references/GSE217792/COMPLETE').exists(),'HLCA_reference_complete':(O/'references/HLCA/COMPLETE').exists()}
 if alive:d['assessment']='RUNNING'
 else:
  qa=read_json(O/'DOWNLOAD_QA.json',{});mf=O/'all67_nonPDX_asset_manifest.tsv';checks={'runner_exit_zero':rc==0,'run_finished_marker':d['run_finished_marker'],'all67_marker':d['all67_marker'],'QA_67_samples':qa.get('total_nonPDX_sample_records')==67,'QA_165_new_objects':qa.get('new_expression_objects')==165,'no_PDX_or_cell_lines':qa.get('PDX_in_manifest')==0 and qa.get('cell_lines_in_manifest')==0,'manifest_exists':mf.exists()}
  if mf.exists():
   a=pd.read_csv(mf,sep='\t',keep_default_na=False);checks['manifest_67_unique_samples']=len(a)==67 and a.scpca_sample_id.is_unique;checks['manifest_source_flags']=a.is_xenograft.astype(str).str.lower().eq('false').all() and a.is_cell_line.astype(str).str.lower().eq('false').all();checks['all_201_objects_present']=all(Path(p).is_file() and Path(p).stat().st_size>0 for c in ['filtered','processed','unfiltered'] for p in a[c])
  d['completion_checks']={k:bool(v) for k,v in checks.items()};d['assessment']='DOWNLOAD_COMPLETE_VERIFIED' if all(checks.values()) else 'STOPPED_REQUIRES_REVIEW'
 return d
while True:
 try:
  d=snapshot();save('LATEST.json',d)
  with (M/'progress.jsonl').open('a') as f:f.write(json.dumps(d,ensure_ascii=False)+'\n')
  with (M/'progress.log').open('a') as f:f.write(f"{d['recorded_at']} {d['assessment']} files={d['verified_complete_files']} verified_GB={d['verified_complete_bytes']/1e9:.3f}/{d['planned_bytes']/1e9:.3f} samples={d['new_samples_with_all_three_objects']}/55 exit={d['exit_code']}\n")
  if d['assessment']!='RUNNING':
   save('FINAL.json',d);(M/'monitor.exit_code').write_text('0\n' if d['assessment']=='DOWNLOAD_COMPLETE_VERIFIED' else '1\n');break
 except Exception as e:
  with (M/'monitor_errors.log').open('a') as f:f.write(datetime.datetime.now(TZ).isoformat()+' '+type(e).__name__+'\n')
 time.sleep(INTERVAL)
