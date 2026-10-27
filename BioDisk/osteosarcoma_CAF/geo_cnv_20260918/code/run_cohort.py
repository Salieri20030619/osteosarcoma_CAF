from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,sys,json,csv,subprocess,datetime,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,obj):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n');t.replace(p)
def main():
 lock=(P/'logs/cohort.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert (P/'INPUTS_READY').exists() and (P/'EXPORT_VALIDATED').exists() and (P/'PREFLIGHT_PASS').exists()
 gate=json.loads((P/'review/annotation_gate.json').read_text());assert gate['cnv_launch_allowed'] is True
 rows=list(csv.DictReader((P/'review/input_summary.tsv').open(),delimiter='\t'));rows=[r for r in rows if r['cnv_eligible']=='True']
 rows.sort(key=lambda r:(r['site']!='lung',int(r['run_cells'])))
 assert rows
 state={'status':'DUAL_CNV_RUNNING','started_at_utc':now(),'pid':os.getpid(),'eligible_samples':len(rows),'planned_jobs':2*len(rows),'methods':['infercnv','copykat'],'parallel_method_lanes':2,'threads_per_job':4,'virtual_memory_limit_per_job_GiB':64,'normal_reference':'same-sample filtered immune/endothelial; no mesenchymal reference','minimum_training_normals':30,'minimum_heldout_calibration_normals':20,'all_retained_cells_assessed':True,'completion_policy':'native subprocess waiting; validate exit codes, markers, parameters and cell-result integrity','biological_interpretation':'no automatic CAF confirmation; FLEX and reference calibration review required'}
 dump(P/'STATUS.json',state);(P/'logs/cohort.pid').write_text(str(os.getpid())+'\n')
 def lane(method):
  records=[]
  for row in rows:
   sm=row['sample'];print('JOB_START',method,sm,now(),flush=True)
   with (P/'logs'/(method+'_'+sm+'.log')).open('a') as f:
    rc=subprocess.run(['bash',str(P/'code/run_cnv_job.sh'),method,sm],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT).returncode
   d=P/method/sm;passed=False;why='';calls={};held={}
   try:
    assert rc==0 and (d/'COMPLETE').exists() and (d/'job.exit_code').read_text().strip()=='0'
    params=json.loads((d/'parameters.json').read_text());assert params['reference_cells']>=30 and params['heldout_cells']>=20
    raw=list(csv.DictReader((P/'inputs'/sm/'cells.tsv').open(),delimiter='\t'));out=list(csv.DictReader((d/'cell_results.tsv').open(),delimiter='\t'))
    assert len(out)==len(raw)==int(row['run_cells'])
    assert [r['tool_id'] for r in out]==[r['tool_id'] for r in raw] and [r['cell_id'] for r in out]==[r['cell_id'] for r in raw]
    for rr in out:
     calls[rr['call']]=calls.get(rr['call'],0)+1
     if rr['cnv_role'].startswith('holdout_'):held[rr['call']]=held.get(rr['call'],0)+1
    passed=True
   except Exception as e:why=repr(e)
   rec={'method':method,'sample':sm,'exit_code':rc,'complete_marker':(d/'COMPLETE').exists(),'output_integrity_pass':passed,'integrity_error':why,'calls':calls,'heldout_calibration_calls':held,'finished_utc':now()}
   records.append(rec);dump(P/'logs'/(method+'_jobs.json'),records);print('JOB_EXIT',method,sm,rc,passed,flush=True)
  return records
 with ThreadPoolExecutor(max_workers=2) as pool:flat=sum(list(pool.map(lane,['infercnv','copykat'])),[])
 dump(P/'logs/cohort_jobs.json',flat)
 good=all(r['output_integrity_pass'] for r in flat);state.update(status='COMPUTATION_COMPLETE_REVIEW_REQUIRED' if good else 'COMPUTATION_FINISHED_WITH_FAILURES',finished_at_utc=now(),successful_jobs=sum(r['output_integrity_pass'] for r in flat),failed_jobs=sum(not r['output_integrity_pass'] for r in flat))
 dump(P/'STATUS.json',state);(P/'logs/cohort.exit_code').write_text('0\n' if good else '2\n');(P/'JOBS_FINISHED').write_text(state['status']+'\n')
 f=P.parent/'PROJECT_STATE.json';project=json.loads(f.read_text());project['geo_cnv_20260918']=state
 if project.get('current_stage')=='geo_cnv_20260918':project['current_status']=state['status'];project['status']=state['status'];project['updated_at_utc']=now()
 dump(f,project);print('GEO_DUAL_CNV_JOBS_FINISHED',state['status'],flush=True)
 return 0 if good else 2
try:sys.exit(main())
except Exception:
 (P/'logs/cohort.exit_code').write_text('1\n');dump(P/'COHORT_FAILED.json',{'time':now(),'traceback':traceback.format_exc()});raise
