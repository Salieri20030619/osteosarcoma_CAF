from pathlib import Path
import json,datetime,shutil,os,csv
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase23_loose_first_20260921'
now=datetime.datetime.now(datetime.timezone.utc);stamp=now.strftime('%Y%m%d_%H%M%S')
processes=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  cmd=(p/'cmdline').read_bytes().replace(b'\x00',b' ').decode(errors='replace');cwd=Path(os.readlink(p/'cwd'))
  if (str(O) in cmd or cwd==O) and any(s in cmd for s in ['prepare.py','cnv_official_inputs.py','cnv_official.R','model.py']):processes.append({'pid':int(p.name),'cwd':str(cwd),'command':cmd})
 except (OSError,ProcessLookupError):pass
jobs=O/'CNV_official_inputs/job_status.tsv'
jobrows=list(csv.DictReader(jobs.open(),delimiter='\t')) if jobs.exists() else []
snapshot={'recorded_utc':now.isoformat(),'conversation_closed_by_user':True,'handoff':str(P/'CAF_HANDOFF_20260921.md'),'completed_input_records':len(list((O/'inputs').glob('*.json'))),'inputs_complete_marker':(O/'INPUTS_COMPLETE').exists(),'background_processes':processes,'CNV_recorded_job_outcomes':jobrows,'CNV_run_finished':(O/'CNV_official_inputs/RUN_FINISHED').exists(),'all_cell_CAF_main_model_cancelled':(O/'logs/MODEL_CANCELLED_BY_SCOPE_CORRECTION.json').exists(),'corrected_CAF_and_myCAF_subclustering_complete':False,'delivered_15_page_report_is_historical_and_needs_correction':True,'latest_unresolved_message_4_is_not_an_instruction':True,'source_of_method':'client original requirements, not prior custom gate scripts'}
(O/'HANDOFF_STATE.json').write_text(json.dumps(snapshot,ensure_ascii=False,indent=2))
state=P/'PROJECT_STATE.json';backup=P/f'PROJECT_STATE.before_conversation_close_{stamp}.json';shutil.copy2(state,backup);v=json.loads(state.read_text())
v['phase23_loose_first_20260921']={'status':'IN_PROGRESS_HANDOFF_CLIENT_METHOD_CORRECTED','handoff_state':str(O/'HANDOFF_STATE.json'),'snapshot':snapshot}
v['current_user_priority']='Conversation ended at user request. On resumption read CAF_HANDOFF_20260921.md and original customer requirements; CNV identity screening -> CAF subtype -> myCAF subclustering. Do not resume cancelled all-cell CAF main model.'
v['latest_handoff']=str(P/'CAF_HANDOFF_20260921.md');tmp=state.with_suffix('.json.tmp');tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2));tmp.replace(state)
print(json.dumps({'handoff_saved':True,'project_state_backup':str(backup),'active_processes':len(processes),'CNV_recorded_outcomes':len(jobrows),'inputs_complete':snapshot['inputs_complete_marker']}))
