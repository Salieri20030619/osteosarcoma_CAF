from pathlib import Path
import json,hashlib,datetime,shutil
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase120_patient_supported_CAF_review_20261026';t=datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
q=json.loads((O/'resources/LOCAL_QA.json').read_text());assert q['html_sha256']==sha(O/'delivery/CAF_patient_support_literature_review.html');assert q['txt_sha256']==sha(O/'delivery/CAF_patient_support_literature_review.txt');assert q['browser_exit_codes']==[0,0]
q['visual_review']='PASS: main Chinese report and all-parent gate tables inspected, source tie-break correction verified';q['verified_at']=t;dump(O/'DELIVERY_RECEIPT.json',q)
a=json.loads((O/'ACCEPTANCE.json').read_text());a['status']='PASS_VERIFIED';a['remaining']=[];a['delivery_hashes_verified']=True;dump(O/'ACCEPTANCE.json',a)
final=dict(status='completed_verified',phase=O.name,completed_at=t,cells=6633,main_cells=6423,parents_reviewed=29,children_reviewed=36,final_candidates=0,priority_exploratory_parent='Lung_04',secondary_exploratory_parent='Bone_07',threshold=dict(evaluable_independent_donors=6,cells_each_side=5,positive_fraction=.75,preferred_donors=10),final_gate_note='Lung_04 passes main-source function and7 donor pairs; identity sensitivity5 donor pairs due control attrition; remains strongest exploratory signal',figures='paused',report=str(O/'delivery/CAF_patient_support_literature_review.html'),local_report='C:/Users/31937/Downloads/CAF_Phase120_patient_review_20261026/CAF_patient_support_literature_review.html',remaining_required_work=[])
dump(O/'FINAL_STATUS.json',final);dump(O/'STATUS.json',final)
active=json.loads((O/'ACTIVE_PROCESS.json').read_text());active['last_pid']=active.pop('pid');active.update(status='completed_verified',pid=None,completed_at=t);dump(O/'ACTIVE_PROCESS.json',active)
dump(O/'resources/EXECUTION_NOTES.json',dict(initial_launch='SCP succeeded; first inline remote launch failed at shell quoting before any analysis. Replaced by Python subprocess wrapper; no analytic changes.',main_native_exit=0,focus_tool_exit=0,report_tool_exit=0,report_revision='Actual source required in display tie-break; statistics untouched',scientific_plots=False))
ptr=R/'CAF_PATIENT_SUPPORTED_CANDIDATE_CURRENT.json';dump(ptr,final)
old=R/'CAF_RESUME_STATUS_CURRENT.json';shutil.copyfile(old,O/'resources/RESUME_before_Phase120_completion.json');x=json.loads(old.read_text());x.update(authoritative_current_file=str(ptr),current_phase=str(O),phase=O.name,status='completed_verified',current_state='completed_verified',active_task='Phase120 patient-supported candidate audit and literature crosscheck complete',latest_task=final,last_completed_task=final,updated_at=t,next_action='Review Lung_04 exploratory matrix evidence and additional patient/control requirements; no active computation; figures remain paused');dump(old,x)
with (O/'events.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(dict(utc=t,stage='final_verified',final_candidates=0,remaining=[]))+'\n')
files=[p for p in O.rglob('*') if p.is_file() and p.name!='PROCESS_RECORD_MANIFEST.json'];dump(O/'PROCESS_RECORD_MANIFEST.json',[dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=sha(p)) for p in files])
assert json.loads(ptr.read_text())['status']=='completed_verified';assert json.loads((O/'ACCEPTANCE.json').read_text())['remaining']==[]
print(json.dumps(final,ensure_ascii=False))
