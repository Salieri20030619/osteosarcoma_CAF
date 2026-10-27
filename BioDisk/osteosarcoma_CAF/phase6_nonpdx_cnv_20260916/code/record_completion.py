"""SYNTAX_CHECK: python -m py_compile code/record_completion.py
RUN_COMMAND: python code/record_completion.py
EXPECTED_SUCCESS_MARKER: PHASE6_COMPLETION_STATE_RECORDED
"""
from pathlib import Path
import json,datetime
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase6_nonpdx_cnv_20260916')
state=json.loads((p/'STATUS.json').read_text());ec=p/'logs/cohort.exit_code'
state['cohort_exit_code']=int(ec.read_text()) if ec.exists() else None
state['cohort_process_running']=False;state['completion_observed_at_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
if state.get('state')=='DUAL_CNV_RUNNING':state['state']='COHORT_EXITED_REVIEW_REQUIRED'
state['next_action']='review actual CNV/reference consistency, failures and candidate sets before advancing biological analysis'
(p/'STATUS.json').write_text(json.dumps(state,indent=2))
project=p.parent/'PROJECT_STATE.json';old=json.loads(project.read_text())
if old.get('current_stage')==p.name:
 old['phase6_completion_observed_at_utc']=state['completion_observed_at_utc'];old['phase6_cohort_exit_code']=state['cohort_exit_code'];old['next_model_review']=state['next_action'];project.write_text(json.dumps(old,ensure_ascii=False,indent=2))
print('PHASE6_COMPLETION_STATE_RECORDED',json.dumps(state),flush=True)
