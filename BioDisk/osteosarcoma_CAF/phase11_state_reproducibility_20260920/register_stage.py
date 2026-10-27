from pathlib import Path
import json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase11_state_reproducibility_20260920'
assert int((O/'exit_code').read_text())==0 and int((O/'stress_exit_code').read_text())==0
assert (O/'ANALYSIS_COMPLETE').exists() and (O/'STRESS_SENSITIVITY_COMPLETE').exists()
assert json.loads((O/'REVIEW_QA.json').read_text())['all_checks_passed']
manifest=json.loads((O/'LOCAL_DELIVERY_SHA256.json').read_text())
for name,digest in manifest.items():
    path=O/name
    assert path.is_file(),name
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,name
state=json.loads((O/'STATUS.json').read_text())
state.update(status='STATE_REPRODUCIBILITY_ATTEMPT_COMPLETE_REVIEWED',targeted_refits=3,visual_QA='PASS',
 report=str(O/'阶段结果.md'),independent_validation=False,validated_subtypes=0,
 findings=['U0 stress-associated program consistent in four eligible records under frozen labels; structural persistence after stress exclusion; donor-specific effect reversal retained',
 'U3 cycle-associated candidate state; joins U2 after cycle exclusion',
 'U7 joins U6-dominated group after stress exclusion',
 'single-donor groups remain unsupported for cross-patient subtype claims'])
(O/'STATUS.json').write_text(json.dumps(state,ensure_ascii=False,indent=2))
backup=P/'PROJECT_STATE.before_phase11_20260920.json'
if not backup.exists():shutil.copy2(P/'PROJECT_STATE.json',backup)
project=json.loads((P/'PROJECT_STATE.json').read_text())
project.update(current_stage=O.name,dynamic_status_file=str(O/'STATUS.json'),
 status=state['status'],current_status=state['status'],updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 latest_report=str(O/'阶段结果.md'),next_action='Design training-isolated heldout donor/study validation of broader candidate programs; do not promote existing U labels to validated subtypes',
 SR1004='Deferred per user 2026-09-20; name issue acknowledged by customer',
 phase11_state_reproducibility_20260920=state)
(P/'PROJECT_STATE.json').write_text(json.dumps(project,ensure_ascii=False,indent=2))
(O/'REVIEW_COMPLETE').write_text('STATE_REPRODUCIBILITY_REVIEW_COMPLETE\n')
print(json.dumps({'status':state['status'],'verified_delivery_files':len(manifest),'original_stage_preserved':True},ensure_ascii=False))
