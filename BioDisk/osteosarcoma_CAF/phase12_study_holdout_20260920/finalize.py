from pathlib import Path
import json,hashlib,shutil,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase12_study_holdout_20260920'
assert int((O/'exit_code').read_text())==0 and (O/'ANALYSIS_COMPLETE').exists()
assert json.loads((O/'REVIEW_QA.json').read_text())['all_checks_passed']
m=json.loads((O/'DELIVERY_SHA256.json').read_text())
for n,h in m.items():assert hashlib.sha256((O/n).read_bytes()).hexdigest()==h,n
s=json.loads((O/'STATUS.json').read_text());s.update(status='STUDY_HOLDOUT_REANALYSIS_COMPLETE_REVIEWED',validated_subtypes=0,training_candidate_pairs=1,heldout_primary_evaluable_candidate_records=0,report=str(O/'留出研究结果.md'),visual_QA='PASS')
(O/'STATUS.json').write_text(json.dumps(s,ensure_ascii=False,indent=2))
backup=P/'PROJECT_STATE.before_phase12_20260920.json'
if not backup.exists():shutil.copy2(P/'PROJECT_STATE.json',backup)
p=json.loads((P/'PROJECT_STATE.json').read_text());p.update(current_stage=O.name,status=s['status'],current_status=s['status'],dynamic_status_file=str(O/'STATUS.json'),latest_report=s['report'],updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase12_study_holdout_20260920=s,next_action='Prioritize continuous patient-level programs and technology-stratified validation; new matched cohorts needed before subtype claims; preserve heldout negative/unevaluable results')
(P/'PROJECT_STATE.json').write_text(json.dumps(p,ensure_ascii=False,indent=2));(O/'REVIEW_COMPLETE').write_text('STUDY_HOLDOUT_REVIEW_COMPLETE\n')
print(json.dumps({'status':s['status'],'verified_delivery_files':len(m)}))
