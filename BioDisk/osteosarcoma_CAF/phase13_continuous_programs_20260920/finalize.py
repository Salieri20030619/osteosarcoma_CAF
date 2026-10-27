from pathlib import Path
import json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase13_continuous_programs_20260920'
assert int((O/'exit_code').read_text())==0
assert (O/'ANALYSIS_COMPLETE').exists() and (O/'SPLIT_SENSITIVITY_COMPLETE').exists()
assert json.loads((O/'REVIEW_QA.json').read_text())['all_checks_passed']
m=json.loads((O/'DELIVERY_SHA256.json').read_text())
for n,h in m.items():assert hashlib.sha256((O/n).read_bytes()).hexdigest()==h,n
s=json.loads((O/'STATUS.json').read_text());s.update(status='CONTINUOUS_PROGRAM_COHERENCE_COMPLETE_REVIEWED',report=str(O/'连续程序结果.md'),visual_QA='PASS',main_findings=['Stress coherence positive across 11 evaluable records and all 25 alternative partitions','Matrix coherence strong in 7 non-FLEX records; weaker and sensitivity-dependent in FLEX','Inflammatory panel has weak and inconsistent joint coherence'],validated_subtypes=0)
(O/'STATUS.json').write_text(json.dumps(s,ensure_ascii=False,indent=2))
backup=P/'PROJECT_STATE.before_phase13_20260920.json'
if not backup.exists():shutil.copy2(P/'PROJECT_STATE.json',backup)
p=json.loads((P/'PROJECT_STATE.json').read_text());p.update(current_stage=O.name,status=s['status'],current_status=s['status'],dynamic_status_file=str(O/'STATUS.json'),latest_report=s['report'],updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase13_continuous_programs_20260920=s,next_action='Freeze continuous stress and matrix program scores for independent clinical/spatial context validation; keep technology strata and inflammatory incoherence; no validated subtype claims')
(P/'PROJECT_STATE.json').write_text(json.dumps(p,ensure_ascii=False,indent=2));(O/'REVIEW_COMPLETE').write_text('CONTINUOUS_PROGRAM_REVIEW_COMPLETE\n');print(json.dumps({'status':s['status'],'verified_files':len(m)}))
