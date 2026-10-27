from pathlib import Path
import json,hashlib,tarfile,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919'
with tarfile.open(P/'phase10_batch_gates_20260919_final_updates.tar.gz',encoding='utf-8') as t:
 for m in t.getmembers():assert (O/m.name).resolve().is_relative_to(O.resolve()) and not m.issym() and not m.islnk()
 t.extractall(O)
man=json.loads((O/'MANIFEST_SHA256.json').read_text())
for name,h in man.items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==h,name
status=json.loads((O/'STATUS.json').read_text());assert status['status']=='BATCH_AND_EXCLUSION_RULES_IMPLEMENTED_REVIEWED'
f=P/'PROJECT_STATE.json';state=json.loads(f.read_text());state.update(current_stage=O.name,status=status['status'],current_status=status['status'],dynamic_status_file=str(O/'STATUS.json'),latest_report=str(O/status['report']),updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),next_action='Continue exact-name SR1004 work using unified raw-expression evidence and measurement coverage; preserve sensitivity sets and technology dependence')
state['batch_gate_remediation_20260919']={'status':status['status'],'sample_records_reviewed':40,'retained_cells_reviewed':255241,'primary_cells':3976,'old_joint_retained':3937,'new_vs_old':39,'old_not_in_new_strict':363,'main_clusters':11,'batch_models':34,'old_results_preserved':True,'technology_dependence_unresolved':True,'validated_biological_subtypes':False,'target_name':'SR1004','report':state['latest_report'],'delivery_files_verified':len(man)}
f.write_text(json.dumps(state,ensure_ascii=False,indent=2));print(json.dumps({'files_verified':len(man),'project_status':state['status'],'main_cells':3976,'report_pages':8},ensure_ascii=True))