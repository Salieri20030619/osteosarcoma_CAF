from pathlib import Path
import tarfile,json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase9_joint_20260919';archive=P/'phase9_joint_20260919_final_updates.tar.gz'
with tarfile.open(archive,encoding='utf-8') as tar:
 for item in tar.getmembers():
  dest=(O/item.name).resolve()
  assert dest.is_relative_to(O.resolve()) and not item.issym() and not item.islnk()
 tar.extractall(O)
m=json.loads((O/'MANIFEST_SHA256.json').read_text());bad=[]
for name,want in m.items():
 p=O/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=want:bad.append(name)
assert not bad,bad
assert (O/'logs/integration.exit_code').read_text().strip()=='0'
assert (O/'logs/review.exit_code').read_text().strip()=='0'
state=json.loads((P/'PROJECT_STATE.json').read_text());state.update(current_stage='phase9_joint_20260919',status='JOINT_INTEGRATION_REVIEW_COMPLETE',updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),dynamic_status_file=str(O/'STATUS.json'),next_action='Assess reference-paper methods using frozen joint integration; review competing identity programs before subtype naming.',latest_report=str(O/'delivery/骨肉瘤CAF_联合整合分析报告_20260919.pdf'))
state['joint_integration_20260919']={'status':'JOINT_INTEGRATION_REVIEW_COMPLETE','cells':4300,'datasets':5,'sample_records':17,'shared_genes':15391,'clusters':12,'batch_variable':'study_technology','source_counts_unchanged':True,'donor_dominated_clusters':[2,3,8,10,11],'within_study_paired_patients_range':[1,2],'inferential_tests_run':False,'report':state['latest_report'],'delivery_manifest_verified_files':len(m)}
(P/'PROJECT_STATE.json').write_text(json.dumps(state,ensure_ascii=False,indent=2))
print(json.dumps({'remote_manifest_verified_files':len(m),'project_state':state['status'],'report':state['latest_report']},ensure_ascii=True))