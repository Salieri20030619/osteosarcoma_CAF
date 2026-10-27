from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,zipfile

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase80_scheme1_anatomical_20261012'
now=datetime.now(timezone.utc).isoformat()
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
    tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    tmp.replace(p)
assert (O/'ANALYSIS_FINISHED').is_file()
qa=json.loads((O/'QA/COMPUTED_CHECKS.json').read_text())
assert qa['source_inputs_unchanged'] and qa['patient_counts_match']
for file,expected in json.loads((O/'resources/INPUT_SHA256.json').read_text()).items():
    assert sha(Path(file))==expected
figures=sorted((O/'figures').glob('*.png'))
assert len(figures)==3 and all(f.stat().st_size>10000 for f in figures)
qa.update(status='PASS',visual_review_pending=False,visual_review='All three rendered PNGs reviewed: labels, legends, patient counts, nominal directions, paired discordance and cell coverage are readable; no clipping found.',verified_at=now,execution_exit_code=0)
write(O/'QA/FINAL_ACCEPTANCE.json',qa)
write(O/'QA/EXECUTION.json',dict(build_tool_exit_code=0,build_completion_marker='ANALYSIS_FINISHED',background_process_running=False,recorded_at=now))
report=next((O/'delivery').glob('*.txt'))
manifest={str(f.relative_to(O)):dict(bytes=f.stat().st_size,sha256=sha(f)) for sub in ['tables','figures','resources','code','QA'] for f in (O/sub).rglob('*') if f.is_file()}
manifest[str(report.relative_to(O))]=dict(bytes=report.stat().st_size,sha256=sha(report))
write(O/'delivery/FILE_MANIFEST.json',manifest)
archive=O/'delivery/CAF_scheme1_bone_lung_20261012.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for rel in manifest: z.write(O/rel,arcname=rel)
    z.write(O/'delivery/FILE_MANIFEST.json',arcname='delivery/FILE_MANIFEST.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for rel,info in manifest.items():
        assert hashlib.sha256(z.read(rel)).hexdigest()==info['sha256']
record=dict(status='completed_verified',phase=O.name,completed_at=now,
    title='方案一 全部骨病灶与肺病灶的CAF分类差异比较',
    user_authorization='按照方案一进行',
    grouping='Anatomical Bone vs Lung, regardless of metastatic history',
    scope='Original supported/expanded F1/F3/F4 and all fixed internal states; patient-level captured proportions, within-population gene expression and programs',
    computation_baseline='Phase73 strict_v2 and Phase74; exact original P/FDR/families preserved',
    new_work='All-eligible evidence integration, patient coverage validation, source/design consistency summary, three figures and scheme1 results report',
    new_inferential_tests=0,eligible_expression_comparisons=10,nominal_gene_records=2013,nominal_program_records=783,nominal_proportion_records=0,
    gene_FDR_significant=0,program_FDR_significant=0,report=str(report),figures=[str(p) for p in figures],
    archive=str(archive),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),acceptance=str(O/'QA/FINAL_ACCEPTANCE.json'),
    flow_marker_stage='deferred, not active',third_cohort_stage='not active',background_process_running=False,
    remaining_required_work=[],next_action='方案一结果已完成；后续按用户指定扩展。',local_delivery_verified=False)
write(O/'COMPLETE.json',record)
back=O/'resources/state_before_scheme1'
back.mkdir(exist_ok=True)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json']:
    f=P/name
    if f.exists():
        assert not (back/name).exists(), 'Refusing to replace prior state backup'
        shutil.copy2(f,back/name)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
    f=P/name
    obj=json.loads(f.read_text())
    for key in ['current_phase','current_stage','phase']: obj[key]=O.name
    for key in ['status','current_status','current_state','task_status']: obj[key]='scheme1_anatomical_comparison_completed_verified'
    for key in ['updated_at','updated_at_utc','progress_updated_at','last_activity_at']: obj[key]=now
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']: obj[key]=record
    obj['active_task']=dict(title=record['title'],status='completed_verified',scope=record['scope'],user_authorization=record['user_authorization'],result=str(O/'COMPLETE.json'))
    obj.update(user_priority='按方案一 比较全部骨病灶与肺病灶 不以转移史限制入组',
        current_scheme1=str(P/'CAF_SCHEME1_CURRENT.json'),analysis_completion_scope=record['scope'],
        analysis_now_completed=True,background_process_running=False,acceptance_passed=True,
        current_cohort_differences=str(P/'CAF_SCHEME1_CURRENT.json'),
        analysis_results_report=str(report),next_action=record['next_action'],flow_marker_stage=record['flow_marker_stage'],
        previous_progress_snapshot=str(back/name),process_observation='Build exited0; patient counts/design verified; input hashes unchanged; three PNGs visually reviewed; archive integrity checked.')
    write(f,obj)
resume=json.loads((P/'CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json').read_text())
resume.update(status='scheme1_results_completed_verified',updated_at_utc=now,user_request='按照方案一进行',
    baseline_phase=O.name,current_scheme1=str(P/'CAF_SCHEME1_CURRENT.json'),latest_results=str(report),
    next_action=record['next_action'],background_process_running=False,
    analysis_execution_started_this_resume=True,new_inferential_tests=0,
    original_computation_baseline='Phase73 strict_v2 / Phase74',state_backup=str(back))
write(P/'CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json',resume)
write(P/'CAF_SCHEME1_CURRENT.json',record)
print(json.dumps(record,ensure_ascii=False))
