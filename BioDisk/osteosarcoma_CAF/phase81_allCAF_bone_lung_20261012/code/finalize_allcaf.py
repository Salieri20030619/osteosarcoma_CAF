from pathlib import Path
import json,hashlib,shutil,zipfile
from datetime import datetime,timezone
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase81_allCAF_bone_lung_20261012'
now=datetime.now(timezone.utc).isoformat()
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def write(f,obj):
    tmp=f.with_suffix(f.suffix+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
assert (O/'logs/analysis.exit_code').read_text().strip()=='0'
qa=json.loads((O/'QA/NUMERIC_ACCEPTANCE.json').read_text());comp=json.loads((O/'COMPUTATION_COMPLETE.json').read_text())
assert qa['patient_alignment_verified'] and qa['max_P_error']<1e-9 and qa['max_BH_error']<1e-9
for info in json.loads((O/'resources/INPUT_MANIFEST.json').read_text()).values():assert sha(Path(info['path']))==info['sha256']
qa.update(status='PASS',visual_review_pending=False,visual_review='All three PNGs inspected, including revised figure1 FDR hit-count label; no clipping or overlap; displayed counts and directions agree with result tables.',accepted_at=now)
write(O/'QA/FINAL_ACCEPTANCE.json',qa)
report=O/'delivery/全部CAF候选_骨肺病灶比较_20261012.txt'
assert report.is_file() and (O/'tables/cross_source_nominal_signal_summary.tsv').is_file()
members=[f for sub in ['tables','figures','resources','code','QA','logs'] for f in (O/sub).rglob('*') if f.is_file()]+[report,O/'COMPUTATION_COMPLETE.json']
manifest={str(f.relative_to(O)):dict(sha256=sha(f),bytes=f.stat().st_size) for f in members}
write(O/'delivery/FILE_MANIFEST.json',manifest)
archive=O/'delivery/CAF_allCAF_bone_lung_20261012.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for f in members:z.write(f,arcname=str(f.relative_to(O)))
    z.write(O/'delivery/FILE_MANIFEST.json',arcname='delivery/FILE_MANIFEST.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for rel,info in manifest.items():assert hashlib.sha256(z.read(rel)).hexdigest()==info['sha256']
record=dict(status='completed_verified',phase=O.name,completed_at=now,title='全部CAF候选细胞的骨病灶与肺病灶比较',
    user_request='先用所有CAF细胞中的两转移灶做比对',interpretation='Approved scheme1 anatomical Bone versus Lung; all F classes pooled within each original candidate layer, not two selected individual samples',
    grouping='Anatomical site regardless of metastatic history',scope='All CAF candidates pooled within original supported1224/expanded8174; source-specific patient-level gene/program/proportion comparisons',
    new_inferential_analysis=True,main_testable_contexts=comp['main_testable_contexts'],nominal_gene_records=comp['main_nominal_gene_records'],gene_FDR_significant=comp['main_FDR_gene_records'],nominal_program_records=comp['main_nominal_program_records'],program_FDR_significant=comp['main_FDR_program_records'],nominal_proportion_records=comp['main_nominal_proportions'],
    FDR_family_policy='Valid P only, separately by source/layer/design/sensitivity/type; prior NA-as-one correction retained as diagnostic; defined before computation',
    report=str(report),archive=str(archive),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),acceptance=str(O/'QA/FINAL_ACCEPTANCE.json'),
    figures=[str(f) for f in sorted((O/'figures').glob('*.png'))],background_process_running=False,flow_marker_stage='deferred',remaining_required_work=[],
    local_delivery_verified=False,next_action='整体CAF骨肺比较已完成；后续扩展按用户新指示。')
write(O/'COMPLETE.json',record)
backup=O/'resources/state_before_allCAF';backup.mkdir(exist_ok=True)
state_names=['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json','CAF_SCHEME1_CURRENT.json']
for name in state_names:
    f=P/name
    if f.exists():
        assert not (backup/name).exists()
        shutil.copy2(f,backup/name)
for name in state_names[:2]:
    f=P/name;d=json.loads(f.read_text())
    for key in ['current_phase','current_stage','phase']:d[key]=O.name
    for key in ['status','current_status','current_state','task_status']:d[key]='allCAF_bone_lung_completed_verified'
    for key in ['updated_at','updated_at_utc','progress_updated_at','last_activity_at']:d[key]=now
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:d[key]=record
    d['active_task']=dict(title=record['title'],status='completed_verified',scope=record['scope'],result=str(O/'COMPLETE.json'))
    d.update(user_priority='全部CAF候选合并后的骨肺病灶部位比较',current_allCAF_comparison=str(P/'CAF_ALLCAF_BONE_LUNG_CURRENT.json'),current_cohort_differences=str(P/'CAF_ALLCAF_BONE_LUNG_CURRENT.json'),analysis_completion_scope=record['scope'],analysis_now_completed=True,background_process_running=False,acceptance_passed=True,analysis_results_report=str(report),next_action=record['next_action'],flow_marker_stage='deferred',previous_progress_snapshot=str(backup/name),local_delivery_verified=False,process_observation='New allCAF computation exited0, valid-P BH predefined, patient alignment and exact P verified, all three figures reviewed.')
    write(f,d)
resume=json.loads((P/'CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json').read_text());resume.update(status='allCAF_bone_lung_completed_verified',updated_at_utc=now,baseline_phase=O.name,latest_results=str(report),scope=record['scope'],user_request=record['user_request'],next_action=record['next_action'],background_process_running=False,current_allCAF_comparison=str(P/'CAF_ALLCAF_BONE_LUNG_CURRENT.json'))
write(P/'CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json',resume)
write(P/'CAF_ALLCAF_BONE_LUNG_CURRENT.json',record);write(P/'CAF_SCHEME1_CURRENT.json',record)
print(json.dumps(record,ensure_ascii=False))
