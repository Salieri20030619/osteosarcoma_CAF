from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,zipfile
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase82_allCAF_crossdataset_markers_20261012';now=datetime.now(timezone.utc).isoformat()
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def write(f,d):
    tmp=f.with_suffix(f.suffix+'.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
assert (O/'logs/screen.exit_code').read_text().strip()=='0'
screen=json.loads((O/'QA/SCREEN_ACCEPTANCE.json').read_text());qa=json.loads((O/'QA/NUMERIC_ACCEPTANCE.json').read_text())
assert screen['primary_membership_reproduced'] and screen['max_raw_patient_expression_error']<1e-9 and qa['max_P_error']<1e-9
for info in json.loads((O/'resources/INPUT_MANIFEST.json').read_text()).values():assert sha(Path(info['path']))==info['sha256']
overlap=pd.read_csv(O/'tables/dataset_identity_overlap_check.tsv',sep='\t').fillna('')
assert overlap.shared_sample_ids.eq('').all() and overlap.shared_recorded_patient_ids.eq('').all()
qa.update(status='PASS',visual_review_pending=False,visual_review='All3 rendered PNGs inspected: all-marker effect plot, evidence agreement matrix, selected patient-expression panels. Labels and legends readable; numeric directions match tables.',raw_patient_expression_checks=screen['raw_patient_expression_checks'],max_raw_patient_expression_error=screen['max_raw_patient_expression_error'],accepted_at=now)
write(O/'QA/FINAL_ACCEPTANCE.json',qa)
report=O/'delivery/整体CAF_跨数据集marker筛选_20261012.txt'
members=[f for sub in ['code','tables','resources','QA','figures','logs'] for f in (O/sub).rglob('*') if f.is_file()]+[report,O/'SCREEN_COMPLETE.json']
manifest={str(f.relative_to(O)):dict(bytes=f.stat().st_size,sha256=sha(f)) for f in members};write(O/'delivery/FILE_MANIFEST.json',manifest)
archive=O/'delivery/CAF_allCAF_crossdataset_markers_20261012.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for f in members:z.write(f,arcname=str(f.relative_to(O)))
    z.write(O/'delivery/FILE_MANIFEST.json',arcname='delivery/FILE_MANIFEST.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for rel,info in manifest.items():assert hashlib.sha256(z.read(rel)).hexdigest()==info['sha256']
record=dict(status='completed_verified',phase=O.name,completed_at=now,title='整体CAF骨肺差异的跨数据集marker筛选',user_request='在这个层面重新寻找跨数据集支持的显著marker',
    scope='All CAF candidates, anatomical Bone vs Lung; same-layer same-direction raw P<0.05 replicated across source datasets',
    baseline='phase81_allCAF_bone_lung_20261012',primary_markers=14,lung_higher=13,bone_higher=1,FDR_supported_primary_markers=0,
    extra_supported17_nominal=qa['extra_supported17_nominal'],both_Fibro_direction_retained=qa['both_Fibro_direction_retained'],paired_both_same=qa['paired_both_same'],cross_layer_supplementary_markers=6,cross_layer_direction_conflicts=1,
    alias_correction='ENSG00000285280: source AL390957.1; frozen canonical name RGS2-AS1',
    all_Fibro_and_pair_directions_retained=[],all_LOPO_directions_retained=qa['all_LOPO_directions_retained'],all_doublet_nominal_retained=qa['all_doublet_nominal_retained'],
    marker_claim='Cross-dataset nominal expression candidates; not FDR-significant, not validated CAF-specific or flow markers',
    report=str(report),marker_table=str(O/'tables/MARKER_EVIDENCE_REGISTRY.tsv'),archive=str(archive),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),acceptance=str(O/'QA/FINAL_ACCEPTANCE.json'),
    new_inferential_tests=0,source_statistics_unchanged=True,background_process_running=False,flow_assay_stage='not started',local_delivery_verified=False,remaining_required_work=[],next_action='跨数据集marker筛选已完成；保留主14候选、跨层补充及矛盾证据。')
write(O/'COMPLETE.json',record)
backup=O/'resources/state_before_marker_screen';backup.mkdir(exist_ok=True)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json']:
    f=P/name;assert not (backup/name).exists();shutil.copy2(f,backup/name)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
    f=P/name;d=json.loads(f.read_text());assert d['current_phase']=='phase81_allCAF_bone_lung_20261012'
    for key in ['current_phase','current_stage','phase']:d[key]=O.name
    for key in ['status','current_status','current_state','task_status']:d[key]='allCAF_crossdataset_markers_completed_verified'
    for key in ['updated_at','updated_at_utc','progress_updated_at','last_activity_at']:d[key]=now
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:d[key]=record
    d.update(active_task=dict(title=record['title'],status='completed_verified',scope=record['scope'],result=str(O/'COMPLETE.json')),
      user_priority=record['user_request'],current_crossdataset_markers=str(P/'CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json'),analysis_results_report=str(report),analysis_completion_scope=record['scope'],analysis_now_completed=True,background_process_running=False,acceptance_passed=True,next_action=record['next_action'],local_delivery_verified=False,previous_progress_snapshot=str(backup/name),process_observation='Canonical feature-aligned cross-source screen; original P/FDR preserved;1362 patient expression records reconstructed;45 exact P checks;3 figures visually reviewed.')
    write(f,d)
f=P/'CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json';d=json.loads(f.read_text());d.update(status='allCAF_crossdataset_markers_completed_verified',updated_at_utc=now,baseline_phase=O.name,latest_results=str(report),scope=record['scope'],user_request=record['user_request'],next_action=record['next_action'],current_crossdataset_markers=str(P/'CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json'));write(f,d)
write(P/'CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json',record)
print(json.dumps(record,ensure_ascii=False))
