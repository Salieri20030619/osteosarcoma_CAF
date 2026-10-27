from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,zipfile
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012';now=datetime.now(timezone.utc).isoformat()
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def write(f,obj):
    tmp=f.with_suffix(f.suffix+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(f)
assert (O/'logs/analysis.exit_code').read_text().strip()=='0'
qa=json.loads((O/'QA/NUMERIC_ACCEPTANCE.json').read_text());lit=json.loads((O/'QA/LITERATURE_REVIEW.json').read_text());result=json.loads((O/'ANALYSIS_COMPLETE.json').read_text())
assert qa['gate_mismatch_cells']==0 and qa['source_effect_P_FDR_unchanged'] and lit['genes_covered']==14
base_manifest=json.loads((P/'phase82_allCAF_crossdataset_markers_20261012/resources/INPUT_MANIFEST.json').read_text())
assert sha(Path(base_manifest['stats']['path']))==base_manifest['stats']['sha256']
table=pd.read_csv(O/'tables/ALL14_LITERATURE_DIRECTION_INTEGRATION.tsv',sep='\t');assert len(table)==14 and table.gene.nunique()==14 and table.reference_urls.notna().all()
qa.update(status='PASS',visual_review_pending=False,visual_review='Both final PNGs inspected; source counts, sign agreement, zero/NA separation and recovery sensitivity are readable and agree with tables.',literature_genes=14,literature_reference_records=lit['reference_records'],literature_terms=lit['search_terms'],accepted_at=now)
write(O/'QA/FINAL_ACCEPTANCE.json',qa)
reports=sorted((O/'delivery').glob('*.txt'));assert len(reports)==2
archive=O/'delivery/CAF_upstream_direction_literature_20261012.zip'
members=[f for sub in ['code','tables','resources','QA','figures','logs'] for f in (O/sub).rglob('*') if f.is_file()]+reports+[O/'ANALYSIS_COMPLETE.json']
manifest={str(f.relative_to(O)):dict(sha256=sha(f),bytes=f.stat().st_size) for f in members};write(O/'delivery/FILE_MANIFEST.json',manifest)
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for f in members:z.write(f,arcname=str(f.relative_to(O)))
    z.write(O/'delivery/FILE_MANIFEST.json',arcname='delivery/FILE_MANIFEST.json')
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for rel,info in manifest.items():assert hashlib.sha256(z.read(rel)).hexdigest()==info['sha256']
record=dict(status='completed_verified_with_documented_upstream_limits',phase=O.name,completed_at=now,title='CAF上游证据核查 跨来源方向与既往文献联合评估',
    user_requests=['检查上游问题；解决不了则比较方向','查候选基因既往文献作为参考之一'],scope='30 source-sample retention audit; official site metadata;384-cell provisional late-CNV recovery sensitivity;14 markers across9-source registry with5 directional sources;14-gene literature review',
    primary_marker_count=14,new_inferential_P_values=0,formal_candidate_pool_unchanged=True,site_labels_changed=False,provisional_recovered_cells=384,independent_patients_added=0,recovery_direction_changes=0,
    unresolved_upstream=['Arm/Leg/Thigh tissue compartment not verified; original supplemental clinical table not obtained','Designated normal-reference CNV calls raise calibration/annotation concerns; no automatic normal/malignant reclassification','384 recovered entry-rule candidates have not undergone final cluster/lineage acceptance'],
    direction_no_opposite_candidates=['NEBL','GRID2','RGS2-AS1','AL096854.1','PEF1-AS1'],stronger_direction_candidate='NEBL:4same,1zero; not proven CAF-specific',
    clear_cross_source_discordance=['RTKN','TAF10','TOR1AIP1'],literature_gene_count=14,literature_reference_records=lit['reference_records'],literature_search_type='targeted reference review, not exhaustive systematic review',
    combined_research_leads=['NEBL:direction and adhesion/migration literature; isoform and CAF-identity limits','AL096854.1:OS univariate candidate association, not final signature member','GRID2:direction plus other-cancer association','RGS2-AS1:direction, limited direct mechanism literature','EBF1:bone/stroma mechanism reference with explicit discordance'],
    reports=[str(f) for f in reports],integrated_marker_table=str(O/'tables/ALL14_LITERATURE_DIRECTION_INTEGRATION.tsv'),direction_table=str(O/'tables/MARKER_DIRECTION_SUMMARY.tsv'),archive=str(archive),archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),acceptance=str(O/'QA/FINAL_ACCEPTANCE.json'),background_process_running=False,local_delivery_verified=False,remaining_required_work=[],next_action='本轮上游核查、方向复核和文献参考已完成；原14候选保留，后续选择依据联合证据。')
write(O/'COMPLETE.json',record)
backup=O/'resources/state_before_upstream_review';backup.mkdir(exist_ok=True)
names=['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json','CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json']
for name in names:
    f=P/name;assert not (backup/name).exists();shutil.copy2(f,backup/name)
for name in names[:2]:
    f=P/name;d=json.loads(f.read_text());assert d['current_phase']=='phase82_allCAF_crossdataset_markers_20261012'
    for key in ['current_phase','current_stage','phase']:d[key]=O.name
    for key in ['status','current_status','current_state','task_status']:d[key]='upstream_direction_literature_review_completed'
    for key in ['updated_at','updated_at_utc','progress_updated_at','last_activity_at']:d[key]=now
    for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:d[key]=record
    d.update(active_task=dict(title=record['title'],status=record['status'],scope=record['scope'],result=str(O/'COMPLETE.json')),user_priority='上游问题与跨来源方向复核；纳入既往文献参考',current_upstream_direction_review=str(P/'CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json'),analysis_results_report=record['reports'],analysis_completion_scope=record['scope'],analysis_now_completed=True,background_process_running=False,acceptance_passed=True,next_action=record['next_action'],local_delivery_verified=False,previous_progress_snapshot=str(backup/name),process_observation='30sample entry rules reproduced;384cell late evidence recovered in separate sensitivity;14 directions unchanged by recovery;14 literature records integrated;2figures verified.')
    write(f,d)
f=P/names[2];d=json.loads(f.read_text());d.update(updated_at_utc=now,status=record['status'],baseline_phase=O.name,latest_results=record['reports'],next_action=record['next_action'],current_upstream_direction_review=str(P/'CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json'));write(f,d)
f=P/names[3];d=json.loads(f.read_text());d.update(upstream_direction_literature_review=str(P/'CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json'),direction_no_opposite_candidates=record['direction_no_opposite_candidates'],cross_source_discordance_warning=record['clear_cross_source_discordance'],updated_at=now);write(f,d)
write(P/'CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json',record)
print(json.dumps(record,ensure_ascii=False))
