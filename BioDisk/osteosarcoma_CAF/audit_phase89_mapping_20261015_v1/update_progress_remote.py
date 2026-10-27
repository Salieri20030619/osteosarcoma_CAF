from pathlib import Path
import json, hashlib, datetime, shutil, subprocess

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=Path(__file__).resolve().parent
assert O.parent==P and O.name.startswith('audit_phase89_mapping_20261015_')
REF=P/'phase89_fixed_state_mapping_20261015'
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(f): return json.loads(Path(f).read_text(encoding='utf-8'))
def sha(f):
    h=hashlib.sha256()
    with Path(f).open('rb') as s:
        for b in iter(lambda:s.read(8388608),b''): h.update(b)
    return h.hexdigest()
def dump(f,o):
    f=Path(f); temp=f.with_name(f.name+'.audit_tmp')
    with temp.open('x',encoding='utf-8',newline='\n') as s:
        s.write(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    assert read(temp)==o
    temp.replace(f)

accept=read(O/'ACCEPTANCE.json')
assert accept['computational_consistency_pass'] and accept['status']=='PASS_EXPLORATORY_REVIEW_ONLY'
expected=read(O/'REVIEWED_PHASE89_HASHES.json')
assert all(sha(REF/n)==h for n,h in expected.items()), 'Reviewed Phase89 files changed'
inputs=read(REF/'resources/INPUT_HASHES.json')['before']
assert all(sha(f)==h for f,h in inputs.items()), 'Registered upstream input changed'
assert sha(P/'audit_expanded_admission_20261015/tables/final_cell_review_decisions.tsv.gz')==inputs[str(P/'audit_expanded_admission_20261015/tables/final_cell_review_decisions.tsv.gz')]
saved=O/'source_evidence'; saved.mkdir(exist_ok=False)
for f in [P/'search_patient_expansion_20261012_v1/sources/NCC_patient_metadata.txt',P/'phase88_expansion_identity_20261014/resources/CELLxGENE_AUTHOR_METADATA_CHECK.json',P/'phase88_expansion_identity_20261014/tables/patient_sample_roster.tsv']:
    shutil.copy2(f,saved/f.name)

ps=subprocess.run(['ps','-u','wangmingyuan','-o','pid,etimes,stat,args','--width','240'],capture_output=True,text=True,check=True)
running=[s for s in ps.stdout.splitlines() if any(t in s for t in ['osteosarcoma_CAF','copykat','infercnv']) and 'update_progress_remote.py' not in s and 'base64' not in s]
# A review cannot silently interfere with another active scientific worker.
if running: raise RuntimeError('Potential project process found: '+repr(running))

backup=O/'previous_state'; backup.mkdir(exist_ok=False)
names=['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_EXPANSION_IDENTITY_CURRENT.json']
original={}
for n in names:
    f=P/n; original[n]=sha(f); shutil.copy2(f,backup/n)
for n in ['CAF_FIXED_STATE_MAPPING_CURRENT.json','CAF_PHASE89_REVIEW_CURRENT.json']:
    assert not (P/n).exists(), 'Do not overwrite a newer mapping progress pointer'

next_action='Retain source-stratified exploratory mappings. Resolve CAF-EXP-01/03/04 before formal expansion; satisfy CAF-MAP-01 before stronger state-identity claims or inference using transferred labels. Await user-forwarded CC review/revision; do not auto-launch repair.'
review={'phase':REF.name,'audit_dir':str(O),'status':accept['status'],'completed_utc':now,'report':str(O/'REVIEW_ZH.txt'),'english_feedback':str(O/'CC_FEEDBACK_EN.txt'),'acceptance':str(O/'ACCEPTANCE.json'),'baseline_phase':'phase86_functional_bone_lung_20261012','baseline_cells':8174,'baseline_states':33,'query_cells':1268,'cell_counts':accept['cell_counts'],'new_formal_patients':0,'new_Bone_Lung_patients':0,'patient_metadata_review':accept['patient_metadata_review'],'issue_status':accept['issue_status'],'next_action':next_action,'analysis_artifacts_modified':False}
dump(P/'CAF_PHASE89_REVIEW_CURRENT.json',review)
mapping=read(REF/'MAPPING_COMPLETE.json')
mapping.update({'original_completion_status':mapping['status'],'status':'MAPPING_REVIEWED_EXPLORATORY_ONLY','review':review,'patient_independence_review':'existing_source_metadata_reviewed_cross_study_and_anatomical_evidence_pending','new_formal_patients':0,'new_Bone_Lung_patients':0,'updated_utc':now})
dump(P/'CAF_FIXED_STATE_MAPPING_CURRENT.json',mapping)
for name in names[:2]:
    state=read(P/name)
    assert sha(P/name)==original[name]
    for k in ['phase','current_phase']: state[k]=REF.name
    for k in ['status','task_status','current_status','current_state']: state[k]='phase89_review_complete_exploratory_only_formal_expansion_pending'
    for k in ['checked_at','updated_at','updated_at_utc','progress_updated_at','last_activity_at','progress_reconciled_at']: state[k]=now
    state.update({'current_stage':'fixed_state_mapping_review_completed','baseline_phase':'phase86_functional_bone_lung_20261012','formal_baseline_cells':8174,'formal_baseline_states':33,'latest_task':review,'last_completed_task':review,'active_task':{'phase':REF.name,'status':'completed','scope':'read-only mapping audit and existing patient-metadata review'},'next_action':next_action,'next_analysis_plan':next_action,'current_mapping_review':str(P/'CAF_PHASE89_REVIEW_CURRENT.json'),'current_fixed_state_mapping':str(P/'CAF_FIXED_STATE_MAPPING_CURRENT.json'),'background_process_running':False,'new_analysis_started':False,'new_inferential_tests':False,'analysis_now_completed':True,'analysis_completion_scope':'Phase89 descriptive mapping review and archived patient metadata only; no formal expansion','acceptance_passed':True,'acceptance_scope':'exploratory review only; no validated identity or clinical admission','acceptance':str(O/'ACCEPTANCE.json'),'report':str(O/'REVIEW_ZH.txt'),'latest_report':str(O/'REVIEW_ZH.txt'),'local_report':'C:\\Users\\31937\\CAF_phase89_review_20261015\\REVIEW_ZH.txt','results_zip':None,'local_delivery_verified':False,'progress_reconciliation_record':str(O/'PROGRESS_UPDATE_RECEIPT.json'),'previous_progress_snapshot':str(backup/name),'process_observation':'No relevant scientific process observed during publication','execution_status':{'scope':'independent audit_mapping.py','exit_code':0,'original_phase89_exit_code_available':False},'followup_recommendations':[next_action],'user_priority':'Continue Phase89 review and reconcile project progress','legacy_result_sections':'Preserved historical scientific results; use current_mapping_review for current work and baseline_phase for formal baseline'})
    dump(P/name,state)
identity=read(P/names[2]); assert sha(P/names[2])==original[names[2]]
identity.update({'status':'computation_complete_expanded_manual_review_complete_mapping_review_complete','manual_review_pending':False,'expanded_manual_review':str(P/'audit_expanded_admission_20261015/ACCEPTANCE.json'),'subsequent_mapping_review':str(P/'CAF_PHASE89_REVIEW_CURRENT.json'),'baseline_admission_automatic':False,'new_formal_CAF_patients':0,'updated_utc':now})
dump(P/names[2],identity)

# Verify the new progress and original scientific artifacts after the metadata-only update.
assert all(sha(REF/n)==h for n,h in expected.items())
assert all(sha(f)==h for f,h in inputs.items())
for n,h in original.items(): assert sha(backup/n)==h
for n in names[:2]:
    s=read(P/n)
    assert s['current_phase']==REF.name and s['baseline_phase']=='phase86_functional_bone_lung_20261012'
    assert s['formal_baseline_cells']==8174 and not s['background_process_running']
outputs={str(f.relative_to(O)):sha(f) for f in O.rglob('*') if f.is_file() and 'input_phase89' not in f.parts}
receipt={'status':'PROGRESS_UPDATED_VERIFIED','utc':now,'audit_dir':str(O),'review_status':accept['status'],'updated_progress_files':{n:sha(P/n) for n in names+['CAF_FIXED_STATE_MAPPING_CURRENT.json','CAF_PHASE89_REVIEW_CURRENT.json']},'state_backup_hashes':original,'reviewed_Phase89_unchanged':True,'four_registered_upstream_hashes_unchanged':True,'scientific_baseline_unchanged':True,'running_scientific_processes':running,'deliverable_hashes':outputs,'new_formal_patients':0,'new_Bone_Lung_patients':0}
dump(O/'PROGRESS_UPDATE_RECEIPT.json',receipt)
print(json.dumps({'status':receipt['status'],'audit_dir':str(O),'review_status':accept['status'],'progress_files':list(receipt['updated_progress_files']),'scientific_baseline_unchanged':True},ensure_ascii=False))
