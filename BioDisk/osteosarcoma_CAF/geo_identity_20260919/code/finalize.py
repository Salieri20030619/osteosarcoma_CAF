from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'geo_identity_20260919'
assert (O/'logs/identity.exit_code').read_text().strip()=='0'
assert (O/'CAF_GEO_identity_review_20260919.md').exists()
assert (O/'VISUAL_QA_PASS.json').exists()
state=json.loads((O/'STATUS.json').read_text())
state.update(status='IDENTITY_REVIEW_COMPLETE_CANDIDATES_CLASSIFIED',figure_visual_QA='PASS',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),next_action='Use versioned supported candidate sets for any separately scoped downstream analysis; retain unresolved cells, FLEX validation and cross-cohort patient mapping boundaries')
(O/'STATUS.json').write_text(json.dumps(state,indent=2))
(O/'COMPLETE').write_text('CAF_IDENTITY_REVIEW_COMPLETE_CANDIDATES_CLASSIFIED\n')
# Record source evidence without modifying previous analysis stages.
srcs=[P/'geo_cnv_20260918/result_review_20260919/REVIEW_STATUS.json',P/'geo_cnv_20260918/result_review_20260919/tables/cell_review.tsv.gz',P/'geo_preprocess_20260918/manifest.tsv',P/'geo_preprocess_20260918/review/lineage_rules.json',P/'phase7_review_states_20260917/code/integrate_identity.py']
rows=[]
for f in srcs:
 with open(f,'rb') as h: rows.append({'path':str(f),'sha256':hashlib.file_digest(h,'sha256').hexdigest()})
pd.DataFrame(rows).to_csv(O/'tables/source_evidence_sha256.tsv',sep='\t',index=False)
paths=[p for name in ['code','tables','figures','objects','logs'] for p in (O/name).rglob('*') if p.is_file() and '__pycache__' not in str(p)]
paths += [O/name for name in ['CAF_GEO_identity_review_20260919.md','STATUS.json','decision_rules.json','gene_panels.json','VISUAL_QA_PASS.json','COMPLETE','README.md','python_environment.json']]
rows=[]
for f in sorted(paths):
 with open(f,'rb') as h: rows.append({'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.file_digest(h,'sha256').hexdigest()})
pd.DataFrame(rows).to_csv(O/'SHA256SUMS.tsv',sep='\t',index=False)
# Preserve a timestamped backup, then sync the main project resume pointer.
main=P/'PROJECT_STATE.json';timestamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
shutil.copy2(main,P/('PROJECT_STATE.before_identity_'+timestamp+'.json'))
project=json.loads(main.read_text())
project.update(updated_at_utc=state['completed_utc'],current_stage='geo_identity_20260919',dynamic_status_file='geo_identity_20260919/STATUS.json',status=state['status'],current_status=state['status'],next_action=state['next_action'],latest_identity_report=str(O/'CAF_GEO_identity_review_20260919.md'))
project['geo_identity_20260919']=state
project['geo_cnv_result_review_20260919']={'status':'COMPLETE_IDENTITY_ACCEPTANCE_OPEN_AT_REVIEW','status_file':'geo_cnv_20260918/result_review_20260919/REVIEW_STATUS.json','subsequent_identity_review':'geo_identity_20260919/STATUS.json'}
project['completed_stages']=list(dict.fromkeys(project.get('completed_stages',[])+['geo_cnv_20260918','geo_identity_20260919']))
tmp=P/'PROJECT_STATE.identity_update.tmp'
tmp.write_text(json.dumps(project,indent=2,ensure_ascii=False),encoding='utf-8');tmp.replace(main)
print(json.dumps(state),flush=True)
