from pathlib import Path
import json,hashlib,shutil,datetime
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase14_client_objectives_20260920'
assert int((O/'exit_code').read_text())==0 and (O/'ANALYSIS_COMPLETE').exists()
assert json.loads((O/'REVIEW_QA.json').read_text())['all_checks_passed']
m=json.loads((O/'DELIVERY_SHA256.json').read_text())
for n,h in m.items():assert hashlib.sha256((O/n).read_bytes()).hexdigest()==h,n
c=pd.read_csv(O/'tables/01_clinical_sample_audit_updated.tsv',sep='\t');assert len(c)==40
assert c.loc[c.metastasis_any.eq('No'),'unified_strict'].sum()==136
s=json.loads((O/'STATUS.json').read_text());s['initial_A_strict_CAF_before_supplement']=s.pop('A_strict_CAF')
s.update(status='CLIENT_OBJECTIVE_REVIEW_COMPLETE_DEPENDENCIES_OPEN',report=str(O/'客户目标逐项推进.md'),clinical_current_table='tables/01_clinical_sample_audit_updated.tsv',A_candidate_CAF_after_crosspublication_update=136,metastasis_label_yes_CAF_GSE162454=260,ScPCA_metadata_pair_groups=4,ScPCA_expression_pairs_available=0,patient3_myCAF=196,patient6_myCAF=89,final_client_goal_achieved=False,steps={'clinical_audit':'complete_with_new_labels_and_pending_lung_timing','customer_panels':'complete; group-specific new subpopulation endpoint unresolved','surface_screen':'preparatory complete; no final metastasis-specific panel'},next_required=['ScPCA paired-expression files or usable download links','OS3 OS4 metastatic organ and timing; patient alias confirmation when possible','biological replication and surface protein validation before final panel'])
(O/'STATUS.json').write_text(json.dumps(s,ensure_ascii=False,indent=2))
backup=P/'PROJECT_STATE.before_phase14_20260920.json'
if not backup.exists():shutil.copy2(P/'PROJECT_STATE.json',backup)
p=json.loads((P/'PROJECT_STATE.json').read_text());p.update(current_stage=O.name,status=s['status'],current_status=s['status'],dynamic_status_file=str(O/'STATUS.json'),latest_report=s['report'],updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase14_client_objectives_20260920=s,active_objective='Customer requirement: metastasis-related CAF/myCAF candidate and flow-sort marker panel; no full-paper scope expansion',next_action='Prioritize clinical/pair-expression dependencies; use 2025 SEC16B S1 and original 2021 S1 cross-mapping for explicitly exploratory any-metastasis comparisons; FAP+GPRC5A is provisional and not lung-general')
(P/'PROJECT_STATE.json').write_text(json.dumps(p,ensure_ascii=False,indent=2));(O/'REVIEW_COMPLETE').write_text('CUSTOMER_OBJECTIVE_PROGRESS_REVIEWED_DEPENDENCIES_OPEN\n');print(json.dumps({'status':s['status'],'verified_files':len(m)}))
