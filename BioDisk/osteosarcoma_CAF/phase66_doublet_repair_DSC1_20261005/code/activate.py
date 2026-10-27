from pathlib import Path
import json,hashlib,datetime,os
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005'
q=json.loads((O/'QA/LOCAL_DELIVERY_QA.json').read_text());assert q['passed']
assert json.loads((O/'QA/ACCEPTANCE.json').read_text())['passed']
m=json.loads((O/'resources/CURRENT_CELL_EVIDENCE_MANIFEST.json').read_text())
assert hashlib.sha256(Path(m['cell_table']).read_bytes()).hexdigest()==m['cell_table_sha256']
assert hashlib.sha256(Path(m['overlay']).read_bytes()).hexdigest()==m['overlay_sha256']
p=B/'CAF_NEXT_STAGE_CURRENT.json';d=json.loads(p.read_text());assert d['current_phase']=='phase65_spatial_structure_20261005'
newmanifest=B/'CAF_CURRENT_CELL_EVIDENCE.json'
if newmanifest.exists():(O/'resources/PREVIOUS_CELL_EVIDENCE_MANIFEST.json').write_bytes(newmanifest.read_bytes())
temp=B/'CAF_CURRENT_CELL_EVIDENCE.tmp.json';temp.write_text(json.dumps(m,indent=2));os.replace(temp,newmanifest)
d.update({'current_phase':O.name,'status':'doublet_metadata_repaired_sensitivity_verified_DSC1_checked_delivered','updated_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'report':str(O/'delivery/CAF双联体证据修复与SR1004核查报告_20261005.txt'),'results_zip':str(O/'delivery/CAF_doublet_repair_DSC1_20261005.zip'),'current_cell_metadata':m['cell_table'],'current_doublet_evidence':m['overlay'],'metadata_loader':m['loader'],'current_evidence_manifest':str(newmanifest),'marker_name_registry':m['marker_name_registry'],'acceptance':str(O/'QA/ACCEPTANCE.json'),'local_delivery_verified':True,'previous_current_reference':str(O/'resources/PREVIOUS_CURRENT.json'),'doublet_repair':{'assessed':8174,'restored_records':5728,'model_positive':26,'additional_positive_restored':15,'master_cells':8174,'sensitivity_retained':8148,'primary_source_effects_changed':0,'secondary_program_rows_lost_inference':132,'program_nominal_membership_changes':0,'gene_nominal_membership_changes':0,'historical_labels_confidence_preserved':True},'SR1004':{'entity':'antibody clone','product':'NBP3-22640','gene':'DSC1','RNA_assayed':8170,'raw_positive_cells':3,'raw_counts_per_positive_cell':1,'states':['F4_C8','F4_C9','F4_C12'],'interpretation':'Very low RNA detection; not supported as stable CAF marker; no protein absence claim.'},'next_action':'Use current_metadata.py/current evidence manifest in all future8174analyses; secondary-source inference uses updated doublet-exclusion N. Keep F2_C3 fibrillar candidate and vascular/osteogenic comparisons; independent high-resolution pathology and collagen deposition validation remain required.'})
tmp=B/'CAF_NEXT_STAGE_CURRENT.phase66.tmp.json';tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2));os.replace(tmp,p)
resolved={'at':d['updated_at'],'resolved':['doublet source metadata propagation','26model-positive-cell retain/exclude sensitivity'],'phase':O.name,'remaining_from_audit':['ambient fraction metadata provenance reconciliation beyond doublet scope','independent pathology and physical collagen measurements']}
(B/'review_intake_20261005_checklist/FOLLOWUP_RESOLUTION_PHASE66.json').write_text(json.dumps(resolved,indent=2))
(O/'PHASE66_DELIVERED.json').write_text(json.dumps({'at':d['updated_at'],'local_delivery_verified':True,'current_metadata_activated':True,'current_project_updated':True,'root_cause_fixed_in_successor_interface':True},indent=2))
print('PHASE66_ACTIVE',m['cell_table'],flush=True)
