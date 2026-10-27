from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase59_four_functional_axes_20261004'
def rd(p):return json.loads(p.read_text(encoding='utf-8-sig'))
local=rd(O/'delivery/LOCAL_DELIVERY.json');arc=rd(O/'delivery/ZIP_ACCEPTANCE.json');acc=rd(O/'delivery/FINAL_ACCEPTANCE.json');old=rd(P/'CAF_NEXT_STAGE_CURRENT.json')
assert local['matches_server_archive'] and local['SHA256']==arc['SHA256'] and local['bytes']==arc['bytes']
assert acc['status']=='complete' and (O/'QA_COMPLETE').exists()
report=Path(acc['report']);assert hashlib.sha256(report.read_bytes()).hexdigest()==local['report_SHA256']
assert old.get('current_phase')=='phase58_F1F5_round1_20261004','Current phase changed; do not overwrite concurrent progress'
new={'project':'osteosarcoma_CAF','current_phase':O.name,'status':'four_axes_F2F4_source_attribution_complete','updated_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Original8174 F1-F5 labels preserved; four overlapping functional axes; F2 contraction/matrix and F4 osteogenic expression sources','report':str(report),'results_zip':arc['archive'],'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'local_delivery':str(O/'delivery/LOCAL_DELIVERY.json'),'local_delivery_verified':True,'previous_current_reference':str(O/'resources/PREVIOUS_CURRENT.json'),'original_F1F5_round1_phase':'phase58_F1F5_round1_20261004','formal_client_delivery_baseline':'phase49 unchanged','new_clustering':False,'CNV_rerun':False,'spatial_eggshell_geometry':'not_started','treatment_resistance_inference':'not_performed','scientific_decision':acc['key_findings'],'next_action':'Use lineage-separated F2/F4 results; no forced four-class labels and no whole-F4 osteogenic-CAF identity claim'}
tmp=P/'CAF_NEXT_STAGE_CURRENT.phase59.tmp.json';tmp.write_text(json.dumps(new,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(P/'CAF_NEXT_STAGE_CURRENT.json');(O/'COMPLETE').write_text('SUCCESS\n');print(json.dumps(new,ensure_ascii=False),flush=True)
