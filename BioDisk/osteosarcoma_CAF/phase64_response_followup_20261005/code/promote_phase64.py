from pathlib import Path
import json,datetime,hashlib,os
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
local=json.loads((O/'QA/LOCAL_DELIVERY_QA.json').read_text());receipt=json.loads((O/'delivery/DELIVERY_RECEIPT.json').read_text())
assert local['passed'] and local['zip_sha256']==receipt['sha256']
assert json.loads((O/'FINAL_RUN_STATUS.json').read_text())['success']
current=json.loads((B/'CAF_NEXT_STAGE_CURRENT.json').read_text())
assert current['current_phase']=='phase63_matrix_resistance_20261005'
current.update({'current_phase':O.name,'status':'clinical_response_validation_and_spatial_source_audit_complete','updated_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'report':receipt['report'],'results_zip':receipt['zip'],'candidate_priorities':str(O/'tables/updated_candidate_decisions.tsv'),'local_delivery_verified':True,'acceptance':str(O/'QA/ACCEPTANCE.json'),'previous_current_reference':str(O/'resources/PREVIOUS_CURRENT.json'),'clinical_validation':{'GSE39055':{'submitted':37,'included':36,'poor':23,'good':13,'ambiguous':1},'GSE21257':{'submitted':53,'included':47,'poor':29,'good':18,'unknown':6},'programs':14,'result':'No program raw P<0.05 in either unadjusted or adjusted models; fixed fibrillar primary not supported as general poor-response association.','primary_differences':[0.197019,-0.039321],'primary_P':[0.399260,0.846810],'single_cohort_nominal_genes':[13,24],'both_cohorts_nominal_genes':0},'next_action':'Independent pathology contours and collagen protein/fiber localization; distinguish fibrillar collagen from typeIV vascular/necrotic-front and osteoid sources. Subsequent functional matrix/drug-response experiments conditional on actual deposition/localization evidence. Eggshell remains side branch.'})
current['conditional_endpoints']['direct_treatment_response']='Two independent public biopsy cohorts: 83 known responses; clinical association analysis complete, no universal collagen-low-response support.'
current['conditional_endpoints']['independent_tumor_boundary']='0 accepted masks after author code and source-data audit; IC Cancer_sweep>0.4 is expression-derived, not independent pathology.'
tmp=B/'CAF_NEXT_STAGE_CURRENT.phase64.tmp.json';tmp.write_text(json.dumps(current,ensure_ascii=False,indent=2));os.replace(tmp,B/'CAF_NEXT_STAGE_CURRENT.json')
receipt['local_verified']=True;receipt['local_QA']=str(O/'QA/LOCAL_DELIVERY_QA.json');(O/'delivery/DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
(O/'PHASE64_COMPLETE.json').write_text(json.dumps({'at':current['updated_at'],'analysis':True,'independent_QA':True,'local_delivery':True,'project_pointer_updated':True,'remaining':'Independent pathology contour, actual physical collagen readouts and controlled functional experiments; BS3 raw spatial download unavailable.'},indent=2))
print('PHASE64_DELIVERED_AND_CURRENT',current['current_phase'],flush=True)
