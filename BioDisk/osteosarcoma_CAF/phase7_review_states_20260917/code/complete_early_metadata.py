"""SYNTAX_CHECK: python -m py_compile code/complete_early_metadata.py
RUN_COMMAND: python code/complete_early_metadata.py
EXPECTED_SUCCESS_MARKER: EARLY_PHASE_METADATA_AUDIT_RECORDED
"""
from pathlib import Path
import json,hashlib,datetime
import pandas as pd
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=r/'phase7_review_states_20260917'
clinical=pd.read_csv(r/'phase4_download_20260914/clinical_source_audit.tsv',sep='\t',encoding='utf-8-sig').fillna('')
# Independently transcribed from the already-present public Supplementary Table 1 image, PDF page 25.
extra={'BC2':('Conventional','5.5*5*3','<90%',50),'BC3':('Conventional','8*6*6','<90%',70),'BC5':('Conventional','8*7.5*6','>=90%',80),'BC6':('Conventional','7*7*4','>=90%',15),'BC10':('Conventional','3.5*3*2','<90%',8),'BC11':('Conventional','20*11*10','<90%',30),'BC16':('Conventional','6*4*2.5','<90%',40),'BC17':('Chondroblastic','2*2*1.5','<90%',40),'BC20':('Chondroblastic','10*8*5','<90%',50),'BC21':('Intraosseous osteosarcoma','5*4*1','<90%',10),'BC22':('Chondroblastic','18*15*12','<90%',20)}
for col in ['pathological_type_raw','tumor_size_cm_raw','necrosis_rate_raw','Ki67_percent_raw','age_sex_source_status']:clinical[col]='not_extracted_from_other_source'
for i,row in clinical.iterrows():
 if row.dataset!='GSE152048':continue
 hist,size,necrosis,ki=extra[row['sample']]
 for col,value in [('pathological_type_raw',hist),('tumor_size_cm_raw',size),('necrosis_rate_raw',necrosis),('Ki67_percent_raw',str(ki))]:clinical.loc[i,col]=value
 clinical.loc[i,'age_sex_source_status']='absent_in_reviewed_Supplementary_Table_1; patient-level GEO mapping remains to verify'
 clinical.loc[i,'sex_raw']='not_reported_in_reviewed_table';clinical.loc[i,'age_raw']='not_reported_in_reviewed_table'
clinical.to_csv(p/'review/clinical_enriched_20260917.tsv',sep='\t',index=False)
manifest=json.loads((r/'phase4_download_20260914/download_manifest.json').read_text());rows=[]
for item in manifest:
 for gsm in item['geo_sample_ids']:
  rows.append({'dataset':item['dataset'],'patient_id':item['patient_id'],'patient_mapping_status':'known_within_dataset','tissue_id':item['dataset']+':'+item['sample'],'analysis_sample_id':item['dataset']+'__'+item['sample'],'submitted_sample_id':gsm,'library_id':'UNRESOLVED','library_mapping_status':'GSM-to-tissue known; independent capture/library not recoverable from tissue matrix alone','expression_input':'phase4_download_20260914/matrices/'+item['dataset']+'/'+item['sample'],'source':'phase4_download_20260914/download_manifest.json'})
old=pd.read_csv(r/'OS_CAF_phase2/code/manifest.tsv',sep='\t')
for row in old.to_dict('records'):
 rows.append({'dataset':'GSE270231','patient_id':row['patient_id'],'patient_mapping_status':'unresolved' if row['patient_id']=='UNCONFIRMED' else 'known_within_dataset','tissue_id':'GSE270231:'+row['short'],'analysis_sample_id':row['sample'],'submitted_sample_id':row['sample'],'library_id':'UNRESOLVED','library_mapping_status':'sample-to-object known; separate capture/library ID not established','expression_input':'OS_CAF_phase2/results/'+row['sample']+'/annotated.h5ad','source':'OS_CAF_phase2/code/manifest.tsv'})
pd.DataFrame(rows).to_csv(p/'review/active_patient_tissue_sample_library_master.tsv',sep='\t',index=False)
notes={'updated_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'master_scope':'24 actually analysed tissue/sample units; one row per submitted GEO sample, not one biological replicate','library_mapping':'unknown fields remain explicit; GSM is not assumed to equal a capture library','clinical_image_source':'phase4_download_20260914/clinical_page_25.png','clinical_image_sha256':hashlib.sha256((r/'phase4_download_20260914/clinical_page_25.png').read_bytes()).hexdigest(),'clinical_source_table':'PMC7730477 Supplementary Table 1, PDF page 25','changes':'pathology, size, necrosis and Ki67 extracted; age/sex absence in reviewed table distinguished from unread data','ScPCA_source_table':'existing local audit retained on Windows; upload denied by automatic approval; not silently synchronized','no_new_clinical_groups_inferred':True}
(p/'review/early_metadata_audit.json').write_text(json.dumps(notes,indent=2));print('EARLY_PHASE_METADATA_AUDIT_RECORDED',len(rows),len(clinical),flush=True)
