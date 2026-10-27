from pathlib import Path
import json,hashlib,zipfile,shutil,datetime
import pandas as pd
from pypdf import PdfReader
O=Path('C:/codex_test/CAF_external_delivery_20260921');D=O/'delivery';bundle=O/'client_package'
f=pd.read_csv(O/'tables/external_primary_associations.tsv',sep='\t');rows=[]
for r in f.to_dict('records'):rows.append({'signature':r['signature'],'cohort':r['cohort'],'validation_type':'external_bulk_outcome','endpoint':r['endpoint'],'performed':True,'status':'performed_no_clear_association_support','P':r['P'],'FDR':r['FDR'],'effect':r['effect'],'CI_low':r['CI_low'],'CI_high':r['CI_high']})
for v in ['external_single_cell_identity','surface_protein_sorting']:rows.append({'signature':'C7_and_C9','cohort':'not_run_this_stage','validation_type':v,'performed':False,'status':'not_performed; do not interpret as failed'})
pd.DataFrame(rows).to_csv(O/'tables/external_validation_evidence_status.tsv',sep='\t',index=False)
for p in [O/'tables/external_validation_evidence_status.tsv',O/'tables/customer_acceptance.tsv']:shutil.copy2(p,bundle/'external_validation'/p.name)
shutil.copy2(O/'code/validate.py',bundle/'methods/validate.py')
pdf=D/'骨肉瘤CAF阶段交付报告.pdf';reader=PdfReader(pdf);assert len(reader.pages)==4 and all(len(p.extract_text())>150 for p in reader.pages)
assert (f.FDR>.05).all() and len(f)==4
status={'status':'EXTERNAL_VALIDATION_COMPLETE_STAGE_DELIVERY_READY','external_bulk_validation':'performed; current data do not support clear C7/C9 outcome association','external_single_cell_identity':'not_performed','protein_sorting_validation':'not_performed','final_client_goal_achieved':False,'client_delivery_type':'stage_results_for_feedback','client_message_sent':False,'PDF_pages':4,'PDF_visual_QA':'all four pages inspected','clinical_cohorts':{'GSE21257':53,'GSE39055_RFS':36},'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');(D/'DELIVERY_QA.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(D/'DELIVERY_QA.json',bundle/'methods/DELIVERY_QA.json')
manifest={str(p.relative_to(bundle)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in bundle.rglob('*') if p.is_file() and p.name!='SHA256.json'}
(bundle/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=D/'CAF阶段性交付包_20260921.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(bundle.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(bundle))
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for k,v in manifest.items():assert hashlib.sha256(z.read(k)).hexdigest()==v
final={'PDF':{'path':pdf.name,'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()},'package':{'path':zpath.name,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'files':len(manifest)}}
(D/'FINAL_SHA256.json').write_text(json.dumps(final,ensure_ascii=False,indent=2),encoding='utf-8');(O/'DELIVERY_COMPLETE').write_text('STAGE_DELIVERY_READY_WITH_NEGATIVE_EXTERNAL_RESULTS\n')
print(json.dumps(status),flush=True);print('Package integrity PASS',len(manifest),'files',flush=True)
