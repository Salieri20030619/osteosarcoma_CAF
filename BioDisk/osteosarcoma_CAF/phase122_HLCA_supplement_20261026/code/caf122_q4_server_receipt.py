from pathlib import Path
import json,hashlib,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';p=O/'delivery/CAF_reviewer_QA_Q4_revision_20261027.pdf';sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='725596aef42ad1dde61390c6f99adc871d689a665515377657c1f90d9f35e6b3';r=dict(status='COMPLETE',file=str(p),sha256=sha,changed_question='Q4 only',unchanged_questions=['Q1','Q2','Q3','Q5'],pages=2,main_report_unchanged=True,verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat());(O/'delivery/Q4_REVISION_RECEIPT.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
for f in [O/'delivery/UNIFIED_CUSTOMER_DELIVERY_COMPLETE.json',R/'LATEST_CAF_INTEGRATED_CUSTOMER_REPORT.json']:
 d=json.loads(f.read_text(encoding='utf-8'));d['server_QA']=str(p);d['QA_sha256']=sha;d['QA_revision']=r;f.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
for f in [O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:
 d=json.loads(f.read_text(encoding='utf-8'));d['integrated_customer_report']['reviewer_QA']=str(p);d['integrated_customer_report']['reviewer_QA_sha256']=sha;f.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(r))
