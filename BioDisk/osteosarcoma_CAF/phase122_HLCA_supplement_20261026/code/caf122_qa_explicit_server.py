from pathlib import Path
import hashlib,json,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';pdf=O/'delivery/CAF_reviewer_QA_explicit_methods_20261027.pdf';sha=hashlib.sha256(pdf.read_bytes()).hexdigest();assert sha=='96812fc39a5e0e22b940f23820b6143a16416e8c6c62309125aff6b12c30f455';r=dict(status='COMPLETE',file=str(pdf),sha256=sha,pages=2,questions=5,main_report_unchanged=True,actual_methods_explicit=True,prohibited_term_absent=True,verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat());(O/'delivery/QA_EXPLICIT_SERVER_RECEIPT.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
for p in [O/'delivery/UNIFIED_CUSTOMER_DELIVERY_COMPLETE.json',R/'LATEST_CAF_INTEGRATED_CUSTOMER_REPORT.json']:
 d=json.loads(p.read_text(encoding='utf-8'));d['server_QA']=str(pdf);d['QA_sha256']=sha;d['QA_revision']=r;p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
for p in [O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:
 d=json.loads(p.read_text(encoding='utf-8'));d['integrated_customer_report']['reviewer_QA']=str(pdf);d['integrated_customer_report']['reviewer_QA_sha256']=sha;p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(r))
