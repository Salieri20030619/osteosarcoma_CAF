from pathlib import Path
import json,hashlib,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';p=O/'delivery/CAF_full_customer_report_integrated_20261027.pdf';sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='e295aa0bc7fde9f4cd5227a9e057270102325e94b886eb334ca12b66c3a7d637'
qa=json.loads((O/'code/FULL_REPORT_INTEGRATION_COMPLETE.json').read_text(encoding='utf-8'));assert qa['final_sha256']==sha;qa.update(server_report=str(p),server_hash_verified=True,server_verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'delivery/FULL_REPORT_INTEGRATION_SERVER_RECEIPT.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');(R/'LATEST_CAF_INTEGRATED_CUSTOMER_REPORT.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
for pointer in [O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:
 d=json.loads(pointer.read_text(encoding='utf-8'));d['integrated_customer_report']=dict(path=str(p),sha256=sha,pages=17,local_original_name_preserved=True,main_HLCA_pages=[12,13,14],supplement_pages=[15,16,17]);pointer.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status='PASS',sha256=sha,pages=17)))
