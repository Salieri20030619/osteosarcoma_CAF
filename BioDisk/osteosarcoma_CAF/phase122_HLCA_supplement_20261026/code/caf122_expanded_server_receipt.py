from pathlib import Path
import hashlib,json,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';pdf=O/'delivery/CAF_reviewer_QA_expanded_reviewed_20261027.pdf';sha=hashlib.sha256(pdf.read_bytes()).hexdigest();assert sha=='611956c510e848d787db4628c36c4dbd0975934405c8d7476a31ffb5eb42e2d6';rec=dict(status='COMPLETE_SPECIFIED_SCOPE',file=str(pdf),sha256=sha,pages=4,questions_retained=11,addressed=['QA-01','QA-02','QA-03','QA-04','QA-05'],unchanged_questions=[1,2,3,5,9],Q7='Only banned generic term replaced; retention pending user decision',Q9='Denominator clarification outside this edit scope',verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat());(O/'delivery/EXPANDED_QA_REVIEWED_RECEIPT.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8')
for p in [R/'LATEST_CAF_INTEGRATED_CUSTOMER_REPORT.json',O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:
 d=json.loads(p.read_text(encoding='utf-8'));d['reviewer_QA_expanded']=rec;p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(rec))
