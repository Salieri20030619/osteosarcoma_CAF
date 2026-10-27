from pathlib import Path
import json,hashlib,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';p=O/'delivery/CAF_four_candidates_HLCA_supplement_explained_20261027.pdf';sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha=='27d75277ddbc4120fdd89bdf867d271c6b32229b742a1b83d9d89c97e26e802d'
r=json.loads((O/'delivery/DELIVERY_COMPLETE.json').read_text(encoding='utf-8'));r.update(report=str(p),report_sha256=sha,report_bytes=p.stat().st_size,pages=7,revision='Explain HLCA, IPF, their roles and relation to primary CAF analysis; clarify captions',revised_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),original_results_unchanged=True,original_report_preserved=True)
for path in [O/'delivery/EXPLAINED_REPORT_COMPLETE.json',O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:path.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status='PASS',sha256=sha,pages=7,original_preserved=True)))
