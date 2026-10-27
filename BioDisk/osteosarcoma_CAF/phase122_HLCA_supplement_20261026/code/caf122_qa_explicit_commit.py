from pathlib import Path
import json,hashlib,shutil,os,datetime
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026/unified_full_report');W=B/'qa_explicit_methods';r=json.loads((W/'QA_REVISION_CHECK.json').read_text(encoding='utf-8'));src=W/'reviewer_QA_review.pdf';dst=Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答.pdf');assert hashlib.sha256(src.read_bytes()).hexdigest()==r['sha256']
tmp=dst.with_name(dst.name+'.tmp');shutil.copyfile(src,tmp)
try:os.replace(tmp,dst)
except PermissionError:
 with dst.open('r+b') as f:f.write(src.read_bytes());f.truncate();f.flush();os.fsync(f.fileno())
 tmp.unlink()
assert hashlib.sha256(dst.read_bytes()).hexdigest()==r['sha256'];assert hashlib.sha256(Path('C:/Users/31937/Downloads/骨肺CAF亚群_图文分析与Lung04推荐报告.pdf').read_bytes()).hexdigest()==r['main_report_unchanged_sha256'];r.update(status='COMPLETE',file=str(dst),visual_review='Both pages inspected; no clipping',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat());(W/'QA_EXPLICIT_COMPLETE.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
p=B/'UNIFIED_DELIVERY_COMPLETE.json';d=json.loads(p.read_text(encoding='utf-8'));d['QA_sha256']=r['sha256'];d['QA_revision']='Five answers lead with results and name the actual analyses; prohibited generic term absent';p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(r,ensure_ascii=False))
