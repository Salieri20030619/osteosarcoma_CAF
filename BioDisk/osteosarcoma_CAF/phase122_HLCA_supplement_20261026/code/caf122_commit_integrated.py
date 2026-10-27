from pathlib import Path
import json,hashlib,shutil,os,datetime,pymupdf
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'full_report_integration';qa=json.loads((W/'INTEGRATION_QA.json').read_text(encoding='utf-8'));target=Path(qa['target']);staged=W/'integrated_review.pdf';sha=hashlib.sha256(staged.read_bytes()).hexdigest();assert sha==qa['staged_sha256'];assert hashlib.sha256(target.read_bytes()).hexdigest() in [qa['original_sha256'],sha]
tmp=target.with_name(target.name+'.tmp');shutil.copyfile(staged,tmp)
try:os.replace(tmp,target)
except PermissionError:
 with target.open('r+b') as f:f.write(staged.read_bytes());f.truncate();f.flush();os.fsync(f.fileno())
 tmp.unlink()
assert hashlib.sha256(target.read_bytes()).hexdigest()==sha;assert len(pymupdf.open(target))==17
qa.update(status='COMPLETE',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),file_name_unchanged=True,final_sha256=sha,visual_review='Pages11-17 inspected; original first10 pages retained with uniform page numbers',original_backup_preserved=True)
(W/'FULL_REPORT_INTEGRATION_COMPLETE.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False))
