from pathlib import Path
import pymupdf,re,json,hashlib,shutil,os,datetime
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'expanded_QA_revision';src=W/'expanded_QA_review.pdf';backup=W/'original_expanded_QA.pdf';dst=Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答_扩充版.pdf');rec=json.loads((W/'EXPANDED_REVISION_QA.json').read_text(encoding='utf-8'))
def questions(path):
 d=pymupdf.open(path);t='\n'.join(p.get_text(sort=True) for p in d);matches=list(re.finditer(r'Q(\d+)\s*｜',t));return {int(m.group(1)):re.sub(r'\s+','',t[m.start():matches[i+1].start() if i+1<len(matches) else len(t)]) for i,m in enumerate(matches)}
old=questions(backup);new=questions(src);assert set(old)==set(new)==set(range(1,12))
for i in [1,2,3,5,9]:assert old[i]==new[i],('unexpected_change',i)
assert hashlib.sha256(src.read_bytes()).hexdigest()==rec['sha256'];tmp=dst.with_name(dst.name+'.tmp');shutil.copyfile(src,tmp)
try:os.replace(tmp,dst)
except PermissionError:
 with dst.open('r+b') as f:f.write(src.read_bytes());f.truncate();f.flush();os.fsync(f.fileno())
 tmp.unlink()
assert hashlib.sha256(dst.read_bytes()).hexdigest()==rec['sha256'];rec.update(status='COMPLETE',file=str(dst),unchanged_questions_verified=[1,2,3,5,9],visual_review='All4 pages inspected; no clipping or overlaps',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat());(W/'EXPANDED_DELIVERY_COMPLETE.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(dict(status='PASS',sha256=rec['sha256'],pages=4,questions=11,unchanged_questions=[1,2,3,5,9])))
