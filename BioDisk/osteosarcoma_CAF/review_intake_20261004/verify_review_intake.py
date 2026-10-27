from pathlib import Path
import csv,json,zipfile,hashlib,datetime
b=Path(r'C:\Users\31937\Downloads\CAF_blind_review_20261003\review_output')
zpath=Path(r'C:\Users\31937\Downloads\CAF盲审与揭盲比对_20261003.zip')
o=Path(r'C:\Users\31937\CAF_resume_20261004')
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
receipt=read(b/'DELIVERY_COMPLETE.json');assert sha(zpath)==receipt['SHA256']
with zipfile.ZipFile(zpath) as z:
 assert len(z.infolist())==receipt['files'] and z.testzip() is None
 for x in z.infolist():
  dest=(b/x.filename).resolve();assert dest.is_relative_to(b.resolve())
  assert dest.stat().st_size==x.file_size
  assert hashlib.sha256(z.read(x)).hexdigest()==sha(dest),x.filename
freeze=read(b/'blind_submission/BLIND_FREEZE.json')
assert sha(b/'blind_submission/BLIND_FREEZE_MANIFEST.tsv')==freeze['freeze_manifest_SHA256']
assert sha(b/'blind_submission/reviewer_response_template.tsv')==freeze['response_SHA256']
fr=rows(b/'blind_submission/BLIND_FREEZE_MANIFEST.tsv')
for x in fr:
 p=b/'blind_submission'/x['path'];assert p.stat().st_size==int(x['bytes']) and sha(p)==x['SHA256'],x['path']
comp=rows(b/'unblinded_comparison/cluster_annotation_comparison.tsv')
assert len(comp)==307 and sum(int(x['cells']) for x in comp)==876648
assert len({(x['source_id'],x['cluster_id']) for x in comp})==307
summary=rows(b/'unblinded_comparison/comparison_summary.tsv')
assert sum(int(x['clusters']) for x in summary)==307
assert sum(int(x['cells']) for x in summary)==876648
fib=rows(b/'unblinded_comparison/fibroblast_candidate_comparison.tsv')
assert len(fib)==15 and sum(int(x['cells']) for x in fib)==91012
result={'checked_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive':str(zpath),'archive_SHA256':sha(zpath),'archive_bytes':zpath.stat().st_size,'zip_CRC_passed':True,'archive_files_identical_to_local_review_output':len(z.infolist()),'frozen_files_verified':len(fr),'freeze':freeze,'clusters':len(comp),'cells':876648,'comparison_summary':summary,'fibroblast_candidate_clusters':len(fib),'fibroblast_candidate_cells':91012,'existing_fibroblast_clusters':sum(x['proposed_lineage']=='Fibroblast_like' for x in fib),'existing_fibroblast_cells':sum(int(x['cells']) for x in fib if x['proposed_lineage']=='Fibroblast_like'),'review_provenance_scope':'Freeze hashes and received metadata verified; reviewer independence and pre-freeze access history are statements in the supplied review, not independently observed in this recovery.','no_formal_labels_or_statistics_changed':True}
(o/'REVIEW_INTAKE_LOCAL_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ['comparison_summary','freeze']},ensure_ascii=True))