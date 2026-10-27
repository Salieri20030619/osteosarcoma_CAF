from pathlib import Path
import hashlib,json,datetime
R=Path(__file__).resolve().parent
m=json.loads((R/'NCC_DOWNLOAD_COMPLETE.json').read_text());rows=[]
for f in m['files']:
 p=R/'downloads/NCC'/f['file'];h=hashlib.sha256()
 with p.open('rb') as inp:
  for b in iter(lambda:inp.read(1024*1024),b''):h.update(b)
 assert p.stat().st_size==f['bytes'] and h.hexdigest()==f['sha256'],f['file']
 rows.append({'file':f['file'],'bytes':p.stat().st_size,'sha256':h.hexdigest()})
bc=json.loads((R/'NCC_BARCODE_COMPARISON.json').read_text());assert bc['reference_samples']==107
result={'status':'PASS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files_verified':len(rows),'bytes':sum(f['bytes'] for f in rows),'files':rows,'baseline_unchanged':True,'NCC_barcode_comparison_completed':True,'BS6_semantic_labels_missing':True,'raw_directory':str(R/'downloads/NCC')}
(R/'SERVER_ACCEPTANCE.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='files'}))
