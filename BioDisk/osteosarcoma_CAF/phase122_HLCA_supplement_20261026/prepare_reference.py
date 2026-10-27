from pathlib import Path
import json,hashlib,datetime,requests,zipfile
W=Path(__file__).parent;D=W/'reference';D.mkdir(exist_ok=True)
items=[('HLCA_reference_model.zip','282060f9e4dc9281e02b1ba56b80cdaa'),('HLCA_reference_model_gene_order_ids_and_symbols.csv','53d367dbf4f30f766060ff92a83cd432')]
records=[]
for name,expected in items:
 url='https://zenodo.org/records/7599104/files/'+name+'?download=1';p=D/name;rec=dict(name=name,url=url,expected_md5=expected)
 try:
  if not p.exists():
   temp=p.with_suffix(p.suffix+'.part');r=requests.get(url,timeout=(10,30),stream=True);r.raise_for_status()
   with temp.open('wb') as f:
    for b in r.iter_content(1048576):f.write(b)
   assert hashlib.md5(temp.read_bytes()).hexdigest()==expected,'Official MD5 mismatch';temp.replace(p)
  data=p.read_bytes();assert hashlib.md5(data).hexdigest()==expected
  if name.endswith('.zip'):
   with zipfile.ZipFile(p) as z:assert z.testzip() is None;rec['archive_members']=z.namelist()
  rec.update(status='verified',bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
 except Exception as ex:rec.update(status='download_failed',error=str(ex))
 records.append(rec)
(W/'REFERENCE_PREPARATION.json').write_text(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),records=records,query_analysis_started=False),indent=2),encoding='utf-8');print(json.dumps(records))
