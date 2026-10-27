from pathlib import Path
import requests,time,json,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
g='GPL14951';url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=self&form=text&view=data'
records=[]
for attempt in range(1,4):
 buf=bytearray();t=time.time();last=0
 try:
  with requests.get(url,stream=True,timeout=(20,60),headers={'Accept-Encoding':'gzip'}) as r:
   r.raise_for_status()
   for c in r.iter_content(65536):
    buf.extend(c)
    if len(buf)-last>5000000:print(attempt,len(buf),round(time.time()-t),flush=True);last=len(buf)
    end=buf.find(b'!platform_table_end')
    if end>=0:buf=buf[:end+len(b'!platform_table_end')];break
    if len(buf)>80000000 or time.time()-t>700:raise RuntimeError('bound exceeded')
  assert b'!platform_table_end' in buf
  p=O/'sources'/f'{g}_platform_only.soft';p.write_bytes(buf)
  records.append({'attempt':attempt,'success':True,'bytes':len(buf),'url':url,'sha256':hashlib.sha256(buf).hexdigest()})
  print('SUCCESS',len(buf),flush=True);break
 except Exception as e:
  (O/'sources'/f'{g}_incomplete_attempt{attempt}.txt').write_bytes(buf)
  records.append({'attempt':attempt,'success':False,'bytes':len(buf),'error':str(e)})
  print('RETRY',attempt,len(buf),str(e)[:150],flush=True)
(O/'resources/GPL14951_RETRY_RECEIPTS.json').write_text(json.dumps(records,indent=2))
assert (O/'sources'/f'{g}_platform_only.soft').exists()
