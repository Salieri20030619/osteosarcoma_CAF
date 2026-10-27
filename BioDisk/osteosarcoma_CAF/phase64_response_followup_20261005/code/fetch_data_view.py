from pathlib import Path
import requests,time,json,concurrent.futures,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
def fetch(g):
 p=O/'sources'/f'{g}_platform_only.soft';url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=self&form=text&view=data'
 if p.exists():return {'gpl':g,'cached':True}
 t=time.time();buf=bytearray()
 with requests.get(url,stream=True,timeout=(15,45)) as r:
  r.raise_for_status()
  for c in r.iter_content(262144):
   buf.extend(c)
   if len(buf)%2097152<262144:print(g,len(buf),'bytes',round(time.time()-t),'sec',flush=True)
   end=buf.find(b'!platform_table_end')
   if end>=0:buf=buf[:end+len(b'!platform_table_end')];break
   if len(buf)>80000000 or time.time()-t>900:raise RuntimeError('bounded fetch exceeded')
 assert b'!platform_table_end' in buf
 p.write_bytes(buf);print('PLATFORM_READY',g,len(buf),flush=True)
 return {'gpl':g,'url':url,'bytes':len(buf),'sha256':hashlib.sha256(buf).hexdigest()}
with concurrent.futures.ThreadPoolExecutor(2) as ex:r=list(ex.map(fetch,['GPL14951','GPL10295']))
(O/'resources'/'PLATFORM_SUCCESS_RECEIPTS.json').write_text(json.dumps(r,indent=2))
for g in ['GPL14951','GPL10295']:
 s=(O/'sources'/f'{g}_platform_only.soft').read_text();i=s.index('!platform_table_begin');print(g,s[i:i+1500],flush=True)
