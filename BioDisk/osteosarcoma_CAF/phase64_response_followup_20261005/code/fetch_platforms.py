from pathlib import Path
import requests,time,json,hashlib,gzip,concurrent.futures
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
def fetch(gpl):
 urls=[f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl[:-3]}nnn/{gpl}/annot/{gpl}.annot.gz',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full']
 records=[]
 for i,url in enumerate(urls):
  try:
   start=time.time()
   with requests.get(url,timeout=(10,25),stream=True) as r:
    r.raise_for_status();length=int(r.headers.get('Content-Length','0'));assert length<80000000,f'oversized response {length}'
    parts=[];size=0
    for chunk in r.iter_content(262144):
     parts.append(chunk);size+=len(chunk)
     if size>80000000 or time.time()-start>120:raise RuntimeError('bounded download size/time exceeded')
   data=b''.join(parts);text=gzip.decompress(data).decode() if i==0 else data.decode()
   assert ('!platform_table_begin' in text or '!platform_table_begin' in text.lower() or '#ID' in text) and len(text)>10000,'not a platform annotation'
   path=O/'sources'/f'{gpl}_annotation.txt';path.write_text(text,encoding='utf-8');records.append({'url':url,'path':str(path),'status':'available','bytes':len(data),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()});print('PLATFORM_READY',gpl,len(text),flush=True);break
  except Exception as e:records.append({'url':url,'status':'failed','error':str(e)});print('PLATFORM_ATTEMPT',gpl,str(e)[:180],flush=True)
 return records
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 result=list(pool.map(fetch,['GPL14951','GPL10295']))
(O/'resources/PLATFORM_FETCH_RECEIPTS.json').write_text(json.dumps(sum(result,[]),indent=2))
assert all((O/'sources'/f'{x}_annotation.txt').exists() for x in ['GPL14951','GPL10295'])
print('PLATFORMS_COMPLETE',flush=True)
