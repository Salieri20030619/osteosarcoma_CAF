from pathlib import Path
import requests,zlib,time,json,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
url='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14951/soft/GPL14951_family.soft.gz'
t=time.time();dec=zlib.decompressobj(31);buf=bytearray();total=0;last=0
with requests.get(url,stream=True,timeout=(20,40),headers={'Range':'bytes=0-20971519'}) as r:
 print('HTTP',r.status_code,r.headers.get('Content-Length'),r.headers.get('Content-Range'),flush=True);r.raise_for_status()
 for chunk in r.iter_content(65536):
  total+=len(chunk);d=dec.decompress(chunk);start=max(0,len(buf)-32);buf.extend(d)
  if total-last>1000000:print(total,len(buf),round(time.time()-t),flush=True);last=total
  end=buf.find(b'!platform_table_end',start)
  if end>=0:
   buf=buf[:end+len(b'!platform_table_end')]
   p=O/'sources/GPL14951_platform_only.soft'
   if not p.exists():p.write_bytes(buf)
   (O/'resources/GPL14951_RANGE_RECEIPT.json').write_text(json.dumps({'url':url,'range_request':'0-20971519','status':r.status_code,'compressed_read':total,'annotation_bytes':len(buf),'sha256':hashlib.sha256(buf).hexdigest()},indent=2))
   print('COMPLETE',len(buf),flush=True);break
  if total>20971520 or len(buf)>120000000 or time.time()-t>500:raise RuntimeError('range cap')
assert (O/'sources/GPL14951_platform_only.soft').exists()
