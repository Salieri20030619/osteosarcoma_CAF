import requests,concurrent.futures,time,zlib,json
from pathlib import Path
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
def probe(item):
 name,url=item
 try:
  t=time.time()
  with requests.get(url,stream=True,timeout=(10,25)) as r:
   c=next(r.iter_content(16384));d=zlib.decompressobj(31).decompress(c) if c[:2]==b'\x1f\x8b' else c
   (O/'sources'/f'probe_{name}.txt').write_bytes(d)
   print(name,r.status_code,len(c),round(time.time()-t),d[:1000].decode(errors='replace'),flush=True)
 except Exception as e:print(name,str(e),flush=True)
with concurrent.futures.ThreadPoolExecutor(4) as ex:list(ex.map(probe,[
 ('gpl14951_soft','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14951/soft/GPL14951_family.soft.gz'),
 ('gpl14951_data','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951&targ=self&form=text&view=data'),
 ('gpl10295_data','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10295&targ=self&form=text&view=data'),
 ('gpl14951_full','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951&targ=self&form=text&view=full')]))
