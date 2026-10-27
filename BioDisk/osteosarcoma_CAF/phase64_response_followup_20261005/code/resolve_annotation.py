from pathlib import Path
import requests,zlib,time,json,concurrent.futures,gzip,csv
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
def fetch(gpl):
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl[:-3]}nnn/{gpl}/soft/{gpl}_family.soft.gz'
 out=B/'sources'/f'{gpl}_platform_only.soft'
 if out.exists(): return [gpl,out.stat().st_size,'cached']
 start=time.time(); dec=zlib.decompressobj(16+zlib.MAX_WBITS); buf=b''; total=0
 with requests.get(url,stream=True,timeout=(20,35)) as r:
  r.raise_for_status()
  for c in r.iter_content(65536):
   total+=len(c); buf+=dec.decompress(c)
   if b'!platform_table_end' in buf:
    buf=buf[:buf.index(b'!platform_table_end')+len(b'!platform_table_end')]
    out.write_bytes(buf); return [gpl,len(buf),total,url]
   if total>30000000 or len(buf)>150000000 or time.time()-start>180: raise RuntimeError('bounded platform stream exceeded')
 raise RuntimeError('no platform end')
with concurrent.futures.ThreadPoolExecutor(2) as ex:
 for f in [ex.submit(fetch,g) for g in ['GPL14951','GPL10295']]:
  try: print(f.result(),flush=True)
  except Exception as e: print(type(e).__name__,str(e),flush=True)
for g in ['GSE39055','GSE21257']:
 with gzip.open(B/'sources'/f'{g}_series_matrix.txt.gz','rt') as f:
  for l in f:
   if l.startswith('!Sample_characteristics') or l.startswith('!Sample_title'):print(g,l[:1600].strip(),flush=True)
url='https://zenodo.org/api/records/19369637/files/Supplementary_Code.Rmd/content'
r=requests.get(url,timeout=60);r.raise_for_status();(B/'sources'/'GSE293065_Supplementary_Code.Rmd').write_bytes(r.content)
print('code',len(r.content))
