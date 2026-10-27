from pathlib import Path
import json,h5py,numpy as np,xml.etree.ElementTree as ET,zipfile,re
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parent
for p in sorted((ROOT/'downloads/NCC').glob('*')):
 if p.suffix not in ['.h5','.h5ad']:continue
 with h5py.File(p) as f:
  print('H5',p.name,'ROOT',list(f.keys()))
  if 'matrix' in f:print('MATRIX',list(f['matrix'].keys()),f['matrix/shape'][:].tolist())
  else:
   for key in f:
    print(key,list(f[key].keys())[:30] if isinstance(f[key],h5py.Group) else f[key].shape,dict(f[key].attrs))
print('METADATA',(ROOT/'downloads/NCC/Patient metadata').read_text())
for name in ['hra_search_all_direct','hra_retry_omix_direct']:
 soup=BeautifulSoup((ROOT/'evidence'/f'{name}.txt').read_text(encoding='utf8'),'html.parser')
 print('PAGE',name)
 print('RELEVANT LINKS',[(a.get_text(' ',strip=True),a.get('href')) for a in soup.find_all('a',href=True) if any(w in a.get('href','') for w in ['browse','search','release'])])
 for s in soup.find_all('script'):
  if s.get('src'): print('SRC',s['src'])
  elif any(w in s.get_text() for w in ['axios','ajax','url','HRA013407']): print(s.get_text()[:9000])
try:
 from pypdf import PdfReader
 p=ROOT/'evidence/hra_supplement_proxy.txt';r=PdfReader(p)
 texts=['PAGE '+str(i+1)+'\n'+page.extract_text() for i,page in enumerate(r.pages)]
 (ROOT/'evidence/hra_supplement_extracted.txt').write_text('\n\n'.join(texts),encoding='utf8')
 print('HRA SUPPLEMENT',len(r.pages),'pages','\n'.join(texts)[:11000])
except ImportError as e:print('PDF extraction import',str(e))
