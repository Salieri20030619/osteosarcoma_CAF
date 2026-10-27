from probe_sources import ROOT,fetch
from bs4 import BeautifulSoup
import json,concurrent.futures,re,pdfplumber
with pdfplumber.open(ROOT/'evidence/hra_supplement_proxy.txt') as pdf:
 texts=['PAGE '+str(i+1)+'\n'+(p.extract_text() or '') for i,p in enumerate(pdf.pages)]
(ROOT/'evidence/hra_supplement_extracted.txt').write_text('\n\n'.join(texts),encoding='utf8')
print('HRA SUPPLEMENT',len(texts),'pages','\n'.join(texts)[:6000])
soup=BeautifulSoup((ROOT/'evidence/hra_omix_list_direct.txt').read_text(encoding='utf8'),'html.parser')
found=[]
for tr in soup.find_all('tr'):
 txt=tr.get_text(' ',strip=True)
 if re.search('osteosarcoma|HRA013407|cdr.2026.09|Shen|Paired single',txt,re.I):
  found.append({'text':txt,'links':[a['href'] for a in tr.find_all('a',href=True)]})
(ROOT/'evidence/OMIX_CANDIDATE_ROWS.json').write_text(json.dumps(found,ensure_ascii=False,indent=2),encoding='utf8');print('OMIX',json.dumps(found,ensure_ascii=False))
soup=BeautifulSoup((ROOT/'evidence/shen_springer_page_proxy.txt').read_text(encoding='utf8'),'html.parser')
links=list(dict.fromkeys(a['href'] for a in soup.find_all('a',href=True) if 'MOESM' in a['href']))
print('REAL SHEN SUP LINKS',links)
jobs=[('shen_real_sup'+str(i+1),u,'proxy') for i,u in enumerate(links)]
jobs += [('hra_ngdc_api','https://ngdc.cncb.ac.cn/search/api/ngdc?q=HRA013407','direct'),('hra_osteo_api','https://ngdc.cncb.ac.cn/search/api/ngdc?q=osteosarcoma','direct')]
ids=json.loads((ROOT/'evidence/shen_geo_esearch_direct.txt').read_text())['esearchresult']['idlist']
jobs.append(('shen_geo_esummary','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&id='+','.join(ids),'direct'))
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'verified_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
