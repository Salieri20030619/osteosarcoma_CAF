from probe_sources import ROOT,fetch
from bs4 import BeautifulSoup
import concurrent.futures,json,requests
jobs=[]
for name in ['hra_ngdc_root_direct','hra_human_download_root_direct']:
 soup=BeautifulSoup((ROOT/'evidence'/f'{name}.txt').read_text(encoding='utf8'),'html.parser')
 print(name,[(a.get_text(' ',strip=True)[:60],a.get('href')) for a in soup.find_all('a',href=True) if any(k in a.get('href','').lower() for k in ['human','gsa','omix'])])
for tag,path in [('home','gsa-human/'),('home_noslash','gsa-human'),('browse','gsa-human/browse'),('entry','gsa-human/browse/HRA013407'),('omix','omix/'),('search','search/?dbId=gsa_human&q=HRA013407')]:
 jobs.append(('hra_retry_'+tag,'https://ngdc.cncb.ac.cn/'+path,'direct'))
for n in (1,2,3):
 ext={1:'docx',2:'xls',3:'xlsx'}[n]
 jobs.append((f'shen_sup{n}',f'https://static-content.springer-cdn.com/esm/art%3A10.1007%2Fs13402-026-01162-0/MediaObjects/13402_2026_1162_MOESM{n}_ESM.{ext}','direct'))
jobs.append(('ncc_methods','https://static-content.springer-cdn.com/esm/art%3A10.1186%2Fs13045-024-01598-7/MediaObjects/13045_2024_1598_MOESM4_ESM.docx','direct'))
jobs.append(('shen_supp_epmc','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12816095/supplementaryFiles','proxy'))
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'followup_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
