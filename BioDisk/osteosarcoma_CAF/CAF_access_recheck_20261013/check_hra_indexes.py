from probe_sources import ROOT,fetch
from bs4 import BeautifulSoup
import json,concurrent.futures
for name in ['hra_retry_search_direct','hra_retry_omix_direct']:
 soup=BeautifulSoup((ROOT/'evidence'/f'{name}.txt').read_text(encoding='utf8'),'html.parser')
 print('\nPAGE',name,soup.get_text(' ',strip=True)[-10000:])
 print('FORMS',[(str(f)[:2000]) for f in soup.find_all('form')])
 print('SCRIPTS',[s.get('src') for s in soup.find_all('script',src=True)])
jobs=[('hra_public_download_index','https://download.cncb.ac.cn/gsa-human/','direct'),('hra_search_all','https://ngdc.cncb.ac.cn/search/all?q=HRA013407','direct'),('hra_omix_browse','https://ngdc.cncb.ac.cn/omix/release','direct'),('shen_supp_pmc','https://pmc.ncbi.nlm.nih.gov/articles/PMC12816095/bin/13402_2026_1162_MOESM1_ESM.docx','proxy'),('ncc_supp_pmc','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11337771/supplementaryFiles','proxy')]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'index_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
