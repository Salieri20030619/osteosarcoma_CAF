from probe_sources import ROOT,fetch
import requests,json,concurrent.futures
from bs4 import BeautifulSoup
s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
r=s.post('https://api.figshare.com/v2/articles/search',json={'search_for':'Osteosarcoma Cells Secrete CXCL14','limit':100},timeout=45)
(ROOT/'evidence/CXCL14_FIGSHARE_SEARCH.json').write_bytes(r.content);r.raise_for_status()
items=r.json();print('FIGSHARE RESULTS',[(i['id'],i['title']) for i in items])
jobs=[]
for i in items:
 if 'Osteosarcoma Cells Secrete CXCL14' in i['title'] and any(k in i['title'] for k in ['Table S1A','Table S2C','Materials and Methods']):jobs.append(('CXCL14_'+str(i['id']),i['url'],'proxy'))
for p in ['PRJCA025405','PRJCA028332','PRJCA029010','PRJCA019918']:
 jobs.append((p+'_all_samples','https://ngdc.cncb.ac.cn/search/api/specific?db=biosample&q='+p+'&length=100','direct'))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'related_material_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
