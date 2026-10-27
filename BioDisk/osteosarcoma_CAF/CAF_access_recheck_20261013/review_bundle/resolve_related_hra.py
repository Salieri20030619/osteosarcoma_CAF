from probe_sources import ROOT,fetch
import json,concurrent.futures
from bs4 import BeautifulSoup
data=json.loads((ROOT/'evidence/hra_osteo_all22_direct.txt').read_text(encoding='utf8'))['result']['data']['data']
ids=['HRA016637','HRA007229','HRA008101','HRA007225','HRA008525','HRA006496']
index=(ROOT/'evidence/hra_public_download_index_direct.txt').read_text(encoding='utf8')
jobs=[]
for r in data:
 if r['id'] not in ids:continue
 project=r['attrs']['bioProject Accession']
 print(r['id'],project,'public_directory_listed',r['id']+'/' in index)
 jobs.append((r['id']+'_project','https://ngdc.cncb.ac.cn/search/api/specific?db=bioproject&q='+project,'direct'))
 samples=r['attrs'].get('bioSample Accession',[])
 if samples:jobs.append((r['id']+'_sample','https://ngdc.cncb.ac.cn/search/api/specific?db=biosample&q='+samples[0],'direct'))
 if r['id']+'/' in index:jobs.append((r['id']+'_public_directory','https://download.cncb.ac.cn/gsa-human/'+r['id']+'/','direct'))
jobs.append(('HRA008101_published_share','https://ngdc.cncb.ac.cn/gsa-human/s/70lNUz7v','direct'))
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'related_hra_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
for r in rows:
 if r.get('status')==200 and r['name'].endswith(('_project','_sample')):
  try:print('RECORD',r['name'],json.loads((ROOT/'evidence'/r['file']).read_text(encoding='utf8'))['result']['data']['data'])
  except Exception:pass
