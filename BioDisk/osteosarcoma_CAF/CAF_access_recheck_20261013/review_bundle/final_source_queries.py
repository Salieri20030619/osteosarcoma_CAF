from probe_sources import ROOT,fetch
import json,concurrent.futures,zipfile,xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
for name in ['hra_OMIX012450','hra_OMIX013436']:
 p=ROOT/'evidence'/f'{name}_direct.txt'
 if p.exists():
  s=BeautifulSoup(p.read_text(encoding='utf8'),'html.parser');t=s.get_text(' ',strip=True)
  (ROOT/'evidence'/f'{name}_extracted.txt').write_text(t,encoding='utf8');print(name,t[:6500])
p=ROOT/'evidence/ncc_real_methods_proxy.txt'
if p.exists() and zipfile.is_zipfile(p):
 with zipfile.ZipFile(p) as z:r=ET.fromstring(z.read('word/document.xml'))
 t='\n'.join(''.join(e.itertext()) for e in r.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}p'))
 (ROOT/'evidence/ncc_methods_extracted.txt').write_text(t,encoding='utf8');print('NCC METHODS',t[:11000])
jobs=[('hra_record','https://ngdc.cncb.ac.cn/search/api/specific?db=hra&q=HRA013407&sort=desc','direct'),('hra_osteo_records','https://ngdc.cncb.ac.cn/search/api/specific?db=hra&q=osteosarcoma&sort=desc&size=100','direct'),('ncc_repository_branches','https://api.github.com/repos/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/branches','proxy'),('ncc_repository_issues','https://api.github.com/repos/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/issues?state=all','proxy'),('shen_title_geo','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%22Aging%20reshapes%20the%20osteosarcoma%22&retmode=json','direct'),('shen_names_geo','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%22P2_pre_T%22%20OR%20%22P15_pre_T%22&retmode=json','direct')]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'final_source_query_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
