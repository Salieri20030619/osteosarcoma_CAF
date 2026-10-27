from probe_sources import ROOT,fetch
import json,concurrent.futures,zipfile,xml.etree.ElementTree as ET,re
from bs4 import BeautifulSoup
import openpyxl,io
p=ROOT/'evidence/shen_real_sup1_proxy.txt'
with zipfile.ZipFile(p) as z:
 root=ET.fromstring(z.read('word/document.xml'))
 text='\n'.join(''.join(e.itertext()) for e in root.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}p'))
 (ROOT/'evidence/shen_sup1_extracted.txt').write_text(text,encoding='utf8')
 print('SHEN DOCX extracted',len(text),'chars; accessions',re.findall(r'(?:GSE|GSM|HRA|CRA|PRJCA)\d+',text))
for n in [2,3]:
 p=ROOT/'evidence'/f'shen_real_sup{n}_proxy.txt'
 print('SUP',n,p.read_bytes()[:8])
 try:
  w=openpyxl.load_workbook(io.BytesIO(p.read_bytes()),read_only=True,data_only=True)
  result={s.title:list(s.values) for s in w}
 except Exception:
  soup=BeautifulSoup(p.read_text(encoding='utf8',errors='replace'),'html.parser');result={'text':soup.get_text(' ',strip=True)}
 (ROOT/'evidence'/f'shen_sup{n}_extracted.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,default=str),encoding='utf8')
 print(json.dumps(result,ensure_ascii=False,default=str)[:1000])
data=json.loads((ROOT/'evidence/shen_geo_esummary_direct.txt').read_text(encoding='utf8'))['result']
rows=[{'accession':data[i].get('accession'),'title':data[i].get('title'),'summary':data[i].get('summary'),'taxon':data[i].get('taxon'),'n_samples':data[i].get('n_samples')} for i in data['uids']]
(ROOT/'evidence/shen_geo_review.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
print('GEO TITLES',[(r['accession'],r['title']) for r in rows])
jobs=[('hra_specific','https://ngdc.cncb.ac.cn/search/specific?db=hra&q=HRA013407','direct'),('hra_OMIX012450','https://ngdc.cncb.ac.cn/omix/release/OMIX012450','direct'),('hra_OMIX013436','https://ngdc.cncb.ac.cn/omix/release/OMIX013436','direct'),('ncc_real_methods','https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13045-024-01598-7/MediaObjects/13045_2024_1598_MOESM4_ESM.docx','proxy'),('ncc_commits','https://api.github.com/repos/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/commits?path=BS6-matrix.h5ad&per_page=20','proxy')]
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:rows=list(ex.map(fetch,jobs))
(ROOT/'literature_next_routes.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
