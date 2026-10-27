from pathlib import Path
import requests,json,hashlib,gzip,csv,io,re,datetime,concurrent.futures
from bs4 import BeautifulSoup
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase64_response_followup_20261005'
for d in ['resources','sources','tables','objects','figures','QA','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def fetch(url,name):
 path=O/'sources'/name
 try:
  r=requests.get(url,timeout=(15,60));r.raise_for_status();data=r.content
  if name.endswith('.gz'):gzip.decompress(data)
  path.write_bytes(data);return {'url':url,'path':str(path),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'status':'downloaded'}
 except Exception as e:return {'url':url,'path':str(path),'status':'failed','error':str(e)}
targets=[
('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39055/matrix/GSE39055_series_matrix.txt.gz','GSE39055_series_matrix.txt.gz'),
('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21257/matrix/GSE21257_series_matrix.txt.gz','GSE21257_series_matrix.txt.gz'),
('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11337771/fullTextXML','Zheng2024_fulltext.xml'),
('https://api.zenodo.org/records/14594702','Zenodo14594702.json')]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:receipts=list(ex.map(lambda v:fetch(*v),targets))
(O/'resources/FETCH_RECEIPTS.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
for r in receipts:
 print('FETCH',r['status'],r['url'],r.get('bytes',r.get('error')),flush=True)
 if r['status']=='downloaded' and r['path'].endswith('matrix.txt.gz'):
  text=gzip.open(r['path'],'rt').read();head=text.split('!series_matrix_table_begin')[0];(O/'sources'/(Path(r['path']).stem+'.metadata.txt')).write_text(head,encoding='utf-8')
  for line in head.splitlines():
   if line.startswith(('!Series_title','!Series_pubmed_id','!Sample_platform_id','!Sample_data_processing')):print(line[:1600])
   if line.startswith('!Sample_characteristics_ch1'):
    vals=next(csv.reader([line],delimiter='\t'))[1:];print('CHARACTERISTIC',sorted(set(vals))[:35])
# Targeted local original-source audit, not generated annotation files.
source_dirs=[P/'phase35_external_validation_prep_20260927/sources',P/'phase37_requested_source_search_20260927/sources']
inventory=[]
for root in source_dirs:
 for path in root.rglob('*'):
  if path.is_file() and path.stat().st_size<30000000 and path.suffix.lower() in ['.pdf','.xlsx','.xls','.txt','.json','.html','.xml','.md'] and not path.name.endswith('.part'):
   inventory.append({'path':str(path),'bytes':path.stat().st_size})
(O/'resources/LOCAL_ORIGINAL_SOURCE_INVENTORY.json').write_text(json.dumps(inventory,indent=2))
article=P/'phase35_external_validation_prep_20260927/sources/GSE293065_article.html'
soup=BeautifulSoup(article.read_text(),'html.parser');plain=soup.get_text(' ',strip=True);(O/'sources/GSE293065_article_plain.txt').write_text(plain,encoding='utf-8')
for heading in ['Data availability','Code availability']:
 h=soup.find(lambda t:t.name in ['h2','h3'] and t.get_text(strip=True)==heading)
 if h:
  sec=h.find_parent('section');print('ARTICLE_SECTION',heading,sec.get_text(' ',strip=True)[:3500])
links=[{'text':a.get_text(' ',strip=True),'href':a.get('href')} for a in soup.find_all('a',href=True) if any(k in a['href'].lower() for k in ['moesm','github','zenodo'])]
(O/'resources/SPATIAL_SOURCE_LINKS.json').write_text(json.dumps(links,ensure_ascii=False,indent=2))
print('SPATIAL_LINKS',json.dumps(links[-25:],ensure_ascii=False))
repo=P/'phase35_external_validation_prep_20260927/sources/Zheng2024_readme.json'
obj=json.loads(repo.read_text());print('ZHENG_REPO_JSON_TYPE',type(obj).__name__)
if isinstance(obj,dict):print({k:obj.get(k) for k in ['name','html_url','download_url','message']})
# Scan supplementary text for endpoint or independent annotation references.
for path in [P/'phase35_external_validation_prep_20260927/sources/GSE293065_supplement.txt',O/'sources/Zheng2024_fulltext.xml']:
 if not path.exists():continue
 txt=BeautifulSoup(path.read_text(),'xml').get_text(' ',strip=True) if path.suffix=='.xml' else path.read_text()
 hits=[]
 for m in re.finditer(r'huvos|necrosis|necrotic|chemoresis|patholog|annotation|boundary|qupath|mask',txt,re.I):hits.append(txt[max(0,m.start()-100):m.end()+180])
 (O/'sources'/(path.stem+'_targeted_excerpts.txt')).write_text('\n\n'.join(hits),encoding='utf-8');print('EXCERPTS',path.name,len(hits));print('\n'.join(hits[:6]))
print('SOURCE_INSPECTION_COMPLETE',flush=True)
