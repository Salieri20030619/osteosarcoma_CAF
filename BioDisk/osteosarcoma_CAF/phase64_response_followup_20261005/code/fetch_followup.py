from pathlib import Path
import requests,json,hashlib,gzip,concurrent.futures,xml.etree.ElementTree as ET,re,io,pandas as pd,shutil
from bs4 import BeautifulSoup
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase64_response_followup_20261005'
def fetch(url,name):
 path=O/'sources'/name
 try:
  if path.exists():data=path.read_bytes()
  else:
   rr=requests.get(url,timeout=(15,75));rr.raise_for_status();data=rr.content
   if name.endswith('.gz'):gzip.decompress(data)
   path.write_bytes(data)
  return {'url':url,'path':str(path),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'status':'available'}
 except Exception as e:return {'url':url,'path':str(path),'status':'failed','error':str(e)}
targets=[


('https://zenodo.org/api/records/14594702','Zenodo14594702.json'),
('https://zenodo.org/api/records/19369637','Zenodo19369637.json'),
('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3706900/fullTextXML','GSE39055_paper.xml'),
('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7730477/fullTextXML','GSE152048_paper.xml'),
('https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-74603-x/MediaObjects/41467_2026_74603_MOESM6_ESM.xlsx','GSE293065_source_data.xlsx')]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 records=[]
 for record in ex.map(lambda v:fetch(*v),targets):records.append(record);print('FETCH',record['status'],Path(record['path']).name,record.get('bytes',record.get('error')),flush=True)
(O/'resources/FOLLOWUP_FETCH_RECEIPTS.json').write_text(json.dumps(records,indent=2))
for gpl in ['GPL14951','GPL10295']:
 path=O/'sources'/f'{gpl}_family.soft.gz'
 if path.exists():
  txt=gzip.open(path,'rt').read();block=txt.split('!platform_table_begin')[1].split('!platform_table_end')[0];f=pd.read_csv(io.StringIO(block.strip()),sep='\t',low_memory=False);print('PLATFORM',gpl,f.shape,list(f.columns));print(f.head(2).to_string(index=False)[:2400])
for num in [14594702,19369637]:
 p=O/'sources'/f'Zenodo{num}.json'
 if p.exists():
  j=json.loads(p.read_text());print('ZENODO',num,j.get('metadata',{}).get('title'));print('ZENODO_FILES',[(f['key'],f.get('size'),f.get('links',{})) for f in j.get('files',[])])
for name in ['Zheng2024_fulltext.xml','GSE39055_paper.xml','GSE152048_paper.xml']:
 path=O/'sources'/name
 if not path.exists():continue
 root=ET.fromstring(path.read_text());plain=' '.join(root.itertext());(O/'sources'/name.replace('.xml','_plain.txt')).write_text(plain)
 hits=[plain[max(0,m.start()-140):m.end()+240] for m in re.finditer(r'huvos|necrosis|necrotic|90%|response|chemoresis|annotation|boundary|BS3|data availability',plain,re.I)]
 (O/'sources'/name.replace('.xml','_excerpts.txt')).write_text('\n\n'.join(hits));print('PAPER',name,'HITS',len(hits));print('\n'.join(hits[:8]))
 links=[(n.tag,n.attrib) for n in root.iter() if n.tag in ['supplementary-material','media','ext-link']];(O/'resources'/name.replace('.xml','_links.json')).write_text(json.dumps(links,indent=2))
 print('SUPPLEMENT_LINKS',[x for x in links if x[0]=='media'][:12])
p=O/'sources/GSE293065_source_data.xlsx'
if p.exists():
 sheets=pd.ExcelFile(p).sheet_names;print('SOURCE_DATA_SHEETS',sheets);ledger=[]
 for sh in sheets:
  f=pd.read_excel(p,sheet_name=sh,header=None);txt=' '.join(f.fillna('').astype(str).to_numpy().ravel());hits=re.findall(r'.{0,80}(?:huvos|necrosis|necrotic|qupath|boundary|mask|annotation|x_pixel|y_pixel).{0,100}',txt,re.I);ledger.append({'sheet':sh,'rows':len(f),'columns':len(f.columns),'targeted_hits':hits[:20]})
 (O/'resources/SPATIAL_SOURCE_DATA_AUDIT.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2));print('SOURCE_DATA_HITS',[r for r in ledger if r['targeted_hits']])
# Record supplementary table1 as audited source data, independently of integrated metadata.
st=P/'phase35_external_validation_prep_20260927/sources/GSE293065_supplement.txt';shutil.copy2(st,O/'sources/GSE293065_supplement_original.txt')
(O/'SOURCE_AUDIT_READY').write_text('SUCCESS\n')
