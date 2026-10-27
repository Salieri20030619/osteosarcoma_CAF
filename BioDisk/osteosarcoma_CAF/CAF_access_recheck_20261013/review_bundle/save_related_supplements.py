from pathlib import Path
import requests,json,hashlib,concurrent.futures,openpyxl,io,pdfplumber
ROOT=Path(__file__).resolve().parent;O=ROOT/'literature/CXCL14';O.mkdir(parents=True,exist_ok=True)
rows=[]
for article in [25515870,25515867,25515846]:
 d=json.loads((ROOT/'evidence'/f'CXCL14_{article}_proxy.txt').read_text(encoding='utf8'))
 for f in d['files']:
  s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
  r=s.get(f['download_url'],timeout=60);r.raise_for_status();p=O/f['name'];p.write_bytes(r.content)
  assert len(r.content)==f['size'] and hashlib.md5(r.content).hexdigest()==f['computed_md5']
  row={'article':article,'title':d['title'],'file':f['name'],'bytes':f['size'],'md5':f['computed_md5'],'url':f['download_url']};rows.append(row)
  print('SAVED',row)
  if p.suffix=='.xlsx':
   w=openpyxl.load_workbook(p,read_only=True,data_only=True);content={sh.title:list(sh.values) for sh in w}
   (O/(p.stem+'_extracted.json')).write_text(json.dumps(content,ensure_ascii=False,default=str,indent=2),encoding='utf8')
   print('TABLE',p.name,json.dumps(content,ensure_ascii=False,default=str)[:5000])
  elif p.suffix=='.pdf':
   with pdfplumber.open(p) as pdf:txt='\n'.join(page.extract_text() or '' for page in pdf.pages)
   (O/(p.stem+'_extracted.txt')).write_text(txt,encoding='utf8');print('METHODS',txt[:4000])
(ROOT/'RELATED_SUPPLEMENTS_COMPLETE.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
for project in ['PRJCA025405','PRJCA028332','PRJCA029010','PRJCA019918']:
 d=json.loads((ROOT/'evidence'/f'{project}_all_samples_direct.txt').read_text(encoding='utf8'))['result']['data']
 print('SAMPLES',project,d['recordsTotal'],[(r['id'],r['attrs'].get('SampleName'),r['description']) for r in d['data']])
