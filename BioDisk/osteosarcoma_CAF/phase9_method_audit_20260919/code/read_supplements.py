from pathlib import Path
import urllib.request,zipfile,xml.etree.ElementTree as ET,json
p=Path(r'C:\Users\31937\caf_method_audit_20260919\resources');p.mkdir(exist_ok=True)
rows=[]
for n in [5,6]:
 url=f'https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12943-025-02516-2/MediaObjects/12943_2025_2516_MOESM{n}_ESM.docx'
 try:
  with urllib.request.urlopen(url,timeout=25) as r:data=r.read()
  f=p/f'MOESM{n}.docx';f.write_bytes(data)
  with zipfile.ZipFile(f) as z:root=ET.fromstring(z.read('word/document.xml'))
  ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};text='\n'.join(''.join(t.itertext()) for t in root.findall('.//w:t',ns));(p/f'MOESM{n}.txt').write_text(text,encoding='utf-8');print('SUPPLEMENT',n,len(data),text[:400]);rows.append({'url':url,'success':True})
 except Exception as e:print('SUPPLEMENT',n,str(e));rows.append({'url':url,'success':False,'error':str(e)})
(p/'supplement_fetch.json').write_text(json.dumps(rows,indent=2))