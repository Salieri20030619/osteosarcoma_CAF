from pathlib import Path
import urllib.request,json,xml.etree.ElementTree as ET
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012/resources/public_sources';O.mkdir(parents=True,exist_ok=True)
urls={'europepmc_xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12318355/fullTextXML','pmc_supplement':'https://pmc.ncbi.nlm.nih.gov/articles/instance/12318355/bin/NIHMS2025396-supplement-1.pdf'}
records=[]
opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
for key,url in urls.items():
    try:
        with opener.open(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=25) as r:data=r.read();ct=r.headers.get('Content-Type','');status=r.status
        ext='.pdf' if data.startswith(b'%PDF') else '.xml' if 'xml' in ct else '.html';f=O/(key+ext);f.write_bytes(data)
        records.append(dict(url=url,status=status,content_type=ct,bytes=len(data),file=str(f)))
        if ext=='.xml':
            tree=ET.fromstring(data)
            for e in tree.iter():
                if 'supplementary' in e.tag:print(e.attrib,ET.tostring(e,encoding='unicode')[:1000])
        if ext=='.pdf':
            try:
                from pypdf import PdfReader
                pages=[p.extract_text() or '' for p in PdfReader(f).pages];(O/(key+'.txt')).write_text('\n\n'.join(pages))
                for i,t in enumerate(pages):
                    if any(s in t for s in ['SARC0105','SARC0103','SARC0129']):print('PAGE',i+1,t[:13000])
            except ImportError:print('PDF saved; local extraction required')
    except Exception as e:records.append(dict(url=url,error=str(e)))
(O/'access_log.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
