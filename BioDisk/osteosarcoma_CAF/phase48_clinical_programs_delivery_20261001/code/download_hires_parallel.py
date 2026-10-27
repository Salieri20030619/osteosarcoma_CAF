from pathlib import Path
import os,sys,subprocess,json,re,gzip,hashlib,datetime,traceback,io
from concurrent.futures import ThreadPoolExecutor,as_completed
Q=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001');P=Q.parent;R=Q/'resources/histology_hires';R.mkdir(exist_ok=True)
if '--start' in sys.argv:
 f=(Q/'logs/hires_bounded_controller.log').open('w');p=subprocess.Popen([sys.executable,__file__,'--controller'],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'hires_parallel_controller.pid').write_text(str(p.pid));print('PARALLEL_DIRECT_DOWNLOAD_STARTED',p.pid)
elif '--controller' in sys.argv:
 p=subprocess.Popen([sys.executable,__file__,'--worker']);(Q/'hires_parallel_worker.pid').write_text(str(p.pid));rc=p.wait();(Q/'hires_download.exit_code').write_text(str(rc));(Q/'hires_parallel_native.exit_code').write_text(str(rc));print('NATIVE_HIRES_EXIT',rc);sys.exit(rc)
else:
 import requests
 import urllib3,socket
 urllib3.util.connection.allowed_gai_family=lambda:socket.AF_INET
 from PIL import Image
 from urllib.parse import urlparse
 import pandas as pd
 for k in list(os.environ):
  if k.lower().endswith('_proxy'):os.environ.pop(k,None)
 os.environ['NO_PROXY']='*';text=(P/'phase41_full_transcriptome_setup_20261001/resources/spatial_source_check/GEO_full_family.txt').read_text();meta=pd.read_csv(P/'phase36_external_validation_20260927/resources/sample_patient_manifest_corrected.tsv',sep='\t');urls={}
 for l in text.splitlines():
  if l.startswith('!Sample_supplementary_file') and l.endswith('tissue_hires_image.png.gz'):
   url=l.split(' = ',1)[1].replace('ftp://','https://');gsm=Path(urlparse(url).path).name.split('_')[0];urls[gsm]=url
 assert set(meta.gsm)<=set(urls)
 def task(row):
  url=urls[row.gsm];assert urlparse(url).hostname=='ftp.ncbi.nlm.nih.gov';dest=R/Path(urlparse(url).path).name;session=requests.Session();session.trust_env=False
  try:
   if not dest.exists():
    tmp=dest.with_suffix(dest.suffix+'.part');failures=0;total=None
    while total is None or (tmp.stat().st_size if tmp.exists() else 0)<total:
     offset=tmp.stat().st_size if tmp.exists() else 0;headers={'Range':f'bytes={offset}-{offset+524287}','Accept-Encoding':'identity'}
     try:
      with session.get(url,headers=headers,timeout=(30,90),stream=True) as response:
       response.raise_for_status();cr=response.headers.get('Content-Range','');m=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',cr)
       if response.status_code==206:
        assert m and int(m.group(1))==offset;total=int(m.group(3));expected=int(m.group(2))-offset+1
       else:
        assert response.status_code==200;offset=0;total=int(response.headers['Content-Length']);expected=total
       data=response.content;assert len(data)==expected
       with tmp.open('ab' if offset else 'wb') as f:f.write(data)
       failures=0
     except Exception:
      failures+=1
      if failures>=5:raise
    content=gzip.decompress(tmp.read_bytes());im=Image.open(io.BytesIO(content));im.load();tmp.replace(dest)
   content=gzip.decompress(dest.read_bytes());im=Image.open(io.BytesIO(content));im.load();return {'sample':row.sample,'gsm':row.gsm,'url':url,'path':str(dest),'bytes':dest.stat().st_size,'SHA256':hashlib.sha256(dest.read_bytes()).hexdigest(),'width':im.width,'height':im.height,'status':'complete','proxy':False}
  except Exception as e:traceback.print_exc();return {'sample':row.sample,'url':url,'status':'failed','error':repr(e),'proxy':False}
 rows=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(task,r) for r in meta.itertuples(index=False)]
  for future in as_completed(futures):
   r=future.result();rows.append(r);print('HIRES_OUTCOME',r['sample'],r['status'],r.get('bytes'),flush=True);(Q/'HIRES_DOWNLOAD_STATUS.json').write_text(json.dumps({'state':'running','done':len(rows),'total':len(meta),'outcomes':rows,'parallel':4,'proxy':False},indent=2))
 rc=int(any(r['status']!='complete' for r in rows));(Q/'HIRES_DOWNLOAD_STATUS.json').write_text(json.dumps({'state':'complete' if rc==0 else 'partial_failure','outcomes':rows,'parallel':4,'proxy':False},indent=2));(Q/'hires_download.exit_code').write_text(str(rc))
 if rc==0:(Q/'HIRES_DOWNLOAD_COMPLETE').write_text('SUCCESS')
 sys.exit(rc)
