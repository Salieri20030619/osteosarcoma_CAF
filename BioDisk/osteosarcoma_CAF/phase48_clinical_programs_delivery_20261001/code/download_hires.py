from pathlib import Path
import os,sys,subprocess,json,re,gzip,hashlib,datetime,traceback,io
Q=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001');P=Q.parent;R=Q/'resources/histology_hires';R.mkdir(exist_ok=True)
if '--start' in sys.argv:
 f=(Q/'logs/hires_controller.log').open('w');p=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(Q/'hires_controller.pid').write_text(str(p.pid));print('HIRES_DIRECT_DOWNLOAD_STARTED',p.pid)
else:
 import requests
 from PIL import Image
 from urllib.parse import urlparse
 import pandas as pd
 session=requests.Session();session.trust_env=False
 for k in list(os.environ):
  if k.lower().endswith('_proxy'):os.environ.pop(k,None)
 os.environ['NO_PROXY']='*';text=(P/'phase41_full_transcriptome_setup_20261001/resources/spatial_source_check/GEO_full_family.txt').read_text();meta=pd.read_csv(P/'phase36_external_validation_20260927/resources/sample_patient_manifest_corrected.tsv',sep='\t');rows=[];rc=0
 for r in meta.itertuples(index=False):
  lines=[l.split(' = ',1)[1] for l in text.splitlines() if l.startswith('!Sample_supplementary_file') and str(r.gsm) in l and l.endswith('tissue_hires_image.png.gz')];assert len(lines)==1;url=lines[0].replace('ftp://','https://');assert urlparse(url).hostname=='ftp.ncbi.nlm.nih.gov';dest=R/Path(urlparse(url).path).name
  try:
   if not dest.exists():
    with session.get(url,timeout=(20,180),stream=True) as response:
     response.raise_for_status();length=int(response.headers.get('Content-Length',0));tmp=dest.with_suffix(dest.suffix+'.part')
     with tmp.open('wb') as f:
      for chunk in response.iter_content(1024*1024):f.write(chunk)
     assert not length or tmp.stat().st_size==length;tmp.replace(dest)
   content=gzip.decompress(dest.read_bytes());im=Image.open(io.BytesIO(content));im.load();rows.append({'sample':r.sample,'gsm':r.gsm,'url':url,'path':str(dest),'bytes':dest.stat().st_size,'SHA256':hashlib.sha256(dest.read_bytes()).hexdigest(),'width':im.width,'height':im.height,'status':'complete','proxy':False});print('HIRES_COMPLETE',r.sample,dest.stat().st_size,flush=True)
  except Exception as e:rows.append({'sample':r.sample,'url':url,'status':'failed','error':repr(e),'proxy':False});rc=1;traceback.print_exc()
  (Q/'HIRES_DOWNLOAD_STATUS.json').write_text(json.dumps({'state':'running','done':len(rows),'total':len(meta),'outcomes':rows},indent=2))
 (Q/'HIRES_DOWNLOAD_STATUS.json').write_text(json.dumps({'state':'complete' if rc==0 else 'partial_failure','outcomes':rows,'proxy':False},indent=2));(Q/'hires_download.exit_code').write_text(str(rc))
 if rc==0:(Q/'HIRES_DOWNLOAD_COMPLETE').write_text('SUCCESS')
 print('HIRES_DOWNLOAD_EXIT',rc,flush=True)
