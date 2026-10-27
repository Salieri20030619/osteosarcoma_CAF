from pathlib import Path
import requests,json,hashlib,gzip,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase19_external_validation_20260921');(O/'inputs').mkdir(parents=True,exist_ok=True)
urls={'GSE39055_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39055/matrix/GSE39055_series_matrix.txt.gz','GPL14951_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14951/soft/GPL14951_family.soft.gz'}
records=[]
for name,url in urls.items():
 dest=O/'inputs'/name;errors=[]
 if not dest.exists():
  for proxy in [None,'http://127.0.0.1:7897']:
   s=requests.Session();s.trust_env=False
   if proxy:s.proxies={'https':proxy,'http':proxy}
   try:
    r=s.get(url,timeout=(15,90));r.raise_for_status();gzip.decompress(r.content);dest.write_bytes(r.content);break
   except Exception as e:errors.append(type(e).__name__)
  else:raise RuntimeError(name+' download failed '+str(errors))
 with gzip.open(dest,'rt',errors='replace') as f:
  text=f.read()
 assert '!series_matrix_table_end' in text if 'matrix' in name else '!platform_table_end' in text
 records.append({'file':name,'url':url,'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'complete_marker_checked':True,'downloaded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 if 'matrix' in name:
  print('\n'.join(x for x in text.splitlines() if x.startswith(('!Sample_characteristics_ch1','!Sample_data_processing','!Sample_title'))),flush=True)
 else:
  lines=text.splitlines();i=lines.index('!platform_table_begin');print('PLATFORM_COLUMNS',lines[i+1],flush=True)
(O/'inputs/download_manifest.json').write_text(json.dumps(records,indent=2));print('DOWNLOAD_COMPLETE',flush=True)
