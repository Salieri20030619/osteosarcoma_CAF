
from pathlib import Path
import subprocess,hashlib,json,gzip,datetime,sys,shutil
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1');S=O.parent/'search_patient_expansion_20261012_v1'
files=[('barcodes','GSM7453693_barcodes.tsv.gz',68105),('features','GSM7453693_features.tsv.gz',318733),('matrix','GSM7453693_matrix.mtx.gz',126169981)]
records=[]
try:
 for kind,name,size in files:
  dest=O/'downloads/GSE234187'/name;url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7453nnn/GSM7453693/suppl/'+name
  existing=S/'sources'/('GSE234187_'+kind+'.gz')
  if kind!='matrix' and existing.exists():shutil.copy2(existing,dest);method='reuse_verified_prior_download'
  else:
   assert shutil.disk_usage(O).free>size+5_000_000_000,'Insufficient free disk'
   part=dest.with_name(dest.name+'.part');method='new_download'
   print('DOWNLOADING',name,flush=True)
   r=subprocess.run(['curl','-fL','--connect-timeout','15','--max-time','1800','--retry','2','--retry-delay','5','--continue-at','-','-sS','-o',str(part),url])
   assert r.returncode==0,f'download exit {r.returncode}'
   assert part.stat().st_size==size,(part.stat().st_size,size);part.rename(dest)
  assert dest.stat().st_size==size
  h=hashlib.sha256()
  with dest.open('rb') as f:
   for b in iter(lambda:f.read(8388608),b''):h.update(b)
  with gzip.open(dest,'rb') as f:
   header=f.read(1024)
   for b in iter(lambda:f.read(8388608),b''):pass
  records.append(dict(file=str(dest),url=url,bytes=size,sha256=h.hexdigest(),gzip_integrity='PASS',method=method))
  print('VERIFIED',name,size,flush=True)
 (O/'resources/GSE234187_DOWNLOAD_COMPLETE.json').write_text(json.dumps(dict(status='complete_verified',completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=records),indent=2))
 (O/'logs/GSE234187_download.exit_code').write_text('0\n')
except BaseException as e:
 (O/'resources/GSE234187_DOWNLOAD_FAILED.json').write_text(json.dumps({'error':str(e),'completed_files':records},indent=2));(O/'logs/GSE234187_download.exit_code').write_text('1\n');raise
