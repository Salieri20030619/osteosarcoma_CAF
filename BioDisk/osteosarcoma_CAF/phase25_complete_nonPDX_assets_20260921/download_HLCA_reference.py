from pathlib import Path
import json,subprocess,hashlib,sys,datetime,zipfile
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase25_complete_nonPDX_assets_20260921');O=P/'references/HLCA';O.mkdir(parents=True,exist_ok=True)
d=json.loads((P/'HLCA_current_official.json').read_text());selected=['HLCA_reference_model.zip','HLCA_reference_model_gene_order_ids_and_symbols.csv','HLCA_full_v1.1_emb.h5ad'];records=[]
for name in selected:
 f=next(x for x in d['files'] if x['key']==name);dest=O/name;part=O/(name+'.part')
 if not dest.exists():
  with (O/'download.log').open('a') as log:rc=subprocess.run(['curl','-q','--fail','--location','--silent','--show-error','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--retry','3','--retry-delay','10','--continue-at','-','--proxy','http://127.0.0.1:7897','--noproxy','','--output',str(part),f['links']['self']],stdout=log,stderr=subprocess.STDOUT).returncode
  (O/(name+'.exit_code')).write_text(str(rc)+'\n')
  if rc:sys.exit(rc)
  part.replace(dest)
 md5=hashlib.md5();sha=hashlib.sha256()
 with dest.open('rb') as inp:
  for b in iter(lambda:inp.read(4*1024*1024),b''):md5.update(b);sha.update(b)
 assert dest.stat().st_size==f['size'] and 'md5:'+md5.hexdigest()==f['checksum']
 if name.endswith('.zip'):
  with zipfile.ZipFile(dest) as z:assert z.testzip() is None
 records.append({'file':str(dest),'bytes':f['size'],'official_md5':f['checksum'],'sha256':sha.hexdigest(),'source':f['links']['self']});(O/'verified_files.json').write_text(json.dumps(records,indent=2));print('HLCA_FILE_VERIFIED',name,flush=True)
(O/'COMPLETE').write_text('OFFICIAL_HLCA_MAPPING_ASSETS_MD5_VERIFIED\n')
