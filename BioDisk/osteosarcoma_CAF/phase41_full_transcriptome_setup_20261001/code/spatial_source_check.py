from pathlib import Path
import os,json,requests,datetime,re,hashlib,urllib.parse,traceback,subprocess,sys
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase41_full_transcriptome_setup_20261001';R=O/'resources/spatial_source_check';R.mkdir(exist_ok=True)
if '--start' in sys.argv:
 assert not (O/'spatial_source_check.pid').exists(),'already launched'
 env=os.environ.copy()
 for k in list(env):
  if k.lower().endswith('_proxy'):env.pop(k,None)
 env['NO_PROXY']='*'
 with (O/'logs/spatial_source_check.log').open('x') as log:
  p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',__file__],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
 (O/'spatial_source_check.pid').write_text(str(p.pid)+'\n');print('DIRECT_BACKGROUND_SPATIAL_SOURCE_CHECK_STARTED',p.pid);raise SystemExit(0)
for k in list(os.environ):
 if k.lower().endswith('_proxy'):os.environ.pop(k,None)
s=requests.Session();s.trust_env=False;out=[];rc=0
urls=[('GEO_supplement_directory.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE293nnn/GSE293065/suppl/'),('GEO_full_family.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293065&targ=all&form=text&view=full')]
try:
 for name,url in urls:
  entry={'file':name,'url':url,'proxy_used':False}
  try:
   r=s.get(url,timeout=(20,90));r.raise_for_status();f=R/name;f.write_bytes(r.content);entry.update(state='verified',bytes=f.stat().st_size,sha256=hashlib.sha256(r.content).hexdigest(),http_status=r.status_code)
   if 'directory' in name:
    hrefs=re.findall(r'href=[\"\x27]([^\"\x27]+)',r.text);candidates=[]
    for href in hrefs:
     u=urllib.parse.urljoin(url,href);p=urllib.parse.urlparse(u)
     if p.hostname!='ftp.ncbi.nlm.nih.gov' or not p.path.startswith('/geo/series/GSE293nnn/GSE293065/suppl/'):continue
     n=Path(p.path).name
     if re.search(r'\.(png|jpe?g|tiff?|pdf|tar|gz|zip)$',n,re.I):candidates.append({'file':n,'url':u,'image_or_document':bool(re.search(r'\.(png|jpe?g|tiff?|pdf)$',n,re.I))})
    entry['supplement_candidates']=candidates
   elif 'recaptcha' in r.text.lower():entry['state']='blocked_captcha_saved';rc=1
   print(entry,flush=True)
  except Exception as e:entry.update(state='failed',error=repr(e));rc=1;print(entry,flush=True)
  out.append(entry);(O/'SPATIAL_SOURCE_CHECK_STATUS.json').write_text(json.dumps({'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),'proxy_used':False,'outcomes':out,'state':'checking'},indent=2))
 (O/'SPATIAL_SOURCE_CHECK_STATUS.json').write_text(json.dumps({'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),'proxy_used':False,'outcomes':out,'state':'complete' if rc==0 else 'completed_with_source_access_failures','image_boundary_not_yet_accepted':True,'large_archives_not_downloaded_without_necessary_content_identification':True},indent=2))
 (O/'SPATIAL_SOURCE_CHECK_FINISHED').write_text('FINISHED\n')
finally:(O/'spatial_source_check.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)