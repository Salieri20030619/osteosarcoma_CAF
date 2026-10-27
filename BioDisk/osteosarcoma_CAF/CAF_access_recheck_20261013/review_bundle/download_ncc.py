from pathlib import Path
import requests,json,hashlib,concurrent.futures,datetime,zipfile
ROOT=Path(__file__).resolve().parent
tree=json.loads((ROOT/'evidence/ncc_repo_tree_proxy.txt').read_text())
OUT=ROOT/'downloads/NCC';OUT.mkdir(parents=True,exist_ok=True)
def download(f):
 url='https://raw.githubusercontent.com/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/'+tree['sha']+'/'+requests.utils.quote(f['path'])
 p=OUT/f['path'];s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
 if not p.exists() or p.stat().st_size!=f['size']:
  for attempt in range(4):
   try:
    s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'} if attempt%2==0 else {}
    with s.get(url,stream=True,timeout=(15,60)) as r:
     r.raise_for_status()
     with p.with_suffix(p.suffix+'.part').open('wb') as out:
      for chunk in r.iter_content(1024*1024):out.write(chunk)
    p.with_suffix(p.suffix+'.part').replace(p)
    break
   except requests.RequestException as e:
    print(json.dumps({'file':f['path'],'attempt':attempt+1,'error':str(e)}),flush=True)
    if attempt==3:raise
 sha=hashlib.sha256();git=hashlib.sha1(('blob '+str(f['size'])+'\0').encode())
 with p.open('rb') as inp:
  for c in iter(lambda:inp.read(1024*1024),b''):sha.update(c);git.update(c)
 assert p.stat().st_size==f['size'] and git.hexdigest()==f['sha'],f['path']
 row={'file':f['path'],'bytes':p.stat().st_size,'sha256':sha.hexdigest(),'git_blob_sha1':git.hexdigest(),'git_blob_verified':True,'url':url}
 if p.suffix=='.zip':
  with zipfile.ZipFile(p) as z:
   assert z.testzip() is None
   row['zip_members']=z.namelist()
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(download,[f for f in tree['tree'] if f['type']=='blob']))
(ROOT/'NCC_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'commit':tree['sha'],'files':rows,'bytes':sum(r['bytes'] for r in rows)},indent=2),encoding='utf8')
