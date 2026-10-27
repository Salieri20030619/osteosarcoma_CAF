from pathlib import Path
import os, json, hashlib, datetime, sys, traceback, subprocess, time
import requests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase41_full_transcriptome_setup_20261001'
for sub in ['resources','logs','code','tables','objects','delivery']: (O/sub).mkdir(parents=True,exist_ok=True)
for key in list(os.environ):
 if key.lower().endswith('_proxy'): os.environ.pop(key,None)
os.environ['NO_PROXY']='*'
S=requests.Session(); S.trust_env=False
base='https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/'
resources=[
 ('hallmark.gmt',base+'h.all.v2025.1.Hs.symbols.gmt','gmt'),
 ('reactome.gmt',base+'c2.cp.reactome.v2025.1.Hs.symbols.gmt','gmt'),
 ('go_bp.gmt',base+'c5.go.bp.v2025.1.Hs.symbols.gmt','gmt'),
 ('curated_cgp.gmt',base+'c2.cgp.v2025.1.Hs.symbols.gmt','gmt'),
 ('tf_targets.gmt',base+'c3.tft.v2025.1.Hs.symbols.gmt','gmt'),
 ('hgnc_complete_set.txt','https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt','hgnc'),
 ('collectri.tsv','https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&fields=sources,references&organisms=9606&format=tsv&license=commercial','network'),
 ('ligand_receptor.tsv','https://omnipathdb.org/interactions?datasets=ligrecextra&genesymbols=yes&fields=sources,references&organisms=9606&format=tsv&license=commercial','network')]
state={'pid':os.getpid(),'started':datetime.datetime.now().astimezone().isoformat(),'proxy_used':False,'requests_trust_env':False,'resources':[]}
def write():
 f=O/'DOWNLOAD_STATUS.json'; tmp=f.with_suffix('.tmp'); tmp.write_text(json.dumps(state,indent=2)); tmp.replace(f)
rc=0
try:
 for name,url,kind in resources:
  entry={'file':name,'url':url,'state':'downloading','proxy':False};state['resources'].append(entry);write()
  f=O/'resources'/name
  try:
   if not f.exists():
    last=None
    for attempt in range(3):
     try:
      part=f.with_suffix(f.suffix+'.part')
      with S.get(url,stream=True,timeout=(20,90)) as r:
       r.raise_for_status(); size=0
       with part.open('wb') as h:
        for chunk in r.iter_content(1024*1024):
         h.write(chunk);size+=len(chunk)
         entry['bytes_received']=size;write()
       entry['http_status']=r.status_code;entry['resolved_url']=r.url;entry['etag']=r.headers.get('ETag');entry['last_modified']=r.headers.get('Last-Modified')
      lines=part.read_text().splitlines()
      if kind=='gmt': assert len(lines)>1 and all(len(x.split('\t'))>=3 for x in lines if x.strip()), 'invalid GMT'
      elif kind=='hgnc': assert 'symbol' in lines[0] and 'ensembl_gene_id' in lines[0]
      else: assert 'source_genesymbol' in lines[0] and 'target_genesymbol' in lines[0] and len(lines)>1, 'invalid network response'
      part.replace(f);last=None;break
     except Exception as e:
      last=e;entry['last_error']=repr(e);entry['attempt']=attempt+1;write()
      if attempt<2:time.sleep(10*(attempt+1))
    if last:raise last
   entry.update(state='verified',bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),rows=sum(1 for _ in f.open()))
   print('VERIFIED',name,entry['bytes'],flush=True)
  except Exception as e:
   entry.update(state='failed',error=repr(e));rc=1;print('FAILED',name,repr(e),flush=True)
  write()
 state['state']='complete' if rc==0 else 'complete_with_resource_failures'
 state['completed']=datetime.datetime.now().astimezone().isoformat();write()
 if rc==0:(O/'DOWNLOAD_COMPLETE').write_text('SUCCESS\n')
 else:(O/'DOWNLOAD_FAILURES_RECORDED').write_text('SEE_DOWNLOAD_STATUS\n')
finally:
 (O/'download.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)