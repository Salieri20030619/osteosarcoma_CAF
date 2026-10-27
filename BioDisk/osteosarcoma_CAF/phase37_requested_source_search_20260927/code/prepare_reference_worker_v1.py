from pathlib import Path
import os,json,time,hashlib,concurrent.futures,requests,threading,tarfile,gzip,subprocess,traceback,sys,fcntl,shutil
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927');R=O/'resources';R.mkdir(exist_ok=True)
lk=(O/'reference.lock').open('w');fcntl.flock(lk,fcntl.LOCK_EX|fcntl.LOCK_NB)
URL='https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz';SIZE=11448662640;MD5='dfd654de39bff23917471e7fcc7a00cd';CHUNK=8*1024*1024
state={'state':'starting','pid':os.getpid(),'source_url':URL,'expected_bytes':SIZE,'completed_ranges':0,'time_started':time.time()};sl=threading.Lock();stop=threading.Event()
def save():
 with sl:
  p=O/'REFERENCE_PROGRESS.json.tmp';p.write_text(json.dumps(state,indent=2));p.replace(O/'REFERENCE_PROGRESS.json')
def monitor():
 while not stop.wait(900):save()
threading.Thread(target=monitor,daemon=True).start();save()
code=0
try:
 source=Path('/mnt/BioDisk/wangmingyuan/miniconda3/pkgs/ngs-tools-1.8.6-pyhdfd78af_0/site-packages/ngs_tools/chemistry/whitelists/10x_version3_whitelist.txt.gz')
 source_sha=hashlib.sha256(source.read_bytes()).hexdigest();pkg=source.parents[4]
 paths=json.loads((pkg/'info/paths.json').read_text())['paths'];expected=next(x['sha256'] for x in paths if x['_path']==str(source.relative_to(pkg)));assert source_sha==expected
 w=R/'10x_version3_whitelist.txt';n=0
 with gzip.open(source,'rt') as f,w.open('w') as out:
  for line in f:
   s=line.strip();assert len(s)==16 and set(s)<=set('ACGT');out.write(s+'\n');n+=1
 assert n==6794880,n
 (R/'whitelist_QA.json').write_text(json.dumps({'source':str(source),'source_package_sha256_passed':True,'sha256_gzip':source_sha,'barcodes':n,'format':'16-base kit whitelist, not a project expression/QC threshold'},indent=2))
 dest=R/'refdata-gex-GRCh38-2020-A.tar.gz';part=R/(dest.name+'.partial');journal=R/'reference.completed_ranges.jsonl'
 if not dest.exists():
  done=set(json.loads(l)['chunk'] for l in journal.read_text().splitlines()) if journal.exists() else set()
  if not part.exists():assert not done;part.touch();os.truncate(part,SIZE)
  assert part.stat().st_size==SIZE;fd=os.open(part,os.O_RDWR);j=journal.open('a');jl=threading.Lock()
  state.update(state='downloading_reference',completed_ranges=len(done),total_ranges=(SIZE+CHUNK-1)//CHUNK);save()
  def one(i):
   start=i*CHUNK;end=min(SIZE,start+CHUNK)-1;errs=[]
   for attempt in range(4):
    s=requests.Session();s.trust_env=False;pos=start;t=time.time()
    try:
     with s.get(URL,headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(15,60)) as response:
      assert response.status_code==206;assert response.headers.get('Content-Range')==f'bytes {start}-{end}/{SIZE}'
      for b in response.iter_content(65536):
       assert pos+len(b)<=end+1;assert os.pwrite(fd,b,pos)==len(b);pos+=len(b)
       if time.time()-t>1800:raise TimeoutError('reference range exceeded30minutes')
      assert pos==end+1
     with jl:os.fsync(fd);j.write(json.dumps({'chunk':i})+'\n');j.flush();os.fsync(j.fileno())
     with sl:state['completed_ranges']+=1
     return
    except Exception as e:errs.append(repr(e));time.sleep(2**attempt)
    finally:s.close()
   raise RuntimeError((i,errs))
  try:
   with concurrent.futures.ThreadPoolExecutor(max_workers=64) as pool:list(pool.map(one,[i for i in range((SIZE+CHUNK-1)//CHUNK) if i not in done]))
  finally:j.close();os.close(fd)
  state['state']='checking_reference_md5';save();h=hashlib.md5()
  with part.open('rb') as f:
   for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
  assert h.hexdigest()==MD5,h.hexdigest();part.replace(dest)
 else:
  h=hashlib.md5()
  with dest.open('rb') as f:
   for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
  assert h.hexdigest()==MD5
 state.update(state='extracting_reference',actual_md5=MD5);save()
 needed=['refdata-gex-GRCh38-2020-A/fasta/genome.fa','refdata-gex-GRCh38-2020-A/genes/genes.gtf','refdata-gex-GRCh38-2020-A/reference.json'];found=[]
 with tarfile.open(dest,'r:gz') as tar:
  for member in tar:
   name=member.name.removeprefix('./')
   if name in needed:
    target=R/name;assert target.resolve().is_relative_to(R.resolve());assert member.isfile();target.parent.mkdir(parents=True,exist_ok=True)
    with tar.extractfile(member) as inp,target.open('wb') as out:
     for b in iter(lambda:inp.read(1024*1024),b''):out.write(b)
    found.append(name)
 assert set(found)==set(needed),found
 ref=R/'refdata-gex-GRCh38-2020-A';idx=R/'STAR_2.7.11b_GRCh38_2020A';idx.mkdir(exist_ok=True)
 while True:
  mem=dict((l.split(':')[0],int(l.split()[1])*1024) for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:'))['MemAvailable']
  if mem>=114*2**30 and shutil.disk_usage(R).free>=80*2**30:break
  state.update(state='waiting_compute_resources',available_memory_bytes=mem);save();time.sleep(900)
 state['state']='building_STAR_index';save()
 cmd=[str(O/'tools/STAR_2.7.11b'),'--runMode','genomeGenerate','--runThreadN','16','--genomeDir',str(idx),'--genomeFastaFiles',str(ref/'fasta/genome.fa'),'--sjdbGTFfile',str(ref/'genes/genes.gtf'),'--sjdbOverhang','150','--limitGenomeGenerateRAM','64000000000','--outFileNamePrefix',str(O/'logs/STAR_index_')]
 with (O/'logs/STAR_index.log').open('w') as out:r=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT)
 assert r.returncode==0,r.returncode
 assert all((idx/n).stat().st_size>0 for n in ['Genome','SA','SAindex','genomeParameters.txt'])
 record={'state':'PASS','source_url':URL,'archive_bytes':SIZE,'archive_md5':MD5,'md5_crosscheck_record':'https://zenodo.org/records/6980009','barcode_package':'ngs-tools1.8.6','barcode_count':n,'STAR_version':'2.7.11b','index_command':cmd,'STAR_native_exit':r.returncode,'time_completed':time.time()}
 (R/'REFERENCE_READY.json').write_text(json.dumps(record,indent=2));state.update(state='complete',completed=time.time())
except Exception as e:
 code=1;state.update(state='error',error=repr(e));(O/'REFERENCE_ERROR.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc()},indent=2));print(traceback.format_exc(),flush=True)
finally:stop.set();save();(O/'reference.exit_code').write_text(str(code)+'\n')
sys.exit(code)
