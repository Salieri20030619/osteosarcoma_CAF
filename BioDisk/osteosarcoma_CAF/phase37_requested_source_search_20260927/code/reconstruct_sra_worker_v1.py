from pathlib import Path
import json,os,time,subprocess,concurrent.futures,re,math,fcntl,ctypes,select,shutil,traceback,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927');D=O/'downloads/GSE237070_SRA';F=O/'downloads/GSE237070_FASTQ';F.mkdir(exist_ok=True);R=O/'resources'
lock=(O/'reconstruct.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
libc=ctypes.CDLL(None,use_errno=True);fd=libc.inotify_init1(os.O_NONBLOCK);assert fd>=0
for folder in [O,D,R]:assert libc.inotify_add_watch(fd,str(folder).encode(),0x8|0x80|0x100)>=0
poll=select.poll();poll.register(fd,select.POLLIN)
state={'pid':os.getpid(),'state':'starting','converted_runs':[],'reconstructed_samples':[],'started':time.time()}
def status(stage,**kw):
 state.update(state=stage,updated=time.time(),**kw);p=O/'RECONSTRUCTION_PROGRESS.json.tmp';p.write_text(json.dumps(state,indent=2));p.replace(O/'RECONSTRUCTION_PROGRESS.json')
def wait_for(path,errors):
 while not path.exists():
  for e in errors:
   if e.exists():raise RuntimeError(f'Prerequisite failed: {e}: {e.read_text()[:1000]}')
  events=poll.poll(900000)
  if events:
   try:os.read(fd,65536)
   except BlockingIOError:pass
def resources(need_memory=False):
 while True:
  mem=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
  if shutil.disk_usage(O).free>=80*2**30 and (not need_memory or mem>=114*2**30):return
  status('waiting_compute_resources',available_memory_bytes=mem);time.sleep(900)
env={k:v for k,v in os.environ.items() if k.lower() not in ['http_proxy','https_proxy','all_proxy']};env['LD_LIBRARY_PATH']='/mnt/BioDisk/wangmingyuan/miniconda3/lib';env['NCBI_SETTINGS']=str(O/'tools/sra_config/user-settings.mkfg')
tool='/mnt/BioDisk/wangmingyuan/miniconda3/pkgs/sra-tools-3.2.1-h4304569_1/bin/fastq-dump';PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
manifest=json.loads((O/'tables/GSE237070_sra_download_manifest.json').read_text());order=json.loads((O/'DOWNLOAD_STARTED.json').read_text())['runs'];lookup={r['run']:r for r in manifest};gorder=list(dict.fromkeys(r['gsm'] for r in manifest))
def convert(acc):
 x=lookup[acc];dest=F/acc;dest.mkdir(exist_ok=True);done=dest/'FASTQ_READY.json'
 if done.exists():return
 wait_for(D/(acc+'.verified.json'),[O/'DOWNLOAD_ERROR.json']);resources();status('converting_SRA_to_FASTQ',current_run=acc)
 n=int(x['read_statistics'][0]['count']);parts=[]
 for i in range(8):
  start=n*i//8+1;end=n*(i+1)//8
  if end>=start:parts.append((i,start,end))
 def one(part):
  i,start,end=part;out=dest/f'part{i:02d}';out.mkdir(exist_ok=True);marker=out/'PART_READY.json'
  if marker.exists():return json.loads(marker.read_text())
  cmd=[tool,'--split-files','--gzip','--origfmt','-N',str(start),'-X',str(end),'--outdir',str(out),str(D/(acc+'.sra'))]
  with (out/'fastq_dump.log').open('w') as log:p=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
  text=(out/'fastq_dump.log').read_text();written=re.findall(r'Written\s+([0-9,]+)\s+spots',text);assert p.returncode==0 and written and int(written[-1].replace(',',''))==end-start+1,(acc,i,p.returncode,text)
  files=[out/(acc+'_'+str(j)+'.fastq.gz') for j in [1,2,3]];assert all(f.exists() and f.stat().st_size>0 for f in files)
  with (out/'gzip_integrity.log').open('w') as log:q=subprocess.run(['gzip','-t']+[str(f) for f in files],stdout=log,stderr=subprocess.STDOUT)
  assert q.returncode==0,(acc,i,'gzip integrity')
  row={'part':i,'start':start,'end':end,'spots':end-start+1,'fastq_dump_exit':p.returncode,'gzip_integrity_exit':q.returncode,'index':str(files[0]),'barcode_UMI':str(files[1]),'RNA':str(files[2])};marker.write_text(json.dumps(row,indent=2));return row
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:result=list(pool.map(one,parts))
 assert sum(x['spots'] for x in result)==n
 done.write_text(json.dumps({'run':acc,'spots':n,'parts':result,'source_MD5_passed':True},indent=2));state['converted_runs'].append(acc)
def align(gsm):
 rows=[r for r in manifest if r['gsm']==gsm]
 if not all((F/r['run']/'FASTQ_READY.json').exists() for r in rows):return
 d=O/'objects/GSE237070_reconstructed'/gsm;d.mkdir(parents=True,exist_ok=True)
 if (d/'MATRIX_QA.json').exists():return
 wait_for(R/'REFERENCE_READY.json',[O/'REFERENCE_ERROR.json']);resources(True);status('aligning_and_counting',current_sample=gsm)
 parts=[p for r in rows for p in json.loads((F/r['run']/'FASTQ_READY.json').read_text())['parts']];out=d/'STAR';out.mkdir(exist_ok=True)
 cmd=[str(O/'tools/STAR_2.7.11b'),'--genomeDir',str(R/'STAR_2.7.11b_GRCh38_2020A'),'--runThreadN','16','--readFilesIn',','.join(p['RNA'] for p in parts),','.join(p['barcode_UMI'] for p in parts),'--readFilesCommand','zcat','--soloType','CB_UMI_Simple','--soloCBwhitelist',str(R/'10x_version3_whitelist.txt'),'--soloCBstart','1','--soloCBlen','16','--soloUMIstart','17','--soloUMIlen','12','--soloStrand','Forward','--soloFeatures','Gene','--soloCellFilter','EmptyDrops_CR','--outSAMtype','None','--outFileNamePrefix',str(out)+'/']
 (d/'RECONSTRUCTION_METHOD.json').write_text(json.dumps({'command':cmd,'software_native_cell_calling':'STARsolo EmptyDrops_CR with default native parameters; no force-cell count or project UMI/gene/MT gate','counting':'Gene (exonic); 10x v3 CB16/UMI12; SRA read2 barcode,read3 RNA','raw_and_native_called_matrices_preserved':True,'source_runs':[r['run'] for r in rows],'biological_unit':'patient','CAF_annotation':'pending'},indent=2))
 with (d/'STAR.log').open('w') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 (d/'STAR.exit_code').write_text(str(p.returncode)+'\n');assert p.returncode==0,(gsm,p.returncode)
 with (d/'matrix_QA.log').open('w') as log:q=subprocess.run([PY,str(O/'code/reconstructed_matrix_qa_v1.py'),gsm],stdout=log,stderr=subprocess.STDOUT)
 (d/'matrix_QA.exit_code').write_text(str(q.returncode)+'\n');assert q.returncode==0,(gsm,'matrixQA',q.returncode);state['reconstructed_samples'].append(gsm)
rc=0
try:
 status('waiting_verified_download')
 for acc in order:
  convert(acc)
  for gsm in gorder:align(gsm)
  status('waiting_next_verified_download')
 assert all((O/'objects/GSE237070_reconstructed'/gsm/'MATRIX_QA.json').exists() for gsm in gorder)
 (O/'RECONSTRUCTION_FINISHED.json').write_text(json.dumps({'samples':gorder,'run_count':len(order),'time':time.time(),'CAF_biological_validation':'pending downstream annotation and patient-level evaluation'},indent=2));status('matrices_ready_biological_validation_pending')
except Exception as e:
 rc=1;status('error',error=repr(e));(O/'RECONSTRUCTION_ERROR.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc()},indent=2));print(traceback.format_exc(),flush=True)
finally:os.close(fd);(O/'reconstruction.exit_code').write_text(str(rc)+'\n')
sys.exit(rc)
