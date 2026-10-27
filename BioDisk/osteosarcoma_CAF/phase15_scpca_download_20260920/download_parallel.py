"""Download selected ZIP_STORED members via official HTTP Range support.
Per-chunk SHA256 resume, ZIP CRC32 and final SHA256 verification.
"""
from range_zip import ROOT,TARGETS,RemoteZip,selected
from pathlib import PurePosixPath
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,zipfile,hashlib,zlib,struct,os,time,threading,datetime
import requests,anndata as ad,numpy as np
CHUNK=1024*1024;WORKERS=8;tls=threading.local()
def session():
    if not hasattr(tls,'session'):
        tls.session=requests.Session();tls.session.trust_env=False
        proxy=os.environ.get('SCPCA_DOWNLOAD_PROXY','http://127.0.0.1:7897')
        tls.session.proxies={'http':proxy,'https':proxy}
    return tls.session
state={'status':'DOWNLOADING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples_requested':TARGETS,'completed_files':[],'chunk_bytes':CHUNK,'workers':WORKERS,'proxy':os.environ.get('SCPCA_DOWNLOAD_PROXY','http://127.0.0.1:7897')}
def status():
    p=ROOT/'STATUS.tmp';p.write_text(json.dumps(state,indent=2));p.replace(ROOT/'STATUS.json')
def file_digest(path):
    h=hashlib.sha256();crc=0
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b);crc=zlib.crc32(b,crc)
    return h.hexdigest(),crc&0xffffffff
verified=[];inventory=[];out=ROOT/'files';out.mkdir(exist_ok=True);status()
for project,samples in TARGETS.items():
    remote=RemoteZip(project,block_size=65536)
    with zipfile.ZipFile(remote) as archive:
        for info in sorted([v for v in archive.infolist() if not v.is_dir() and selected(v,project)],key=lambda v:v.header_offset):
            name=PurePosixPath(info.filename)
            assert not name.is_absolute() and '..' not in name.parts
            assert info.compress_type==zipfile.ZIP_STORED and info.file_size==info.compress_size
            dest=out/project/str(name);dest.parent.mkdir(parents=True,exist_ok=True)
            state.update(current_project=project,current_file=info.filename,current_file_total_bytes=info.file_size);status()
            reused=False
            if dest.exists() and dest.stat().st_size==info.file_size:
                digest,crc=file_digest(dest);reused=crc==info.CRC
            if not reused:
                remote.seek(info.header_offset);header=remote.read(30);assert header[:4]==b'PK\x03\x04'
                fname_len,extra_len=struct.unpack_from('<HH',header,26);offset=info.header_offset+30+fname_len+extra_len
                part=dest.with_name(dest.name+'.part');ledger=dest.with_name(dest.name+'.chunks.json');completed={}
                if part.exists() and ledger.exists():
                    old=json.loads(ledger.read_text())
                    if old.get('crc32')==info.CRC and old.get('size')==info.file_size:
                        with part.open('rb') as f:
                            for k,h in old.get('completed',{}).items():
                                start=int(k)*CHUNK;f.seek(start);b=f.read(min(CHUNK,info.file_size-start))
                                if hashlib.sha256(b).hexdigest()==h:completed[k]=h
                fd=os.open(part,os.O_CREAT|os.O_RDWR,0o600);os.ftruncate(fd,info.file_size)
                def worker(k):
                    start=k*CHUNK;n=min(CHUNK,info.file_size-start);lo=offset+start;hi=lo+n-1
                    for attempt in range(5):
                        try:
                            r=session().get(remote.cfg['url'],headers={'Range':f'bytes={lo}-{hi}'},timeout=(15,45))
                            if r.status_code!=206:raise RuntimeError('HTTP '+str(r.status_code))
                            if r.headers.get('Content-Range')!=f'bytes {lo}-{hi}/{remote.size}' or len(r.content)!=n:raise RuntimeError('range mismatch')
                            data=r.content;written=os.pwrite(fd,data,start)
                            if written!=n:raise RuntimeError('short local write')
                            return str(k),hashlib.sha256(data).hexdigest(),n
                        except Exception as e:
                            if attempt==4:raise RuntimeError(f'{project} range {lo}-{hi}: {type(e).__name__}') from None
                            time.sleep(min(2**attempt,10))
                def checkpoint():
                    temp=ledger.with_name(ledger.name+'.tmp');temp.write_text(json.dumps({'crc32':info.CRC,'size':info.file_size,'completed':completed}));temp.replace(ledger)
                    state['current_file_verified_chunks']=len(completed);state['current_file_total_chunks']=(info.file_size+CHUNK-1)//CHUNK;status()
                try:
                    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
                        pending=[pool.submit(worker,k) for k in range((info.file_size+CHUNK-1)//CHUNK) if str(k) not in completed]
                        for task in as_completed(pending):
                            k,h,n=task.result();completed[k]=h
                            if len(completed)%8==0:checkpoint()
                    os.fsync(fd);checkpoint()
                finally:os.close(fd)
                digest,crc=file_digest(part)
                assert crc==info.CRC and part.stat().st_size==info.file_size,'CRC/size verification failed'
                part.replace(dest)
            rec={'project':project,'archive_member':info.filename,'path':str(dest),'bytes':info.file_size,'zip_crc32':info.CRC,'sha256':digest,'reused':reused}
            verified.append(rec);state['completed_files'].append({'project':project,'name':info.filename,'bytes':info.file_size});status();print('FILE_VERIFIED',project,info.filename,info.file_size,flush=True)
            if dest.suffix=='.h5ad':
                a=ad.read_h5ad(dest,backed='r');sample=next(s for s in samples if s in name.parts)
                rec2={'project':project,'sample':sample,'file':info.filename,'cells':a.n_obs,'genes':a.n_vars,'layers':list(a.layers.keys()),'obs_columns':list(a.obs.columns),'var_columns':list(a.var.columns),'read_check':'PASS'}
                matrix=a.layers['counts'] if 'counts' in a.layers else a.X
                ix=np.unique(np.linspace(0,a.n_obs-1,min(50,a.n_obs)).astype(int));chunk=matrix[ix,:];values=chunk.data if hasattr(chunk,'toarray') else np.asarray(chunk)
                assert np.isfinite(values).all() and (values>=0).all();rec2.update(matrix_checked='counts' if 'counts' in a.layers else 'X',sampled_rows=len(ix),sample_values_integer=bool(np.allclose(values,np.round(values),atol=1e-6)))
                a.file.close();inventory.append(rec2);print('H5AD_READABLE',sample,rec2['cells'],rec2['genes'],flush=True)
            (ROOT/'verified_files.json').write_text(json.dumps(verified,indent=2));(ROOT/'matrix_inventory.json').write_text(json.dumps(inventory,indent=2))
for project,samples in TARGETS.items():
    for sample in samples:
        for kind in ['unfiltered','filtered','processed']:
            assert sum(r['sample']==sample and '_'+kind+'_' in r['file'] for r in inventory)==1
state.update(status='DOWNLOAD_COMPLETE_VERIFIED',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),matrix_files=len(inventory),verified_file_count=len(verified),samples=9,total_bytes=sum(r['bytes'] for r in verified));status();(ROOT/'DOWNLOAD_COMPLETE').write_text('DOWNLOAD_COMPLETE_VERIFIED\n');print('DOWNLOAD_COMPLETE',len(verified),len(inventory),state['total_bytes'],flush=True)
