"""Download official GTF references with validated HTTP ranges.
SYNTAX_CHECK: python -m py_compile download_coordinates.py
RUN_COMMAND: python download_coordinates.py
EXPECTED_SUCCESS_MARKER: REFERENCE_DOWNLOAD_COMPLETE
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import urllib.request,gzip,os,signal,subprocess,time,json,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914')
R=P/'resources'; log=[]
(P/'logs/coordinate_download.pid').write_text(str(os.getpid()))
refs=[('84','GRCh38'),('75','GRCh37')]
tasks=[]
for rel,build in refs:
    name=f'Homo_sapiens.{build}.{rel}.gtf.gz';url=f'https://ftp.ensembl.org/pub/release-{rel}/gtf/homo_sapiens/{name}'
    with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':'bytes=0-0'}),timeout=60) as q:
        cr=q.headers.get('Content-Range','');assert q.status==206 and cr.startswith('bytes 0-0/'),(q.status,cr)
        size=int(cr.split('/')[-1]);assert len(q.read())==1
    # Only stop our earlier slow curl downloads after range support has been verified.
    processes=subprocess.check_output(['ps','-u',str(os.getuid()),'-o','pid=,comm=,args='],text=True)
    for line in processes.splitlines():
        v=line.strip().split(None,2)
        if len(v)==3 and v[1]=='curl' and str(R/(name+'.part')) in v[2]:
            os.kill(int(v[0]),signal.SIGTERM);log.append({'stopped_previous_task_curl_pid':int(v[0]),'name':name})
    seg=R/(name+'.segments');seg.mkdir(exist_ok=True)
    for k,start in enumerate(range(0,size,2*1024*1024)):
        end=min(start+2*1024*1024,size)-1;tasks.append((name,url,size,k,start,end,seg/f'{k:03}.part'))
def fetch(t):
    name,url,size,k,start,end,p=t;need=end-start+1
    if p.exists() and p.stat().st_size==need:return t
    last=None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'}),timeout=90) as q:
                assert q.status==206 and q.headers.get('Content-Range')==f'bytes {start}-{end}/{size}'
                data=q.read();assert len(data)==need
            p.write_bytes(data);print('SEGMENT_READY',name,k,len(data),flush=True);return t
        except Exception as e:last=e
    raise RuntimeError((name,k,str(last)))
with ThreadPoolExecutor(max_workers=12) as ex:
    for fut in as_completed([ex.submit(fetch,t) for t in tasks]):fut.result()
for name in sorted({t[0] for t in tasks}):
    tt=sorted((t for t in tasks if t[0]==name),key=lambda t:t[3]);dst=R/name;tmp=R/(name+'.assembling')
    with tmp.open('wb') as f:
        for t in tt:f.write(t[-1].read_bytes())
    assert tmp.stat().st_size==tt[0][2]
    with gzip.open(tmp,'rb') as f:
        while f.read(1024*1024):pass
    tmp.replace(dst);log.append({'url':tt[0][1],'file':name,'bytes':dst.stat().st_size,'sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
(R/'reference_download.json').write_text(json.dumps(log,indent=2))
(P/'logs/REFERENCE_DOWNLOAD_COMPLETE').write_text(time.strftime('%Y-%m-%dT%H:%M:%S'))
print('REFERENCE_DOWNLOAD_COMPLETE')
