from pathlib import Path
from collections import OrderedDict
import io,json,requests,zipfile,time,os
ROOT=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase15_scpca_download_20260920')
TARGETS={'SCPCP000017':['SCPCS000602','SCPCS000603','SCPCS000609','SCPCS000610','SCPCS000611'],
         'SCPCP000023':['SCPCS000794','SCPCS000795','SCPCS000812','SCPCS000813']}
class RemoteZip(io.RawIOBase):
    def __init__(self,project,block_size=4*1024*1024):
        self.project=project;self.cfg=json.loads((ROOT/(project+'_signed.json')).read_text());self.size=self.cfg['size'];self.pos=0;self.block=block_size;self.cache=OrderedDict();self.session=requests.Session();self.transferred=0
        self.session.trust_env=False
        proxy=os.environ.get('SCPCA_DOWNLOAD_PROXY','http://127.0.0.1:7897')
        self.session.proxies={'http':proxy,'https':proxy}
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
        if self.pos<0:raise ValueError('negative seek')
        return self.pos
    def read(self,n=-1):
        n=max(0,self.size-self.pos) if n<0 else min(n,max(0,self.size-self.pos));parts=[]
        while n:
            key=self.pos//self.block;start=key*self.block;end=min(start+self.block,self.size)-1
            if key not in self.cache:
                for attempt in range(4):
                    try:
                        r=self.session.get(self.cfg['url'],headers={'Range':f'bytes={start}-{end}'},timeout=(20,120))
                        if r.status_code!=206:raise RuntimeError(f'HTTP {r.status_code}; project {self.project}; byte range {start}-{end}')
                        expected=f'bytes {start}-{end}/{self.size}'
                        if r.headers.get('Content-Range')!=expected or len(r.content)!=end-start+1:raise RuntimeError('Range integrity mismatch')
                        self.cache[key]=r.content;self.transferred+=len(r.content);break
                    except Exception as exc:
                        if attempt==3:raise RuntimeError(f'Range request failed {self.project} {start}-{end}: {type(exc).__name__}') from None
                        time.sleep(2**attempt)
                while len(self.cache)>2:self.cache.popitem(last=False)
            self.cache.move_to_end(key);offset=self.pos-start;take=min(n,len(self.cache[key])-offset);parts.append(self.cache[key][offset:offset+take]);self.pos+=take;n-=take
        return b''.join(parts)
def selected(info,project):
    name=info.filename
    if any(p.endswith('_bulk') for p in Path(name).parts):return False
    return any(sample in name for sample in TARGETS[project]) or (not any(p.startswith('SCPCS') for p in Path(name).parts) and name.lower().endswith(('.tsv','.md','.txt','.json')))
if __name__=='__main__':
    manifest=[]
    for project in TARGETS:
        f=RemoteZip(project,block_size=65536)
        with zipfile.ZipFile(f) as z:
            for info in z.infolist():
                if not info.is_dir() and selected(info,project):
                    row=dict(project=project,name=info.filename,uncompressed_bytes=info.file_size,compressed_bytes=info.compress_size,crc32=info.CRC);manifest.append(row);print(json.dumps(row),flush=True)
        print('PROJECT_TRANSFER_ESTIMATE',project,sum(r['compressed_bytes'] for r in manifest if r['project']==project),flush=True)
    (ROOT/'selected_file_manifest.json').write_text(json.dumps(manifest,indent=2))
