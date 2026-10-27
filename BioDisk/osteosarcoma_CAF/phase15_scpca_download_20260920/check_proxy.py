from pathlib import Path
import os,json,time,requests,concurrent.futures,urllib.parse
root=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase15_scpca_download_20260920')
def safe_proxy(v):
    u=urllib.parse.urlsplit(v);return f'{u.scheme}://{u.hostname}:{u.port}'
for path in Path('/proc').glob('[0-9]*'):
    try:
        if path.stat().st_uid!=os.getuid():continue
        cmd=(path/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
        if 'download_parallel.py' not in cmd or '/bin/python' not in cmd:continue
        env=dict(v.split(b'=',1) for v in (path/'environ').read_bytes().split(b'\0') if b'=' in v)
        proxies={k.decode():safe_proxy(v.decode()) for k,v in env.items() if k.decode().lower() in ['http_proxy','https_proxy','all_proxy']}
        print('DOWNLOAD_PROCESS',path.name,'PROXY_ENV',json.dumps(proxies),flush=True)
    except (OSError,ValueError):pass
cfg=json.loads((root/'SCPCP000017_signed.json').read_text());url=cfg['url'];size=1024*1024
def check(mode):
    s=requests.Session();s.trust_env=False
    if mode=='proxy':s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
    start=time.monotonic();n=0
    try:
        with s.get(url,headers={'Range':f'bytes=0-{size-1}'},stream=True,timeout=(6,12)) as r:
            code=r.status_code
            if code!=206:return {'mode':mode,'HTTP':code}
            for block in r.iter_content(65536):
                n+=len(block)
                if time.monotonic()-start>18:break
        elapsed=time.monotonic()-start
        return {'mode':mode,'HTTP':code,'bytes':n,'seconds':round(elapsed,3),'MiB_per_second':round(n/1048576/elapsed,3),'full_1MiB':n==size}
    except Exception as e:return {'mode':mode,'bytes':n,'seconds':round(time.monotonic()-start,3),'error':type(e).__name__}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:result=list(pool.map(check,['direct','proxy']))
(root/'proxy_benchmark.json').write_text(json.dumps(result,indent=2))
for row in result:print('BENCHMARK',json.dumps(row),flush=True)
print('STATUS',json.loads((root/'STATUS.json').read_text()).get('current_file'))
