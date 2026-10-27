from pathlib import Path
import json, os, requests

old = Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase15_scpca_download_20260920')
new = old / 'additional23'
new.mkdir(exist_ok=True)
code = (old / 'range_zip.py').read_text()
start = code.index('TARGETS=')
end = code.index('class RemoteZip', start)
code = code[:start] + "TARGETS={'SCPCP000023':['SCPCS000785','SCPCS000786','SCPCS000789']}\n" + code[end:]
code = code.replace(str(old), str(new))
(new / 'range_zip.py').write_text(code)
code = (old / 'download_parallel.py').read_text().replace('samples=9,', 'samples=sum(map(len,TARGETS.values())),')
(new / 'download_parallel.py').write_text(code)
code = (old / 'verify_download.py').read_text().replace('len(inv)==27', 'len(inv)==3*sum(map(len,TARGETS.values()))').replace("'samples':9", "'samples':sum(map(len,TARGETS.values()))").replace("'expression_objects':27", "'expression_objects':len(inv)")
(new / 'verify_download.py').write_text(code)
cfg = json.loads((old / 'SCPCP000023_signed.json').read_text())
session = requests.Session()
session.trust_env = False
session.proxies = {'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
r = session.get('https://api.scpca.alexslemonade.org/v1/ccdl-datasets/'+cfg['dataset_id']+'/', headers={'API-KEY':(old/'.api_token').read_text().strip()}, timeout=45)
if r.status_code != 200:
    raise RuntimeError('Official API HTTP '+str(r.status_code))
d = r.json()
url = d.get('download_url') or (d.get('computed_file') or {}).get('download_url')
assert url, 'No download URL returned'
cfg['url'] = url
protected = new / 'SCPCP000023_signed.json'
protected.write_text(json.dumps(cfg)); os.chmod(protected,0o600)
script = (old / 'run_download.sh').read_text().replace(str(old),str(new))
script = script.replace('exit "$code"', 'if [ "$code" -ne 0 ]; then exit "$code"; fi\n"$py" verify_download.py > verification.log 2>&1\nqa=$?\nprintf \'%s\\n\' "$qa" > qa_exit_code\nexit "$qa"')
(new / 'run_download.sh').write_text(script)
print('READY', str(new), '3 samples; explicit proxy; fresh official URL')
