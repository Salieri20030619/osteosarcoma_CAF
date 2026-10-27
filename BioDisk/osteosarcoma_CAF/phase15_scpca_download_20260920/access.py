from pathlib import Path
import requests,json,os
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase15_scpca_download_20260920');base='https://api.scpca.alexslemonade.org/v1/';s=requests.Session()
tokenfile=O/'.api_token'
if tokenfile.exists():token=tokenfile.read_text().strip()
else:
    r=s.post(base+'tokens/',json={'email':'YOUR_EMAIL@example.invalid','is_activated':True},timeout=30)
    print('TOKEN_HTTP',r.status_code,flush=True);r.raise_for_status();d=r.json();token=d['id'];tokenfile.write_text(token);os.chmod(tokenfile,0o600)
print('TOKEN_AVAILABLE',bool(token),flush=True)
d=json.loads((O/'samples_public_probe.json').read_text());print('SAMPLE_KEYS',list(d),flush=True)
for k in d:
    if 'dataset' in k or 'computed' in k:print('SAMPLE_DOWNLOAD_FIELD',k,json.dumps(d[k]),flush=True)
datasets=json.loads((O/'ccdl-datasets_public_probe.json').read_text())['results'];results=[]
for project in ['SCPCP000017','SCPCP000023']:
    row=next(x for x in datasets if x.get('ccdl_project_id')==project and x['format']=='ANN_DATA' and not x['ccdl_is_merged'])
    r=s.get(base+'ccdl-datasets/'+row['id']+'/',headers={'API-KEY':token},timeout=30);r.raise_for_status();result=r.json()
    url=result.get('download_url') or (result.get('computed_file') or {}).get('download_url')
    if not url:
        r=s.get(base+'computed-files/'+str(row['computed_file']['id'])+'/',headers={'API-KEY':token},timeout=30);r.raise_for_status();url=r.json().get('download_url')
    if not url:print('NO_DOWNLOAD_URL',project,flush=True);continue
    protected=O/(project+'_signed.json');protected.write_text(json.dumps({'url':url,'dataset_id':row['id'],'size':row['computed_file']['size_in_bytes']}));os.chmod(protected,0o600)
    rr=s.get(url,headers={'Range':'bytes=0-3'},stream=True,timeout=30)
    content=next(rr.iter_content(4),b'');rr.close()
    rec={'project':project,'dataset_id':row['id'],'expected_bytes':row['computed_file']['size_in_bytes'],'HTTP':rr.status_code,'range':rr.headers.get('Content-Range'),'zip_magic':content.hex()};results.append(rec);print('DOWNLOAD_ACCESS',json.dumps(rec),flush=True)
(O/'access_status.json').write_text(json.dumps(results,indent=2))
