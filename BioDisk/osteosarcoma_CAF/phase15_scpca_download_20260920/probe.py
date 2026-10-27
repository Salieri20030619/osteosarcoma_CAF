from pathlib import Path
import requests,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase15_scpca_download_20260920');O.mkdir(exist_ok=True)
s=requests.Session()
base='https://api.scpca.alexslemonade.org/v1/'
for resource in ['projects/SCPCP000017/','samples/SCPCS000602/','tokens/','ccdl-datasets/?limit=200']:
    try:
        r=s.options(base+resource,timeout=25) if resource=='tokens/' else s.get(base+resource,timeout=25)
        r.raise_for_status();d=r.json();(O/(resource.split('/')[0]+'_public_probe.json')).write_text(json.dumps(d,indent=2))
        if resource.startswith('ccdl'):
            for x in d.get('results',[]):
                if x.get('ccdl_project_id') in ['SCPCP000017','SCPCP000023']:print('READY',json.dumps({k:x.get(k) for k in ['id','format','ccdl_project_id','ccdl_is_merged','is_succeeded','estimated_size_in_bytes','computed_file']}),flush=True)
        elif resource=='tokens/':print('TOKEN_SCHEMA',json.dumps(d.get('actions',{})),flush=True)
        else:print(resource,json.dumps(d),flush=True)
    except Exception as e:print(type(e).__name__,str(e).split('?')[0],flush=True)
