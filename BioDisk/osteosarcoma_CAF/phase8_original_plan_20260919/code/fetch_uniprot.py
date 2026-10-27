from pathlib import Path
import requests,json,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase8_original_plan_20260919')
genes=json.loads((O/'resources/surface_candidates_frozen.json').read_text())['surface_candidates']+['S100A4','ACTA2','TAGLN']
query='('+ ' OR '.join('gene_exact:'+g for g in genes)+') AND organism_id:9606 AND reviewed:true'
url='https://rest.uniprot.org/uniprotkb/search'
r=requests.get(url,params={'query':query,'format':'json','size':100},timeout=(20,120));r.raise_for_status()
data=r.json();assert data.get('results')
(O/'resources/uniprot_surface_evidence.json').write_text(json.dumps(data,indent=2))
(O/'resources/uniprot_source.json').write_text(json.dumps({'url':r.url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'requested_genes':genes},indent=2))
print('UNIPROT_DOWNLOADED',len(data['results']),flush=True)
