from pathlib import Path
import json,urllib.request,concurrent.futures
p=Path('geo_cnv_20260918/sources');tree=json.loads((p/'author_repository_tree.json').read_text());sha=tree['sha']
paths=['06_major_celltype_annotation/01_cell_annotation_integrated_object.ipynb','06_major_celltype_annotation/02_cell_annotation_integrated_object_eda.ipynb','07_cnv_detection/CopyKat/01_sample_list_copykat_integration.r','07_cnv_detection/inferCNV/02_infercnv_sample_list.r']
def get(path):
 data=urllib.request.urlopen('https://raw.githubusercontent.com/khanhnguyendok/osteosarcoma/'+sha+'/'+path,timeout=60).read();f=p/Path(path).name;f.write_bytes(data)
 if path.endswith('ipynb'):
  nb=json.loads(data);code='\n\n'.join(''.join(c.get('source',[])) for c in nb['cells'] if c['cell_type']=='code');f.with_suffix('.code.txt').write_text(code)
 return path,len(data)
with concurrent.futures.ThreadPoolExecutor(4) as ex:
 for r in ex.map(get,paths):print(r,flush=True)
print('AUTHOR_FETCH_COMPLETE',sha)
