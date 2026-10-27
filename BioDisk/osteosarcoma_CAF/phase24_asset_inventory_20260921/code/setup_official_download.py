from pathlib import Path
import json,os,requests,pandas as pd,py_compile
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');old=P/'phase15_scpca_download_20260920';O=P/'phase25_complete_nonPDX_assets_20260921';O.mkdir(exist_ok=True)
meta=pd.read_csv(P/'phase24_asset_inventory_20260921/tables/nonPDX_missing_matrices_55.tsv',sep='\t',keep_default_na=False)
assert len(meta)==55 and meta.scpca_sample_id.is_unique
assert meta.is_xenograft.astype(str).str.lower().eq('false').all() and meta.is_cell_line.astype(str).str.lower().eq('false').all()
TARGETS={p:g.scpca_sample_id.tolist() for p,g in meta.groupby('scpca_project_id')}
s=requests.Session();s.trust_env=False
r=s.get('https://api.scpca.alexslemonade.org/v1/ccdl-datasets/?limit=200',timeout=30);r.raise_for_status();datasets=r.json()['results'];safe=[]
for project in TARGETS:
 matches=[d for d in datasets if d.get('ccdl_project_id')==project and d.get('format')=='ANN_DATA' and not d.get('ccdl_is_merged') and d.get('ccdl_modality')=='SINGLE_CELL'];assert len(matches)==1
 ident=matches[0]['id'];r=s.get('https://api.scpca.alexslemonade.org/v1/ccdl-datasets/'+ident+'/',headers={'API-KEY':(old/'.api_token').read_text().strip()},timeout=30)
 if r.status_code!=200:raise RuntimeError('Official dataset API status '+str(r.status_code))
 d=r.json();cf=d.get('computed_file') or {};url=d.get('download_url') or cf.get('download_url');assert url
 cfg={'dataset_id':ident,'size':int(cf['size_in_bytes']),'url':url};p=O/(project+'_signed.json');p.write_text(json.dumps(cfg));p.chmod(0o600)
 safe.append({'project':project,'dataset_id':ident,'archive_bytes':cfg['size'],'requested_nonPDX_samples':len(TARGETS[project])})
code=(old/'range_zip.py').read_text().replace(str(old),str(O));start=code.index('TARGETS=');end=code.index('class RemoteZip',start);code=code[:start]+'TARGETS='+repr(TARGETS)+'\n'+code[end:]
code=code.replace("os.environ.get('SCPCA_DOWNLOAD_PROXY','http://127.0.0.1:7897')","os.environ.get('SCPCA_DOWNLOAD_PROXY','')")
code=code.replace("self.session.proxies={'http':proxy,'https':proxy}","self.session.proxies={'http':proxy,'https':proxy} if proxy else {}")
(O/'range_zip.py').write_text(code)
code=(old/'download_parallel.py').read_text().replace("os.environ.get('SCPCA_DOWNLOAD_PROXY','http://127.0.0.1:7897')","os.environ.get('SCPCA_DOWNLOAD_PROXY','')")
code=code.replace("tls.session.proxies={'http':proxy,'https':proxy}","tls.session.proxies={'http':proxy,'https':proxy} if proxy else {}")
code=code.replace('samples=9,','samples=sum(map(len,TARGETS.values())),')
needle="rec2={'project':project,'sample':sample"
insert="""assert set(a.obs['is_xenograft'].astype(str).str.lower())=={'false'}, 'PDX/source unresolved: '+sample
                assert set(a.obs['is_cell_line'].astype(str).str.lower())=={'false'}, 'Cell line/source unresolved: '+sample
                assert set(a.obs['sample_id'].astype(str))=={sample}, 'Sample identity mismatch'
                rec2={'project':project,'sample':sample"""
assert needle in code;code=code.replace(needle,insert)
(O/'download_parallel.py').write_text(code)
for f in ['range_zip.py','download_parallel.py']:py_compile.compile(str(O/f),doraise=True)
(O/'official_download_plan.json').write_text(json.dumps({'targets':TARGETS,'datasets':safe,'nonPDX_required':True,'files_per_sample':['filtered_rna.h5ad','processed_rna.h5ad','unfiltered_rna.h5ad'],'project_extra_QC_thresholds':False,'purpose':'asset acquisition only; no CAF statistical analysis','integrity':'HTTP Range length and content-range, ZIP CRC32, SHA256, object source metadata'},indent=2))
print(json.dumps(safe,indent=2));print('SETUP_COMPLETE',O)
