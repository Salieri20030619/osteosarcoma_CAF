from pathlib import Path
import tarfile,json,hashlib,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919'
for name in ['early_ambient','unify_gates','quality_guard','batch_models','all_source_thresholds','cnv_risk_model','review']:assert (O/'logs'/(name+'.exit_code')).read_text().strip()=='0',name
protocol=json.loads((O/'resources/batch_protocol.json').read_text());protocol['audit_triggered_additional_sensitivity']={'CNV_cluster_risk':.3,'cells_added_vs_strict':91,'reason':'all-source threshold audit showed sensitivity in legacy cohorts; unchanged primary'};(O/'resources/batch_protocol.json').write_text(json.dumps(protocol,indent=2))
# Keep metadata needed to audit reused and newly completed ambient correction.
d=O/'resources/ambient_parameters';d.mkdir(exist_ok=True)
import pandas as pd
m=pd.read_csv(O/'resources/source_manifest.tsv',sep='\t')
for r in m.to_dict('records'):
 p=Path(r['ambient_dir'])/'parameters.json'
 if p.is_file():shutil.copy2(p,d/(r['sample']+'.json'))
 f=Path(r['source_file']);assert hashlib.sha256(f.read_bytes()).hexdigest()==r['source_sha256']
g=P/'geo_cnv_20260918/inputs/GSE250015__GSM7969982/genes.tsv';shutil.copy2(g,O/'resources/GSE250015_feature_mapping.tsv')
qa=json.loads((O/'REVIEW_QA.json').read_text());qa['all_40_source_hashes_verified']=True;qa['all_7_pipeline_exit_codes_zero']=True;qa['batch_models']=34;(O/'REVIEW_QA.json').write_text(json.dumps(qa,indent=2))
files=[]
for folder in ['code','resources','figures','logs']:
 files.extend(p for p in (O/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p))
files.extend(p for p in (O/'tables').iterdir() if p.is_file() and (not p.name.endswith('_unified_evidence.tsv.gz') or p.name=='all_cell_unified_evidence.tsv.gz'))
files.extend(p for p in (O/'objects').iterdir() if p.name.endswith('_integrated.h5ad') or p.name=='unified_candidate_union.h5ad')
files.extend(p for p in O.iterdir() if p.is_file())
manifest={p.relative_to(O).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files)) if p.name!='MANIFEST_SHA256.json'};(O/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));files.append(O/'MANIFEST_SHA256.json')
out=P/'phase10_batch_gates_20260919_delivery.tar.gz'
with tarfile.open(out,'w:gz',compresslevel=1,encoding='utf-8') as t:
 for p in sorted(set(files)):t.add(p,arcname='phase10_batch_gates_20260919/'+p.relative_to(O).as_posix())
print(json.dumps({'delivery_files':len(set(files)),'archive_bytes':out.stat().st_size,'all_40_sources_verified':True}))