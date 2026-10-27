from pathlib import Path
import shutil,json,hashlib,subprocess,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';A=R/'phase25_complete_nonPDX_assets_20260921/references/HLCA'
for d in ['code','resources','logs','reference']:(O/d).mkdir(parents=True,exist_ok=True)
for name,md5 in [('HLCA_reference_model.zip','282060f9e4dc9281e02b1ba56b80cdaa'),('HLCA_reference_model_gene_order_ids_and_symbols.csv','53d367dbf4f30f766060ff92a83cd432')]:
 p=A/name;assert hashlib.md5(p.read_bytes()).hexdigest()==md5;dest=O/'reference'/name
 if dest.exists() and hashlib.md5(dest.read_bytes()).hexdigest()!=md5:dest.rename(O/'resources'/(name+'.interrupted_transfer'))
 if not dest.exists():shutil.copyfile(p,dest)
shutil.copyfile(R/'caf122_preflight.py',O/'code/01_preflight.py')
rec=dict(at=datetime.datetime.now(datetime.timezone.utc).isoformat(),reference_cache=str(A),historical_mapping=str(R/'phase32_lung_sensitivity_20260926'),correction='Historical8174-cell HLCA mapping found after server recovery; previously unavailable local records did not establish current four-candidate coverage. Reassess coverage and map NCI rather than assuming no historical work.',coverage_policy='No unsupported95percent hard gate. Use source-available canonical features for model input; report missingness explicitly, assess uncertainty and low-coverage sources. Common-gene expression comparison remains separate.')
(O/'resources/RECOVERY_AND_REFERENCE_REUSE.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
with (O/'logs/01_preflight.log').open('w',encoding='utf-8') as log:
 p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python','-u',str(O/'code/01_preflight.py')],stdout=log,stderr=subprocess.STDOUT);(O/'ACTIVE_PROCESS.json').write_text(json.dumps(dict(stage='raw_query_preflight',pid=p.pid)),encoding='utf-8');rc=p.wait()
(O/'resources/PREFLIGHT_RUN_RECEIPT.json').write_text(json.dumps(dict(exit_code=rc,marker=(O/'QUERY_PREFLIGHT.json').exists())),encoding='utf-8');print((O/'logs/01_preflight.log').read_text()[-10000:]);raise SystemExit(rc)
