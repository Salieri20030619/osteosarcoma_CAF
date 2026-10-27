from pathlib import Path
import json,hashlib,zipfile,shutil,datetime
ROOT=Path(__file__).resolve().parent;D=ROOT.parent/'CAF_新增数据入口复核_20261013'
for name in ['NCC_DOWNLOAD_COMPLETE.json','NCC_MATRIX_ACCEPTANCE.json','NCC_BARCODE_COMPARISON.json','RECHECK_FINDINGS.json']:
 shutil.copy2(ROOT/name,D/name)
for name in ['NCC_MATRIX_AND_PATIENT_INVENTORY.tsv','NCC_BEST_ORIGINAL_BARCODE_OVERLAP.tsv','ADDITIONAL_HRA_SAMPLE_RECORDS.tsv','OFFICIAL_HRA_OSTEOSARCOMA_22_RECORDS.tsv']:
 shutil.copy2(ROOT/'tables'/name,D/name)
files=[]
for p in ROOT.rglob('*'):
 if not p.is_file() or 'downloads' in p.relative_to(ROOT).parts or '__pycache__' in p.parts:continue
 if p.suffix=='.png':continue
 files.append(p)
report=D/'CAF_NCC_HRA及相关文献复核报告_20261013.txt'
assert report.exists() and (ROOT/'NCC_BARCODE_COMPARISON.json').exists()
manifest=[];zpath=D/'CAF_NCC_HRA_source_access_recheck_20261013.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(files):
  b=p.read_bytes();name=str(p.relative_to(ROOT)).replace('\\','/');z.writestr(name,b);manifest.append({'path':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
 z.write(report,report.name)
 z.writestr('BUNDLE_MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
receipt={'package':zpath.name,'bytes':zpath.stat().st_size,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'source_files':len(files),'contains_raw_NCC_matrices':False,'raw_data_local_directory':str(ROOT/'downloads/NCC'),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(D/'PACKAGE_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf8');print(json.dumps(receipt,ensure_ascii=False))
