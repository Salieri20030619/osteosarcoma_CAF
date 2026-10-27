from pathlib import Path
import json,hashlib,zipfile,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase62_F1F5_internal_20261005')
assert (O/'VERIFIED').exists() and (O/'REPORT_COMPLETE').exists()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
files=sorted(p for folder in ['code','resources','tables','figures','QA','logs','objects'] for p in (O/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p));files.append(O/'delivery/F1-F5内部亚群与置信度报告_20261005.txt');files.extend(O/f for f in ['prepare.exit_code','gsea.exit_code','ANALYSIS_COMPLETE','PREPARE_COMPLETE','GSEA_COMPLETE','REPORT_COMPLETE','VERIFIED'])
manifest=[{'file':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files];mp=O/'delivery/FILE_MANIFEST.json';mp.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');files.append(mp)
archive=O/'delivery/CAF_F1F5_internal_confidence_20261005.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in manifest:assert hashlib.sha256(z.read(r['file'])).hexdigest()==r['sha256']
receipt={'status':'PASS','at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive':str(archive),'bytes':archive.stat().st_size,'sha256':sha(archive),'files':len(files),'all_entry_hashes_verified':True,'report_sha256':sha(O/'delivery/F1-F5内部亚群与置信度报告_20261005.txt'),'summary_sha256':sha(O/'tables/F1F5_subgroup_summary.tsv')};(O/'delivery/PACKAGE_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt,indent=2),flush=True)
