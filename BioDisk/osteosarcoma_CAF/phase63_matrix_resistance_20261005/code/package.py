from pathlib import Path
import json,hashlib,zipfile,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase63_matrix_resistance_20261005')
assert (O/'VERIFIED').exists() and (O/'INTEGRATION_COMPLETE').exists()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(8388608),b''):h.update(z)
 return h.hexdigest()
files=sorted(p for folder in ['code','resources','tables','objects','figures','QA','logs'] for p in (O/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p));files += sorted((O/'delivery').glob('*.txt'));files.append(O/'delivery/REPORT_FACTS.json')
files += sorted(O.glob('*.exit_code'));files += [p for p in O.iterdir() if p.is_file() and (p.name.endswith('_COMPLETE') or p.name=='VERIFIED')]
manifest=[{'file':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files];m=O/'delivery/FILE_MANIFEST.json';m.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');files.append(m)
archive=O/'delivery/CAF_matrix_resistance_round1_20261005.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in manifest:assert hashlib.sha256(z.read(r['file'])).hexdigest()==r['sha256']
receipt={'status':'PASS','archive':str(archive),'sha256':sha(archive),'bytes':archive.stat().st_size,'files':len(files),'all_entry_hashes_verified':True,'plan_sha256':sha(O/'delivery/CAF基质耐药目标导向执行计划_20261005.txt'),'report_sha256':sha(O/'delivery/CAF基质耐药目标导向首轮结果_20261005.txt'),'candidates_sha256':sha(O/'tables/goal_oriented_candidate_priorities.tsv')}
(O/'delivery/PACKAGE_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2),flush=True)
