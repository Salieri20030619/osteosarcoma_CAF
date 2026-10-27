from pathlib import Path
import json,hashlib,zipfile,csv,io,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');R=P/'phase49_audit_revision_20261002';D=R/'delivery/client';D.mkdir(exist_ok=True)
old=P/'phase48_clinical_programs_delivery_20261001/delivery/client/CAF全转录组与空间基质研究_完整统计包_20261002.zip'
archive=D/'CAF全转录组与空间基质研究_修订统计包_20261002_v2.zip'
assert (R/'COMPUTATION_VALIDATED').read_text().strip()=='SUCCESS'
assert (R/'surface/SURFACE_INDEPENDENT_R_QA_COMPLETE').read_text().strip()=='SUCCESS'
qa=json.loads((D/'REPORT_VISUAL_QA.json').read_text());assert qa['all_pages_visual_passed']
entries={};historical={'CAF客户结果报告.pdf','FINAL_ACCEPTANCE.json','README_阅读说明.txt','CONDITIONAL_ENDPOINTS.tsv'}
# Reuse verified original content, explicitly preserving old root delivery as history.
with zipfile.ZipFile(old) as zin:
 for info in zin.infolist():
  if info.filename=='FILE_MANIFEST.tsv':continue
  name='history_v1/'+info.filename if info.filename in historical else info.filename
  entries[name]=('old',info.filename)
# Include every file in previously omitted, documented clinical revision subdirectories.
required=[]
for sub in ['clinical_site_review','clinical_reporting_extension','clinical_model_sensitivity']:
 for f in (P/'phase44_patient_mechanism_20261001'/sub).rglob('*'):
  if f.is_file() and '__pycache__' not in f.parts:
   name=f.relative_to(P).as_posix();entries[name]=('file',str(f));required.append(name)
# All revised non-object deliverables, independent of first-level subdirectory names.
for f in R.rglob('*'):
 if not f.is_file():continue
 rel=f.relative_to(R)
 if any(x in rel.parts for x in ['objects','client','pdf_render','__pycache__','evaluation','internal']):continue
 if f.suffix in ['.zip','.part'] or f.name.startswith('evaluation_'):continue
 entries[f.relative_to(P).as_posix()]=('file',str(f))
for f in D.iterdir():
 if f.is_file() and f.suffix not in ['.zip'] and f.name not in ['FILE_MANIFEST.tsv','PACKAGE_INTEGRITY.json']:
  entries[f.name]=('file',str(f))
assert len(required)==42 and all(n in entries for n in required)
with (D/'CURRENT_RESULT_INDEX.tsv').open() as source:
 for row in csv.DictReader(source,delimiter='	'):
  path=row['current_path'];assert any(n.startswith(path) for n in entries) if path.endswith('/') else path in entries,path
records=[]
with zipfile.ZipFile(old) as zin,zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as zout:
 for name,(kind,loc) in sorted(entries.items()):
  data=zin.read(loc) if kind=='old' else Path(loc).read_bytes();zout.writestr(name,data);records.append({'path':name,'bytes':len(data),'SHA256':hashlib.sha256(data).hexdigest()})
 out=io.StringIO();w=csv.DictWriter(out,fieldnames=['path','bytes','SHA256'],delimiter='	',lineterminator='\n');w.writeheader();w.writerows(records);zout.writestr('FILE_MANIFEST.tsv',out.getvalue());(D/'FILE_MANIFEST.tsv').write_text(out.getvalue())
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in records:
  data=z.read(r['path']);assert len(data)==r['bytes'] and hashlib.sha256(data).hexdigest()==r['SHA256']
 assert all(n in z.namelist() for n in required)
 with zipfile.ZipFile(old) as zin:
  for n in zin.namelist():
   if n=='FILE_MANIFEST.tsv':continue
   dst='history_v1/'+n if n in historical else n;assert zin.read(n)==z.read(dst)
check={'at':datetime.datetime.now().astimezone().isoformat(),'state':'revised_package_accepted','archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_SHA256':hashlib.sha256(archive.read_bytes()).hexdigest(),'files':len(records),'CRC_passed':True,'all_entries_SHA_and_size_verified':True,'missing_clinical_files_added':42,'all_original_entries_preserved_bytewise':True,'current_result_index_included':'CURRENT_RESULT_INDEX.tsv' in entries};(D/'PACKAGE_INTEGRITY.json').write_text(json.dumps(check,indent=2,ensure_ascii=False));(R/'REVISED_PACKAGE_COMPLETE').write_text('SUCCESS');print(json.dumps(check,ensure_ascii=False),flush=True)
