from pathlib import Path
import json,zipfile,csv,io,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase49_audit_revision_20261002');D=R/'delivery/client';p=json.loads((D/'PACKAGE_INTEGRITY.json').read_text());a=Path(p['archive'])
h=hashlib.sha256()
with a.open('rb') as f:
 for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
assert h.hexdigest()==p['archive_SHA256'] and a.stat().st_size==p['archive_bytes']
with zipfile.ZipFile(a) as z:
 entries=list(csv.DictReader(io.StringIO(z.read('FILE_MANIFEST.tsv').decode('utf-8')),delimiter='\t'));assert len(entries)==p['files']
 assert len(set(x['path'] for x in entries))==len(entries)
 for r in entries:
  b=z.read(r['path']);assert len(b)==int(r['bytes']) and hashlib.sha256(b).hexdigest()==r['SHA256']
 names=set(z.namelist());expected=[]
 for sub in ['clinical_site_review','clinical_reporting_extension','clinical_model_sensitivity']:
  root=R.parent/'phase44_patient_mechanism_20261001'/sub
  for f in root.rglob('*'):
   if f.is_file() and '__pycache__' not in f.parts:expected.append(f.relative_to(R.parent).as_posix())
 assert len(expected)==42 and all(n in names for n in expected)
 for row in csv.DictReader(io.StringIO(z.read('CURRENT_RESULT_INDEX.tsv').decode('utf-8')),delimiter='\t'):
  name=row['current_path'];assert any(n.startswith(name) for n in names) if name.endswith('/') else name in names,name
out={'archive_SHA256':p['archive_SHA256'],'files':len(entries),'all_file_hashes_verified':True,'all_42_required_files_present':True,'all_current_index_targets_present':True,'builder_exit_status':'not captured by initial detached launcher; this independent verifier exit captured','acceptance_basis':'independent archive verification plus completed builder marker and original build receipt'}
(R/'evaluation/DELIVERED_ARCHIVE_QA.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
