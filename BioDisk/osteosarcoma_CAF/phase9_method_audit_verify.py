from pathlib import Path
import tarfile,json,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase9_method_audit_20260919';J=P/'phase9_joint_20260919'
with tarfile.open(P/'phase9_method_audit_20260919_final.tar.gz',encoding='utf-8') as t:
 for m in t.getmembers():assert (O/m.name).resolve().is_relative_to(O.resolve()) and not m.issym() and not m.islnk()
 t.extractall(O)
manifest=json.loads((O/'MANIFEST_SHA256.json').read_text())
for name,want in manifest.items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==want,name
original=json.loads((J/'MANIFEST_SHA256.json').read_text())
for name in ['objects/joint_CAF_integrated.h5ad','PROTOCOL.json','STATUS.json']:assert hashlib.sha256((J/name).read_bytes()).hexdigest()==original[name],name
print(json.dumps({'audit_files_verified':len(manifest),'production_object_protocol_status_unchanged':True,'main_project_state_not_modified':True,'audit_status':json.loads((O/'STATUS.json').read_text())['status']}))