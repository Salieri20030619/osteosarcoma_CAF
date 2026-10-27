from pathlib import Path
import json,hashlib,pandas as pd,zipfile,re,datetime
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'review_intake_20261005_checklist'
r=pd.read_csv(O/'tables/checklist_22_item_review.tsv',sep='\t',encoding='utf-8-sig',dtype=str);assert len(r)==22 and r.iloc[:,0].is_unique
src=(O/'sources/extracted_text.txt').read_text();ids=re.findall(r'^\[\d+\]\s+(\d\.\d)\s+',src,re.M);assert set(ids)==set(r.iloc[:,0]) and len(ids)==22
rec=json.loads((O/'sources/source_receipt.json').read_text());assert hashlib.sha256((O/'sources/client_checklist_original.docx').read_bytes()).hexdigest()==rec['sha256']
assert hashlib.sha256((O/'delivery/review_report.txt').read_bytes()).hexdigest()=='91dae80636cb1bf811ae7747d32d6d294c1e4a504daefc103b6470d0fffa1461'
q=json.loads((O/'resources/AUDIT_ACCEPTANCE.json').read_text());assert q['passed'] and q['upstream_assessment_coverage']==8174 and q['additional_positive_calls']==15
assert json.loads((B/'CAF_NEXT_STAGE_CURRENT.json').read_text())['current_phase']=='phase65_spatial_structure_20261005'
actions=[
 {'priority':1,'action':'Restore all8174 doublet evidence fields from source with method/threshold provenance','status':'recommended_not_implemented','reason':'5728 downstream records missing;15 additional model positives'},
 {'priority':1,'action':'Prespecified retain/exclude26model-positive-cell sensitivity, including F2_C1','status':'recommended_not_run','reason':'No direction or P-value impact has yet been assessed; F2_C3 source calls0/38'},
 {'priority':2,'action':'Reconcile ambient fraction metadata with full corrected matrix provenance','status':'recommended_not_implemented','reason':'2446legacy fraction entries versus8174actual corrected count records'},
 {'priority':2,'action':'High-resolution pathology/independent contours and collagen deposition measurements','status':'requires_materials','reason':'26lowres images,0accepted contours; transcriptomic association does not establish barrier'}
]
(O/'resources/FOLLOWUP_REGISTER.json').write_text(json.dumps(actions,indent=2))
(O/'REVIEW_COMPLETE.json').write_text(json.dumps({'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'items':22,'content_audit_complete':True,'read_only_provenance_audit_complete':True,'data_or_models_changed':False,'active_phase_unchanged':True,'word_layout_review_complete':False},indent=2))
files=[p for p in O.rglob('*') if p.is_file() and p.suffix!='.zip' and '__pycache__' not in p.parts]
manifest=[{'path':str(p.relative_to(O)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(files)]
(O/'resources/ARCHIVE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
p=O/'delivery/CAF_client_checklist_audit_20261005.zip'
with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files+[O/'resources/ARCHIVE_MANIFEST.json']:z.write(f,arcname=O.name+'/'+str(f.relative_to(O)))
with zipfile.ZipFile(p) as z:assert z.testzip() is None;n=len(z.infolist())
receipt={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'files':n}
(O/'delivery/ARCHIVE_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
