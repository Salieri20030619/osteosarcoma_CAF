from pathlib import Path
import json,csv,zipfile,hashlib,datetime,os
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');d=r/'review_intake_20261004'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(8388608),b''):h.update(chunk)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def rows(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
proof=read(d/'REVIEW_INTAKE_LOCAL_VALIDATION.json');archive=d/'CAF_blind_review_comparison_20261003.zip'
assert sha(archive)==proof['archive_SHA256']
current=r/'CAF_NEXT_STAGE_CURRENT.json';before=current.read_bytes();c=json.loads(before)
assert c['phase']=='phase53_de_novo_reannotation_20261002'
assert proof['freeze']['input_zip_SHA256']==c['latest_research_followup']['packages']['blind']['SHA256']
b=d/'received_review'
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for f in z.infolist():
  p=(b/f.filename).resolve();assert p.is_relative_to(b.resolve())
  if f.is_dir():p.mkdir(parents=True,exist_ok=True);continue
  content=z.read(f)
  if p.exists():assert p.read_bytes()==content
  else:
   p.parent.mkdir(parents=True,exist_ok=True)
   with p.open('xb') as stream:stream.write(content)
freeze=read(b/'blind_submission/BLIND_FREEZE.json')
assert sha(b/'blind_submission/BLIND_FREEZE_MANIFEST.tsv')==freeze['freeze_manifest_SHA256']
fr=rows(b/'blind_submission/BLIND_FREEZE_MANIFEST.tsv')
for x in fr:assert sha(b/'blind_submission'/x['path'])==x['SHA256']
table=rows(b/'unblinded_comparison/cluster_annotation_comparison.tsv')
assert len(table)==307 and sum(int(x['cells']) for x in table)==876648
mixed={'D01_R10_K002','D01_R10_K034','D04_R10_K007','D06_R10_K022','D07_R10_K027','D09_R10_K010'}
expansion={'D02_R10_K012','D05_R10_K018','D07_R10_K024','D07_R10_K041'}
fields=['source_id','dataset','cluster_id','original_cluster','cells','proposed_lineage','reviewed_lineage','malignancy_assessment','recommended_action']
followup={'state':'pending_scientific_followup_not_started','priority_1_mixed_expression':[ {k:x[k] for k in fields} for x in table if x['cluster_id'] in mixed],'priority_1_fibroblast_CNV_conflict':[ {k:x[k] for k in fields} for x in table if x['cluster_id']=='D01_R10_K017'],'priority_2_fibroblast_expansion':[ {k:x[k] for k in fields} for x in table if x['cluster_id'] in expansion],'other_open_items':['Other fibroblast/mesenchymal CNV conflicts and 14 mesenchymal boundary clusters','D08 common-feature/full-library difference: 122334 counts; raw H5AD attribution not yet independently reconciled','Alternative CNV reference sensitivity and state-dominated subclusters'],'preserve':['Frozen blind submission','Original provisional annotation','phase49 formal statistical and report baseline'],'formal_CAF_identity_confirmed':False,'review_CNV_heuristics_not_adopted_as_cell_filters':True}
assert len(followup['priority_1_mixed_expression'])==6 and len(followup['priority_2_fibroblast_expansion'])==4
dump(d/'FOLLOWUP_REGISTER.json',followup)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
receipt={'at':now,'state':'blind_review_and_unblinded_comparison_received_verified','archive':str(archive),'SHA256':proof['archive_SHA256'],'archive_files':66,'frozen_files_verified':len(fr),'clusters':307,'cells':876648,'comparison_report':str(b/'unblinded_comparison'/('\u63ed\u76f2\u6bd4\u5bf9\u62a5\u544a.txt')),'followup_register':str(d/'FOLLOWUP_REGISTER.json'),'local_verification':str(d/'REVIEW_INTAKE_LOCAL_VALIDATION.json'),'review_provenance_scope':proof['review_provenance_scope'],'formal_annotations_replaced':False,'downstream_statistics_recomputed':False,'next_action':'Resolve mixed-expression and fibroblast/CNV identity conflicts from raw sample-level evidence before adjudicating expansion candidates.'}
dump(d/'REVIEW_INTAKE_ACCEPTANCE.json',receipt)
backup=d/'CAF_NEXT_STAGE_CURRENT.before_review_intake.json'
if backup.exists():assert backup.read_bytes()==before
else:
 with backup.open('xb') as f:f.write(before)
c['at']=now;c['state']='phase53_blind_review_compared_followup_pending'
f=c['latest_research_followup'];f['independent_blind_review_performed']=True
f['review_evidence_basis']='User-supplied frozen blind review and separate unblinded comparison; receipt includes provenance scope.'
f['next_action']=receipt['next_action'];f['blind_review_intake']=receipt
c['latest_review_intake']=str(d/'REVIEW_INTAKE_ACCEPTANCE.json')
assert current.read_bytes()==before
staging=d/'CURRENT_STATE.pending.json';dump(staging,c);os.replace(staging,current)
assert read(current)['latest_review_intake']==str(d/'REVIEW_INTAKE_ACCEPTANCE.json')
(d/'REVIEW_INTAKE_COMPLETE').write_text('SUCCESS\n')
print(json.dumps(receipt,ensure_ascii=True))