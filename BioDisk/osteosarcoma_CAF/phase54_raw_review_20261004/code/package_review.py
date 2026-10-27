from pathlib import Path
import json,hashlib,zipfile,datetime,os
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase54_raw_review_20261004'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
qa=json.loads((O/'FINAL_DATA_QA.json').read_text());qa['figure_visual_QA']='Both final PNGs inspected: labels, legends, numeric display and clipping passed; raw and fractional >0 values shown alongside >=1 sensitivity.';qa['finalize_exit_code']=int((O/'finalize.exit_code').read_text());assert qa['finalize_exit_code']==0;dump(O/'FINAL_DATA_QA.json',qa)
current=P/'CAF_NEXT_STAGE_CURRENT.json';before=current.read_bytes();c=json.loads(before);backup=O/'resources/PREVIOUS_CURRENT.json';assert not backup.exists();backup.write_bytes(before)
report=next((O/'delivery').glob('*报告*.txt'))
accept={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':'review_complete_scientific_identity_limits_retained','samples':107,'input_cells':876648,'reviewed_clusters':42,'reviewed_cells':127395,'raw_counts_and_correction_review_complete':True,'CNV_review_complete':True,'D08_count_discrepancy_resolved':True,'mixed_primary_risks':{'supported':5,'sparse_signal_only':1},'four_expansion_candidates':'Retain differentiated provisional evidence; no automatic confirmed CAF expansion.','confirmed_CAF_identity':False,'formal_labels_changed':False,'formal_statistics_recomputed':False,'report':str(report),'cluster_decisions':str(O/'tables/CLUSTER_REVIEW_DECISIONS.tsv'),'cell_overlay':str(O/'tables/cell_review_decisions.tsv.gz'),'full_307_cluster_identity_review_claimed':False,'remaining':['Malignant/nonmalignant and clone identity require independent evidence; RNA CNV alone does not establish it.','Sparse/low-depth, mixed and mesenchymal boundary identities remain unresolved.','Other state-led and osteochondrogenic subtype disagreements were retained outside the 42-cluster raw deep-review scope.']}
dump(O/'delivery/FINAL_ACCEPTANCE.json',accept);(O/'REVIEW_COMPLETE').write_text('SUCCESS\n')
files=[]
for sub in ['tables','cells','figures','code','resources','logs','delivery']:
 files.extend(p for p in (O/sub).rglob('*') if p.is_file() and p.suffix!='.zip' and '__pycache__' not in p.parts)
files.extend(p for p in O.iterdir() if p.is_file() and (p.suffix=='.json' or p.name.endswith('.exit_code') or p.name.endswith('_COMPLETE')))
files=sorted(set(files));manifest=[{'relative_path':str(p.relative_to(O)),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in files]
m=O/'delivery/PACKAGE_FILE_MANIFEST.tsv';m.write_text('relative_path\tbytes\tSHA256\n'+'\n'.join(f"{x['relative_path']}\t{x['bytes']}\t{x['SHA256']}" for x in manifest)+'\n')
zpath=O/'delivery/CAF_raw_review_20261004.zip'
with zipfile.ZipFile(zpath,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
 z.write(m,'PACKAGE_FILE_MANIFEST.tsv')
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for e in manifest:
  assert z.getinfo(e['relative_path']).file_size==e['bytes']
  with z.open(e['relative_path']) as f:
   h=hashlib.sha256()
   for b in iter(lambda:f.read(8388608),b''):h.update(b)
  assert h.hexdigest()==e['SHA256']
receipt={'archive':str(zpath),'bytes':zpath.stat().st_size,'files':len(files)+1,'SHA256':sha(zpath),'CRC_and_all_entry_hashes_verified':True};dump(O/'delivery/PACKAGE_RECEIPT.json',receipt)
assert current.read_bytes()==before
c['at']=accept['at'];c['state']='phase54_raw_review_complete_identity_limits_retained';c['phase']=O.name;c['previous_current_reference']=str(backup)
c['latest_research_followup']={'question':'Post-unblinding raw count, ambient/doublet and RNA-CNV review','acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'report':str(report),'packages':{'review':receipt},'reviewed_clusters':42,'reviewed_cells':127395,'input_samples':107,'input_cells':876648,'scientific_decision':'Mixed risks stratified; D08 count discrepancy closed; candidate fibroblast expression retained separately from unconfirmed CAF identity.','independent_blind_review_performed':True,'this_stage_is_post_unblinding_followup':True,'formal_labels_replaced':False,'next_action':'Use the phase54 risk/decision overlay for interpretation. Do not promote all candidate cells to confirmed CAF; unresolved identity questions remain explicit.','scope_completed':True,'remaining':accept['remaining']}
pending=O/'resources/CURRENT_STATE.pending.json';dump(pending,c);os.replace(pending,current)
(O/'package.exit_code').write_text('0');(O/'PACKAGES_COMPLETE').write_text('SUCCESS\n');print(json.dumps(receipt,ensure_ascii=True))