from pathlib import Path
import json,hashlib,zipfile,datetime,os
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase55_annotation_freeze_20261004';B=P/'phase56_patient_programs_20261004';O=P/'phase57_downstream_update_20261004'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8388608),b''):h.update(c)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
checks={}
for p,n in [(A,'freeze'),(B,'pseudobulk'),(B,'programs'),(O,'downstream'),(O,'report'),(O/'independent_QA','independent_R_QA')]:
 code=int((p/(n+'.exit_code')).read_text());assert code==0;checks[n]=code
for p,n in [(A,'ANNOTATION_FREEZE_COMPLETE'),(B,'PSEUDOBULK_COMPLETE'),(B,'PROGRAMS_COMPLETE'),(O,'DOWNSTREAM_COMPLETE'),(O/'independent_QA','INDEPENDENT_R_QA_COMPLETE'),(O,'REPORT_READY')]:assert (p/n).is_file()
annotation=json.loads((A/'delivery/FINAL_ACCEPTANCE.json').read_text());assert sha(A/'tables/all_cell_annotation_v2.tsv.gz')==annotation['annotation_SHA256']
current=P/'CAF_NEXT_STAGE_CURRENT.json';before=current.read_bytes();c=json.loads(before);previous=O/'resources/PREVIOUS_CURRENT.json';assert not previous.exists();previous.write_bytes(before)
objects=[{'path':str(p),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in sorted((B/'objects').iterdir()) if p.is_file()];pd.DataFrame(objects).to_csv(O/'delivery/SERVER_OBJECT_MANIFEST.tsv',sep='\t',index=False)
report=next((O/'delivery').glob('*报告*.txt'));text=report.read_text(encoding='utf-8-sig');text+='\n交付体积说明：伪bulk及逐患者全基因效应NPZ对象保留sgs10，SERVER_OBJECT_MANIFEST.tsv列出每个对象的路径、字节数和SHA256。压缩包包含全细胞注释、全部统计表、逐患者程序效应、代码、图件和验收记录，不重复打包原始H5AD或大型NPZ。\n';report.write_text(text,encoding='utf-8-sig')
(B/'PROGRESS.json').write_text(json.dumps({'state':'patient_programs_complete','completed_samples':107,'total_samples':107,'programs_exit_code':0}))
qaP=pd.read_csv(O/'independent_QA/independent_R_P_QA.tsv',sep='\t');qaF=pd.read_csv(O/'independent_QA/independent_R_FDR_QA.tsv',sep='\t')
accept={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':'requested_four_step_sequence_complete','stages':checks,'input_cells':876648,'sources':9,'samples':107,'clusters':307,'candidate_scope_counts':annotation['scope_counts'],'programs':84,'program_statistics_rows':9408,'gene_contexts_evaluable':28,'gene_contexts_unassessable':8,'statistical_QA':{'all_program_P_and_sampled_gene_P_rows':int(qaP.rows.sum()),'max_P_error':float(qaP.max_P_error.max()),'FDR_families':len(qaF),'FDR_rows':int(qaF.rows.sum()),'max_FDR_error':float(qaF.max_FDR_error.max())},'visual_QA':'All three final PNGs inspected: axes, source labels, symbols, NA vs zero, legends and clipping passed.','annotation':str(A/'tables/all_cell_annotation_v2.tsv.gz'),'report':str(report),'current_result_index':str(O/'tables/DOWNSTREAM_RESULT_INDEX.tsv'),'server_objects_manifest':str(O/'delivery/SERVER_OBJECT_MANIFEST.tsv'),'large_objects_in_zip':False,'CNV_rerun':False,'original_versions_overwritten':False,'confirmed_CAF_identity_claimed':False,'key_findings':['Core matrisome mean direction positive in seven evaluable agreed11 sources; strength depends on source and correction.','SCPCP000017 core matrisome corrected P=0.0549, despite preserved positive direction.','Added D05 candidates do not meet nominal P<0.05 for core matrisome or the three-gene program.','Existing fixed spatial statistics retained; attribution updated to RNA/program association, not CAF-specific causation.']}
dump(O/'delivery/FINAL_ACCEPTANCE.json',accept);(O/'SEQUENCE_COMPLETE').write_text('SUCCESS\n')
files=[]
for stage in [A,B,O]:
 for p in stage.rglob('*'):
  if not p.is_file() or '__pycache__' in p.parts or 'objects' in p.relative_to(stage).parts or p.suffix=='.zip' or p.name.endswith('.pid'):continue
  if p.name in ['PACKAGE_RECEIPT.json','PACKAGE_FILE_MANIFEST.tsv','PACKAGES_COMPLETE','package.exit_code']:continue
  files.append(p)
files=sorted(set(files));manifest=[{'relative_path':str(p.relative_to(P)),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in files]
m=O/'delivery/PACKAGE_FILE_MANIFEST.tsv';pd.DataFrame(manifest).to_csv(m,sep='\t',index=False)
zpath=O/'delivery/CAF_reviewed_candidate_analysis_20261004.zip'
with zipfile.ZipFile(zpath,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
 for p in files:z.write(p,str(p.relative_to(P)))
 z.write(m,'PACKAGE_FILE_MANIFEST.tsv')
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for item in manifest:
  assert z.getinfo(item['relative_path']).file_size==item['bytes']
  h=hashlib.sha256()
  with z.open(item['relative_path']) as f:
   for chunk in iter(lambda:f.read(8388608),b''):h.update(chunk)
  assert h.hexdigest()==item['SHA256']
receipt={'archive':str(zpath),'bytes':zpath.stat().st_size,'files':len(files)+1,'SHA256':sha(zpath),'CRC_and_all_entry_SHA256_verified':True,'server_objects':len(objects),'server_object_bytes':sum(x['bytes'] for x in objects),'server_object_manifest':str(O/'delivery/SERVER_OBJECT_MANIFEST.tsv')};dump(O/'delivery/PACKAGE_RECEIPT.json',receipt)
assert current.read_bytes()==before
c['at']=accept['at'];c['phase']=O.name;c['state']='phase55_57_annotation_patient_programs_downstream_complete';c['previous_current_reference']=str(previous);c['current_result_index']=str(O/'tables/DOWNSTREAM_RESULT_INDEX.tsv');c['latest_analysis_report']=str(report);c['latest_statistics_archive']=receipt;c['formal_client_delivery_baseline']='phase49 formal report/archive preserved; phase55-57 is the current reviewed-candidate analysis supplement.'
c['latest_research_followup']={'question':'Sequential reviewed annotation, patient distribution/programs, sensitivity and downstream interpretation','state':'complete','phases':[A.name,B.name,O.name],'annotation':str(A/'tables/all_cell_annotation_v2.tsv.gz'),'report':str(report),'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'packages':{'analysis':receipt},'scientific_decision':accept['key_findings'],'formal_CAF_identity_confirmed':False,'CNV_rerun':False,'scope_completed':True,'next_action':'Use the new candidate annotation and patient-based results with the risk flags and source-specific sensitivity conclusions; previous F3/F4 outputs remain historical operational-label results.'}
tmp=O/'resources/CURRENT_STATE.pending.json';dump(tmp,c);os.replace(tmp,current);(O/'package.exit_code').write_text('0');(O/'PACKAGES_COMPLETE').write_text('SUCCESS\n');print(json.dumps(receipt,ensure_ascii=True))