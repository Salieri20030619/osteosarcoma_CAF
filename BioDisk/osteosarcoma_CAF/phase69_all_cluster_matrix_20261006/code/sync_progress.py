
from pathlib import Path
import json,hashlib,zipfile,datetime,shutil
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase69_all_cluster_matrix_20261006'
receipt=json.loads((O/'delivery/DELIVERY_RECEIPT.json').read_text())
z=O/'delivery/results.zip';assert hashlib.sha256(z.read_bytes()).hexdigest()==receipt['zip_sha256']
with zipfile.ZipFile(z) as a:
 assert a.testzip() is None
 manifest=json.loads(a.read('FILE_MANIFEST.json'))
 for r in manifest:assert hashlib.sha256(a.read(r['path'])).hexdigest()==r['sha256']
 for r in manifest:
  rel=r['path'];dest=None
  if rel.startswith('analysis/tables/'):dest=O/'tables'/Path(rel).name
  elif rel.startswith('analysis/resources/'):dest=O/'resources'/Path(rel).name
  elif rel.startswith('delivery/'):dest=O/'delivery'/Path(rel).name
  elif rel.startswith('figures/'):dest=O/rel
  elif rel in ['analyze.py','plot_figures.py','build_report.py','report_helpers.py','report_style_reference.docx']:dest=O/'code'/rel
  elif rel in ['FINAL_ACCEPTANCE.json','README.txt']:dest=O/'delivery'/rel
  if dest:
   dest.parent.mkdir(exist_ok=True,parents=True);dest.write_bytes(a.read(rel))
qa=json.loads((O/'delivery/FINAL_ACCEPTANCE.json').read_text())
assert qa['passed'] and json.loads((O/'QA/EXECUTION_STATUS.json').read_text())['exit_code']==0
assert (O/'ANALYSIS_COMPLETE').read_text().strip()=='SUCCESS'
report=O/'delivery'/Path(receipt['report'].replace('\\','/')).name
assert hashlib.sha256(report.read_bytes()).hexdigest()==receipt['report_sha256']
# Keep a human-readable, dated archive name as well as the ASCII transfer name.
named_zip=O/'delivery'/Path(receipt['zip'].replace('\\','/')).name
shutil.copy2(z,named_zip)
at=datetime.datetime.now(datetime.timezone.utc).isoformat()
record={'title':'全部33簇和所有来源的胶原与基质表达核查','phase':O.name,'status':'completed_verified_delivered','completed_at':at,'report':str(report),'local_report':receipt['report'],'results_zip':str(named_zip),'local_results_zip':receipt['zip'],'pages':6,'figures':4,'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'remaining_required_work':[],'scope':'33簇、116来源-簇组合；继承完整P/FDR，54条名义上调记录，新纳入9条非主要来源信号；方向敏感性和谱系核查；无新增推断性检验。'}
summary={'subgroups':33,'source_subgroup_pairs':116,'baseline_gene_rows_reproduced':151340,'nominal_rows':54,'nominal_subgroups':7,'new_nonprimary_rows':9,'new_subgroups':['F1_originalC4','F4_C1','F4_C9','F5_originalC12'],'FDR_significant_rows':0,'new_priority_candidates':['F4_C9 COL6A3','F1_originalC4 ADAMTS14'],'sensitivity_findings':['F4_C9 COL6A3 Fibro N7 positive; min5 N1','F1_originalC4 ADAMTS14 Fibro N8 positive; min5 N4 positive','F4_C1 SDC2 Fibro N3 negative and min5 N2 negative','F5_originalC12 S18 has no Fibro cells; cannot attribute signal directly to CAF'],'low_replication_candidates':['F4_C5','F4_C8','F1_originalC10','F3_C0','F2_C5'],'nonfibro_descriptive':['F4_C12'],'important_boundary':'Upregulation is relative to other same-patient candidate cells; no upregulation does not imply no expression/contribution; RNA not actual matrix deposition.'}
next_action='本轮全簇核查已完成。后续胶原与基质验证并列关注F2-C3、F2-C1、F4-C9/COL6A3、F1-C4/ADAMTS14；F4-C5/F4-C8等保留低重复数候选。继续使用Phase66元数据入口；耐药及抑制机制暂缓。'
for name in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
 path=B/name;d=json.loads(path.read_text())
 if d.get('last_completed_task',{}).get('phase')!=O.name:d['previous_completed_task']=d.get('last_completed_task')
 if 'collagen_gene_stage_summary' in d:
  if d['collagen_gene_stage_summary'].get('phase')!=O.name:d['phase68_primary_source_gene_summary']=d['collagen_gene_stage_summary']
 d.update({'current_phase':O.name,'status':'all_cluster_matrix_verified_delivered','task_status':'completed_verified_delivered','active_task':None,'last_completed_task':record,'last_completed_analysis':record,'latest_task':record,'current_state':'all_cluster_matrix_verified_delivered','updated_at':at,'last_activity_at':at,'progress_updated_at':at,'report':str(report),'results_zip':str(named_zip),'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'local_delivery_verified':True,'next_action':next_action,'current_collagen_gene_evidence':str(O/'tables/nominal_upregulated_all_sources_with_sensitivity.tsv'),'candidate_priorities':str(O/'tables/all_33_subgroup_evidence_and_priority.tsv'),'primary_candidates':['F2_C3 COL3A1 COL5A1 COL5A2 PCOLCE LOXL2','F2_C1 COL4A1 COL4A2 PXDN LOXL2','F4_C9 COL6A3','F1_originalC4 ADAMTS14'],'secondary_candidates':['F4_C5','F4_C8','F1_originalC10','F3_C0','F2_C5'],'conditional_candidates':['F2_C0 gene-level only','F4_C1 SDC2 lineage sensitive','F5_originalC12 source/lineage-specific','F4_C12 predominantly Osteo'],'collagen_gene_stage_summary':dict(summary,phase=O.name),'all_cluster_matrix_summary':summary,'analysis_now_completed':True,'acceptance_passed':True,'execution_status':json.loads((O/'QA/EXECUTION_STATUS.json').read_text()),'process_observation':'Native SSH analysis exited 0; recorded exit status, completion marker, output integrity and final report visual QA verified.','followup_recommendations':{'text':next_action,'status':'future_research_recommendations_not_active_task'},'deferred_objectives':['耐药机制验证','抑制相关机制验证'],'resume_record':str(O/'PHASE69_DELIVERED.json')})
 if 'phase' in d:d['phase']=O.name
 tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n');tmp.replace(path)
delivered={'at':at,'record':record,'receipt':receipt,'summary':summary,'all_three_progress_entries_updated':True,'metadata_entry_unchanged':'CAF_CURRENT_CELL_EVIDENCE.json -> phase66'}
(O/'PHASE69_DELIVERED.json').write_text(json.dumps(delivered,ensure_ascii=False,indent=2)+'\n')
for name in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
 d=json.loads((B/name).read_text());assert d['current_phase']==O.name and d['active_task'] is None and d['last_completed_task']['pages']==6
print(json.dumps({'passed':True,'phase':O.name,'files_verified':len(manifest),'progress_entries':3,'report':str(report)},ensure_ascii=False))
