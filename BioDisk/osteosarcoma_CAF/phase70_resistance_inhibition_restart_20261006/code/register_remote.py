from pathlib import Path
import json,hashlib,datetime,shutil,socket
B=Path(__file__).resolve().parents[1];root=B.parent
manifest=json.loads((B/'resources/TRANSFER_MANIFEST.json').read_text(encoding='utf-8'))
for rel,sha in manifest.items():assert hashlib.sha256((B/rel).read_bytes()).hexdigest()==sha,rel
accept=json.loads((B/'delivery/FINAL_ACCEPTANCE.json').read_text(encoding='utf-8'))
assert accept['status']=='passed'
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
states=['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','PROJECT_STATE.json']
for name in states:
 assert hashlib.sha256((root/name).read_bytes()).digest()==hashlib.sha256((B/'resources'/name).read_bytes()).digest(), 'Concurrent state change: '+name
prior=json.loads((root/states[0]).read_text(encoding='utf-8'))
report=Path(prior['latest_report_revision']['report'])
assert hashlib.sha256(report.read_bytes()).hexdigest()==accept['original_client_report_sha256']
next_action='已恢复耐药与抑制研究并完成Phase70首轮临床基因模型。下一步沿Phase66数据入口核查CAF发送端、肿瘤接收端ABCG2/PI3K-AKT、免疫TGFβ/IFNγ-FAS模块；先固定模块与患者层面统计家族。尚未启动下一步后台作业。'
record={'phase':B.name,'status':'first_pass_completed_verified','completed_at':now,'memo':str(B/'delivery/耐药与抑制研究_首轮结果与后续路线_20261006.txt'),'results_zip':str(B/'delivery/CAF耐药抑制重启_首轮结果_20261006.zip'),'acceptance':str(B/'delivery/FINAL_ACCEPTANCE.json'),'next_action':next_action,'remaining_research_work':['CAF-tumor receiver mechanism analysis','immune suppression mechanism analysis','experimental inhibition validation'],'background_process_running':False}
(B/'PHASE70_RESTART_COMPLETED.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
backup=B/'resources/state_before_restart';backup.mkdir(exist_ok=True)
for name in states:
 shutil.copy2(root/name,backup/name)
 d=json.loads((root/name).read_text(encoding='utf-8'))
 d['previous_completed_task']=d.get('latest_task');d['previous_user_priority']=d.get('user_priority')
 d.update({'current_phase':B.name,'status':'resistance_research_resumed_first_pass_completed','current_state':'resistance_research_resumed_first_pass_completed','updated_at':now,'last_activity_at':now,'user_priority':'已提交胶原基质与GSEA合并报告；重新研究耐药与抑制机制','active_research_focus':['CAF基质来源与化疗病理应答','肿瘤细胞接收端药物应答机制','免疫抑制及解除机制','基质干预的候选靶点与因果验证'],'deferred_objectives':[],'deferred_objectives_policy':'2026-10-06用户明确恢复耐药与抑制研究；旧阴性及敏感性分析保留。','next_action':next_action,'active_task':None,'latest_task':record,'latest_research_results':record,'last_completed_analysis':record,'analysis_now_completed':True,'analysis_completion_scope':'Phase70 first-pass models only; overall resistance/inhibition mechanism research remains open','acceptance_passed':True,'process_observation':'Local analysis exited 0; completion marker and outputs verified; no background analysis job running.','execution_status':{'exit_code':0,'at':now,'host':'local Windows analysis; results synced to configured sgs10'},'resume_record':str(B/'PHASE70_RESTART_COMPLETED.json'),'formal_client_delivery_baseline':str(report),'client_delivery_confirmation':{'confirmed_by':'user','statement':'这部分已经提交给客户，现在重新开始研究耐药抑制','recorded_at':now,'sha256':accept['original_client_report_sha256'],'actual_submission_time':'not independently recorded'},'research_candidate_priorities':['F2_C1 ACAN: one-cohort clinical nominal signal','F2_C1/F2_C3 COL6A2: one-cohort nominal; second cohort opposite direction','LOXL2: mechanistic intervention hypothesis; no clinical response association','ABCG2 / PI3K-AKT tumor receiver and TGFβ / IFNγ-FAS immune mechanisms: literature-directed pending analysis']})
 for key in ['current_stage','phase']:
  if key in d:d[key]=B.name
 for key in ['current_status','task_status']:
  if key in d:d[key]='resistance_research_resumed_first_pass_completed'
 if 'followup_recommendations' in d:d['followup_recommendations']=[next_action]
 tmp=root/(name+'.phase70.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(root/name)
for name in states:
 d=json.loads((root/name).read_text(encoding='utf-8'));assert d['current_phase']==B.name and d['deferred_objectives']==[] and d['latest_report_revision']['report']==str(report)
print(json.dumps({'status':'synced_verified_state_updated','host':socket.gethostname(),'phase':B.name,'files_verified':len(manifest),'report_preserved':True,'background_process_running':False},ensure_ascii=False))
