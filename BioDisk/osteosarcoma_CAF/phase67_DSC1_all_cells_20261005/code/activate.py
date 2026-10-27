from pathlib import Path
import json,datetime,shutil
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase67_DSC1_all_cells_20261005'
q=json.loads((O/'QA/ACCEPTANCE.json').read_text());assert q['passed']
local=json.loads((O/'QA/LOCAL_DELIVERY_QA.json').read_text());assert local['passed']
receipt=json.loads((O/'delivery/DELIVERY_RECEIPT.json').read_text());r=json.loads((O/'resources/SUMMARY.json').read_text())
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
task={'title':'SR1004对应靶点DSC1全细胞表达核查','phase':O.name,'status':'completed_verified_delivered','completed_at':now,'scope':'9个来源、107份样本、876648个source-called细胞，原始计数、覆盖、来源/样本/谱系/候选池内外分布；描述性分析。','summary':r,'report':receipt['report'],'results_zip':receipt['zip'],'acceptance':str(O/'QA/ACCEPTANCE.json'),'execution_status':str(O/'QA/EXECUTION_STATUS.json'),'remaining_required_work':[]}
for name in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
 p=B/name;d=json.loads(p.read_text())
 d.update({'current_phase':O.name,'status':'DSC1_all_cells_verified_delivered','task_status':'DSC1_all_cells_verified_delivered','active_task':None,'last_completed_task':task,'updated_at':now,'next_action':'本次全细胞DSC1核查已完成并交付。对照客户湿实验仍需具体厂家货号、实验方法及样本/细胞范围；未开展蛋白验证。后续8174细胞分析沿用Phase66修复数据入口。'})
 if name=='CAF_NEXT_STAGE_CURRENT.json':
  d['previous_completed_task']=json.loads((O/'resources/PREVIOUS_CAF_NEXT_STAGE_CURRENT.json').read_text()).get('last_completed_task')
  d.update({'report':receipt['report'],'results_zip':receipt['zip'],'acceptance':str(O/'QA/ACCEPTANCE.json'),'local_delivery_verified':True,'all_cell_DSC1':r})
  d['SR1004']['all_cell_followup']=r
 if name=='PROJECT_STATE.json':
  d.update({'current_stage':O.name,'current_status':'DSC1_all_cells_verified_delivered','latest_task':task,'updated_at_utc':now})
 if name=='CAF_RESUME_STATUS_CURRENT.json':
  d.update({'phase':O.name,'current_state':'DSC1_all_cells_verified_delivered','latest_task':task,'checked_at':now,'execution_status':json.loads((O/'QA/EXECUTION_STATUS.json').read_text()),'acceptance_passed':True,'new_analysis_started':True,'analysis_now_completed':True,'process_observation':'Native SSH analysis session finished exit 0; recorded exit 0 and completion marker and output acceptance verified.'})
  dated=B/('CAF_RESUME_STATUS_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'.json');d['resume_record']=str(dated);dated.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
 p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
(O/'PHASE67_DELIVERED.json').write_text(json.dumps({'at':now,'local_delivery_verified':True,'progress_entries_synchronized':True,'metadata_evidence_source':'Phase66 unchanged'},indent=2))
for name in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
 d=json.loads((B/name).read_text());assert d['current_phase']==O.name and d['active_task'] is None
print('PHASE67 delivered; all 3 current-progress entries agree.')
