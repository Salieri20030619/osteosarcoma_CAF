from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase86_functional_bone_lung_20261012'
def dump(x,p):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
r=json.loads((O/'DELIVERY_READY.json').read_text());assert r['QA']['status']=='PASS'
assert (O/'logs/analysis.exit_code').read_text().strip()=='0'
visual=dict(status='PASS',reviewed_figures=['01_functional_bone_lung.png','02_patient_scores.png'],review='Both local rendered figures inspected: labels, statistics, clinical denominators and patient-level display legible; no clipping.',checked_at=datetime.now(timezone.utc).isoformat())
dump(visual,O/'QA/VISUAL_QA.json')
r['status']='completed_verified';r['visual_QA']=visual;r['completed_at']=datetime.now(timezone.utc).isoformat();r['local_delivery_verified']=False
for k in ['archive','archive_sha256','archive_bytes']:r.pop(k,None)
dump(r,O/'delivery/SUMMARY.json')
files=[q for folder in ['tables','resources','QA','figures','delivery','code'] for q in (O/folder).rglob('*') if q.is_file() and q.suffix not in ['.zip','.pyc'] and q.name!='FILE_MANIFEST.json']
mf={str(q.relative_to(O)):dict(bytes=q.stat().st_size,sha256=sha(q)) for q in files};dump(mf,O/'delivery/FILE_MANIFEST.json')
archive=O/'delivery/CAF_functional_bone_lung_20261012_verified.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for name in mf:z.write(O/name,name)
 z.write(O/'delivery/FILE_MANIFEST.json','delivery/FILE_MANIFEST.json')
r.update(archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size)
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
for n in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json','CAF_ALL14_RESEARCH_CURRENT.json']:
 path=P/n;shutil.copy2(path,O/'state_backups'/f'{stamp}_{n}');obj=json.loads(path.read_text())
 obj['functional_bone_lung_analysis']=r;obj['clinical_comparison_framework']='Bone versus Lung anatomical clinical groups; source datasets analyzed separately; A/B1/B2 kept as prior context'
 if n in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
  for k in ['current_phase','current_stage','phase']:obj[k]=O.name
  for k in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:obj[k]=r
  obj['active_task']=dict(phase=O.name,status='completed',description='临床骨—肺前提下的功能状态分析、表面候选关联与QA已完成')
  obj['next_analysis_plan']='现有预设功能程序未检出骨—肺名义差异；保留全部14个候选及功能状态标志。复发/耐药和蛋白流式验证尚未完成，不混入本轮结论。'
 if n=='CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json':
  obj['new_state_analysis_status']='completed_phase86';obj['status']='objective_and_current_analysis_completed';obj['new_flow_panel_status']='not_validated'
 if isinstance(obj.get('current_research_objective'),dict):obj['current_research_objective']['new_state_analysis_status']='completed_phase86'
 dump(obj,path)
dump(r,P/'CAF_FUNCTIONAL_BONE_LUNG_CURRENT.json');dump(r,O/'COMPLETED.json')
print(json.dumps(r,ensure_ascii=False))
