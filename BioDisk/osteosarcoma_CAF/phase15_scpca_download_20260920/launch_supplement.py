from pathlib import Path
import subprocess, json, datetime, shutil, os
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase16_scpca_supplement_20260920'
for folder in ['code','logs']: (O/folder).mkdir(parents=True,exist_ok=True)
script=P/'phase15_scpca_download_20260920/preprocess_supplement.py'
shutil.copy2(script,O/'code/preprocess.py')
py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
subprocess.run([py,'-m','py_compile',str(O/'code/preprocess.py')],check=True)
runner=O/'code/run.sh'
runner.write_text('#!/usr/bin/env bash\nset -uo pipefail\ncd '+str(O)+'\n'+py+' -u code/preprocess.py > logs/preprocess.log 2>&1\nrc=$?\nprintf \'%s\\n\' "$rc" > logs/preprocess.exit_code\nexit "$rc"\n')
log=(O/'logs/launcher.log').open('a')
proc=subprocess.Popen(['bash',str(runner)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(O/'logs/runner.pid').write_text(str(proc.pid))
backup=P/'PROJECT_STATE.before_phase16_20260920.json'
if not backup.exists():shutil.copy2(P/'PROJECT_STATE.json',backup)
state=json.loads((P/'PROJECT_STATE.json').read_text())
state.update(current_stage=O.name,status='SCPCA_SUPPLEMENT_QC_RUNNING',current_status='SCPCA_SUPPLEMENT_QC_RUNNING',dynamic_status_file=str(O/'STATUS.json'),updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
state[O.name]={'samples':12,'pipeline':'QC then lineage/ambient/dual-CNV acceptance then client panels','pid':proc.pid,'strict_CAF_confirmed':None,'cross_study_independent_N':'unresolved','user_authorized':'2026-09-20 supplemental analysis'}
state['next_action']='Review ScPCA QC and lineages; ambient correction and dual CNV; client-directed paired descriptive analysis after identity acceptance'
(P/'PROJECT_STATE.json').write_text(json.dumps(state,ensure_ascii=False,indent=2))
print('LAUNCHED',proc.pid,str(O))
