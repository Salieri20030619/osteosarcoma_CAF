from pathlib import Path
import subprocess,json,datetime,os,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase45_LR_expression_20261001');R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript';rc=1
try:
 assert (O/'LR_expression.exit_code').read_text().strip()=='0' and (O/'COMPLETE').read_text().strip()=='SUCCESS'
 with (O/'logs/independent_R_v2.log').open('x') as log:
  p=subprocess.Popen(['/usr/bin/time','-v',R,'--vanilla',str(O/'code/stageD_LR_QA_v2.R')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT);(O/'QA_RESUME_STATUS.json').write_text(json.dumps({'pid':os.getpid(),'active_pid':p.pid,'state':'running','at':datetime.datetime.now().astimezone().isoformat(),'reason':'Python True/False CSV values explicitly parsed as R logical; previous failed QA log retained'},indent=2));code=p.wait()
 (O/'independent_R_v2.exit_code').write_text(str(code)+'\n');assert code==0 and (O/'INDEPENDENT_R_QA_COMPLETE').read_text().strip()=='SUCCESS'
 (O/'PIPELINE_COMPLETE').write_text('SUCCESS\n');(O/'delivery/QA_ACCEPTED_STATUS.json').write_text(json.dumps({'state':'accepted','expression_exit_code':0,'independent_R_v2_exit_code':0,'original_failed_controller_status_retained':True,'current_exit_code':'qa_resume.exit_code','at':datetime.datetime.now().astimezone().isoformat()},indent=2));(O/'QA_RESUME_STATUS.json').write_text(json.dumps({'state':'complete','at':datetime.datetime.now().astimezone().isoformat(),'expression_not_recomputed':True},indent=2));print('PATIENT_LR_PIPELINE_ACCEPTED_AFTER_QA_REPAIR',flush=True);rc=0
except Exception:traceback.print_exc();(O/'QA_RESUME_STATUS.json').write_text(json.dumps({'state':'failed','error':traceback.format_exc()},indent=2))
finally:(O/'qa_resume.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)