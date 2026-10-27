from pathlib import Path
import os,subprocess,json,datetime,traceback,hashlib
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase44_patient_mechanism_20261001');R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
state={'pid':os.getpid(),'started':datetime.datetime.now().astimezone().isoformat(),'status':'running','reason':'fix double close of already-managed R connection; previous failure log retained; statistics not recomputed','outcomes':[]}
rc=1
try:
 for script,label,marker in [('stageD_independent_QA_v2.R','independent_R_QA_v2',O/'INDEPENDENT_R_QA_COMPLETE'),('stageD_reporting_QA.R','reporting_extension_R_QA',O/'clinical_reporting_extension/INDEPENDENT_R_QA_COMPLETE')]:
  with (O/'logs'/f'{label}.log').open('x') as log:
   child=subprocess.Popen(['/usr/bin/time','-v',R,'--vanilla',str(O/'code'/script)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT);state['active_job']=label;state['active_pid']=child.pid;(O/'QA_RESUME_STATUS.json').write_text(json.dumps(state,indent=2));code=child.wait()
  (O/f'{label}.exit_code').write_text(str(code)+'\n');state['outcomes'].append({'job':label,'exit_code':code,'success_marker':marker.exists()});(O/'QA_RESUME_STATUS.json').write_text(json.dumps(state,indent=2))
  assert code==0 and marker.read_text().strip()=='SUCCESS',label
  print('ACCEPTED',label,flush=True)
 state.update(status='complete',completed=datetime.datetime.now().astimezone().isoformat(),active_job=None,active_pid=None);(O/'QA_RESUME_STATUS.json').write_text(json.dumps(state,indent=2));(O/'CORE_COMPLETE').write_text('SUCCESS\n');(O/'delivery/CORE_ACCEPTED_STATUS.json').write_text(json.dumps({'state':'accepted','raw_statistics_exit_codes':[0,0],'independent_R_exit':0,'reporting_extension_R_exit':0,'failure_history':'CORE_PIPELINE_STATUS.json and independent_R_QA.log preserved','resume_status':'QA_RESUME_STATUS.json','scope':'stability, TF signed target response, clinical fractions; communications and final C2 synthesis pending'},indent=2));print('STAGE_D_CORE_ACCEPTED_AFTER_QA_REPAIR',flush=True);rc=0
except Exception:traceback.print_exc();state['status']='failed';(O/'QA_RESUME_STATUS.json').write_text(json.dumps(state,indent=2))
finally:(O/'qa_resume.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)