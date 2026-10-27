from pathlib import Path
import os,sys,subprocess,json,datetime,traceback
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase46_spatial_functional_20261001';E=P/'phase47_external_candidate_integration_20261001';R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
def status(state,**kw):(E/'PIPELINE_STATUS.json').write_text(json.dumps({'state':state,'at':datetime.datetime.now().astimezone().isoformat(),**kw},indent=2))
if '--start' in sys.argv:
 f=(E/'logs/final_controller.log').open('w');job=subprocess.Popen([sys.executable,__file__],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL);(E/'controller.pid').write_text(str(job.pid));print('FINAL_VALIDATION_STARTED',job.pid)
else:
 rc=1
 try:
  jobs=[('external_supervised',E,[sys.executable,str(E/'code/external.py')]),('source_projection',E,[sys.executable,str(E/'code/source_projection.py')]),('spatial_sensitivity',O,[sys.executable,str(O/'code/sensitivity.py')]),('QA_prepare',O,[sys.executable,str(O/'code/prepare_QA.py')]),('independent_R',O,[R,'--vanilla',str(O/'code/independent_QA.R')])]
  for name,out,cmd in jobs:
   status('running',step=name)
   with (out/'logs'/f'{name}.log').open('w') as f:
    job=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT);(out/f'{name}.pid').write_text(str(job.pid));c=job.wait()
   (out/f'{name}.exit_code').write_text(str(c));assert c==0,(name,c)
  for path in [O/'SPATIAL_COMPLETE',O/'SENSITIVITY_COMPLETE',O/'INDEPENDENT_R_QA_COMPLETE',E/'EXTERNAL_COMPLETE',E/'SOURCE_PROJECTION_COMPLETE',E/'INDEPENDENT_R_QA_COMPLETE']:assert path.read_text().strip()=='SUCCESS'
  (E/'PIPELINE_COMPLETE').write_text('SUCCESS');(O/'PIPELINE_COMPLETE').write_text('SUCCESS');status('complete',steps=[r[0] for r in jobs]);rc=0
 except Exception as e:status('failed',error=repr(e));traceback.print_exc()
 (E/'controller.exit_code').write_text(str(rc));print('FINAL_CONTROLLER_EXIT',rc,flush=True)
