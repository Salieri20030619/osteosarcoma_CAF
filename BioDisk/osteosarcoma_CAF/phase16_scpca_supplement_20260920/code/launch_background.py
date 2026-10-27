from pathlib import Path
import subprocess,fcntl,json,datetime
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase16_scpca_supplement_20260920')
lock=(O/'logs/background.lock').open('a')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:
 print('ALREADY_RUNNING');raise SystemExit(0)
py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
for name in ['prepare_followup','finalize_followup','background_pipeline']:
 subprocess.run([py,'-m','py_compile',str(O/'code'/f'{name}.py')],check=True)
subprocess.run(['bash','-n',str(O/'code/run_cnv_job.sh')],check=True)
fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
with (O/'logs/background.log').open('a') as f:
 p=subprocess.Popen([py,'-u',str(O/'code/background_pipeline.py')],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
(O/'logs/background.pid').write_text(str(p.pid)+'\n')
print('BACKGROUND_LAUNCHED',p.pid,'native subprocess completion chaining; SSH-independent')
