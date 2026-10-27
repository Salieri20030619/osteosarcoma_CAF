from pathlib import Path
import os,sys,subprocess,fcntl,json,datetime,traceback
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase25_complete_nonPDX_assets_20260921');os.chdir(O)
f=(O/'runner.lock').open('w');fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
(O/'RUN_STARTED').write_text(str(os.getpid())+'\n');rc=1
try:
 for script,log in [('range_zip.py','archive_index.log'),('download_parallel.py','download.log'),('finalize_official_download.py','verification.log')]:
  with (O/log).open('a') as output:rc=subprocess.run([sys.executable,'-u',script],stdout=output,stderr=subprocess.STDOUT).returncode
  (O/(script+'.exit_code')).write_text(str(rc)+'\n')
  if rc:break
except Exception as e:
 (O/'RUNNER_ERROR').write_text(type(e).__name__);rc=1
finally:
 (O/'exit_code').write_text(str(rc)+'\n');(O/'RUN_FINISHED').write_text(datetime.datetime.now(datetime.timezone.utc).isoformat()+'\n')
sys.exit(rc)
