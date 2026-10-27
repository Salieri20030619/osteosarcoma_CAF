from pathlib import Path
import os,select,time,subprocess,sys,ctypes,platform,errno
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_identity_20260919')
pid=int((O/'logs/pipeline.pid').read_text())
if not (O/'logs/pipeline.exit_code').exists():
 fd=-1
 if hasattr(os,'pidfd_open'):
  try:fd=os.pidfd_open(pid)
  except ProcessLookupError:pass
 elif platform.machine()=='x86_64':
  libc=ctypes.CDLL(None,use_errno=True);fd=libc.syscall(434,pid,0)
  if fd<0 and ctypes.get_errno() not in [errno.ESRCH,errno.ENOSYS]:raise OSError(ctypes.get_errno())
 if fd>=0:
  print('NATIVE_PIDFD_WAIT',pid,flush=True)
  select.select([fd],[],[]);os.close(fd)
 else:
  while Path('/proc/'+str(pid)).exists() and not (O/'logs/pipeline.exit_code').exists():time.sleep(900)
assert (O/'logs/pipeline.exit_code').exists(),'Process ended without exit record'
assert (O/'logs/pipeline.exit_code').read_text().strip()=='0','Ambient pipeline failed'
assert (O/'AMBIENT_COMPLETE').exists()
with open(O/'logs/identity.log','w') as f:
 rc=subprocess.call(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/integrate.py')],stdout=f,stderr=subprocess.STDOUT)
(O/'logs/identity.exit_code').write_text(str(rc)+'\n')
print('IDENTITY_EXIT',rc,flush=True)
sys.exit(rc)
