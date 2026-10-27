
from pathlib import Path
import ctypes,os,subprocess,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1');marker=O/'resources/OsteoCAR_patient_prim_DOWNLOAD_COMPLETE.json'
lib=ctypes.CDLL('libc.so.6',use_errno=True);fd=lib.inotify_init1(os.O_CLOEXEC);assert fd>=0
wd=lib.inotify_add_watch(fd,str(O/'resources').encode(),0x00000008|0x00000080);assert wd>=0
print('WAITING_FOR_PRIMARY_FILE_COMPLETION_EVENT',flush=True)
while not marker.exists():os.read(fd,4096)
os.close(fd);receipt=json.loads(marker.read_text());assert receipt['official_md5_verified'] is True
assert (O/'logs/qs2_install_fixed.exit_code').read_text().strip()=='0'
E='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2';env=os.environ.copy();env.update(PATH=E+'/bin:'+env.get('PATH',''),CONDA_PREFIX=E,LD_LIBRARY_PATH=E+'/lib:'+env.get('LD_LIBRARY_PATH',''))
r=subprocess.run([E+'/bin/Rscript','--vanilla',str(O/'code/extract_OsteoCAR_prim.R')],env=env);(O/'logs/OsteoCAR_prim_extraction.exit_code').write_text(str(r.returncode)+'\n')
