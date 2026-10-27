from pathlib import Path
import os,subprocess,json,hashlib,datetime,fcntl
p=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');o=p/'phase95_NCI_expanded_reanalysis_20261016'
lock=(o/'export.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
f=p/'patient_expansion_download_investigation_20261013_v1/downloads/OsteoCAR/patient_mets.qs2';h=hashlib.sha256()
with f.open('rb') as r:
 for b in iter(lambda:r.read(8*1024*1024),b''):h.update(b)
assert h.hexdigest()=='37e10c55855af596ecb19666b4e9fb0775a6b4bb5c36d0339eb88000fc132856'
(o/'resources/SOURCE_SHA256.json').write_text(json.dumps({'path':str(f),'sha256':h.hexdigest(),'matches_original_download':True}),encoding='utf-8')
e='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2';env=os.environ.copy();env.update(PATH=e+'/bin:'+env.get('PATH',''),CONDA_PREFIX=e,LD_LIBRARY_PATH=e+'/lib:'+env.get('LD_LIBRARY_PATH',''),OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
r=subprocess.run([e+'/bin/Rscript','--vanilla',str(o/'code/export_inputs.R')],env=env)
(o/'logs/export.exit_code').write_text(str(r.returncode),encoding='utf-8')
(o/'resources/EXPORT_JOB_FINISHED.json').write_text(json.dumps({'exit_code':r.returncode,'completion_marker':(o/'resources/EXPORT_COMPLETE').exists(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}),encoding='utf-8')
