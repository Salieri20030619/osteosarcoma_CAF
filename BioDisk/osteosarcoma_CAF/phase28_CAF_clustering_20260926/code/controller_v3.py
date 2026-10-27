from pathlib import Path
import os,select,subprocess,json,traceback,datetime,ctypes
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');b=P/'phase27_identity_review_20260926';o=P/'phase28_CAF_clustering_20260926'
def waitpid(pid):
 try:
  if hasattr(os,'pidfd_open'):fd=os.pidfd_open(pid)
  else:
   lib=ctypes.CDLL(None,use_errno=True);fd=lib.syscall(434,int(pid),0)
   if fd<0:
    if ctypes.get_errno()==3:return
    raise OSError(ctypes.get_errno(),'pidfd_open')
  select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
def status(s,**kw):(o/'CONTROLLER_STATUS_v3.json').write_text(json.dumps(dict(status=s,at=datetime.datetime.now().astimezone().isoformat(),**kw),indent=2))
rc=1
try:
 status('WAITING_FOR_REVIEW_AND_VERIFICATION')
 waitpid(int((b/'runner.pid').read_text()))
 waitpid(int((o/'controller.pid').read_text()))
 raw_rc=int((b/'exit_code').read_text());raw_status=json.loads((b/'STATUS.json').read_text())
 assert (b/'REVIEW_COMPLETE').exists()
 if raw_rc!=0:
  assert 'multiple values for keyword argument' in raw_status.get('error','') and 'at' in raw_status['error'],raw_status
 import pandas as pd,anndata as ad
 qa=json.loads((b/'delivery/IDENTITY_QA.json').read_text());tab=pd.read_csv(b/'tables/sample_review_counts.tsv',sep='\\t')
 a=ad.read_h5ad(b/'objects/CAF_candidates_counts.h5ad',backed='r')
 assert len(tab)==107 and tab['sample'].is_unique and tab.input_cells.sum()==876648
 assert a.n_obs==qa['CAF_candidates']==tab.CAF_candidates.sum() and a.obs_names.is_unique
 assert a.obs.CNV_evidence.isin(['single_method_no_alteration','concordant_no_alteration']).all() and not a.obs.known_doublet.any()
 for r in tab.to_dict('records'):
  df=pd.read_csv(b/'tables'/(r['sample']+'_review.tsv.gz'),sep='\\t',usecols=['cell_id','CAF_candidate'])
  assert len(df)==r['input_cells'] and df.cell_id.is_unique and int(df.CAF_candidate.sum())==r['CAF_candidates']
 a.file.close()
 acceptance={'status':'IDENTITY_REVIEW_ACCEPTED','all107_cell_counts_and_unique_ids':'PASS','candidate_object_and_CNV_identity_integrity':'PASS','original_process_exit_code':raw_rc,'original_final_status':raw_status,'reporting_error_repaired':raw_rc!=0,'scientific_outputs_recomputed':False}
 (b/'delivery/ACCEPTANCE_QA.json').write_text(json.dumps(acceptance,indent=2))
 (b/'REVIEW_ACCEPTED').write_text('ALL107_IDENTITY_OUTPUTS_VERIFIED\\n');(b/'acceptance.exit_code').write_text('0\\n')
 (b/'STATUS.before_acceptance.json').write_text(json.dumps(raw_status,indent=2));(b/'STATUS.json').write_text(json.dumps(acceptance,indent=2))
 status('RUNNING_CAF_CLUSTERING',review_QA=qa)
 with (o/'logs/cluster.log').open('x') as log:
  r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(o/'code/cluster.py')],stdout=log,stderr=subprocess.STDOUT)
 rc=r.returncode
 assert rc==0 and (o/'CLUSTERING_COMPLETE').exists() and (o/'delivery/CLUSTER_QA.json').exists()
 status('CAF_CLUSTERING_COMPLETE_ANNOTATION_REVIEW_REQUIRED')
except Exception as e:status('CONTROLLER_FAILED',error=str(e));traceback.print_exc()
finally:(o/'controller_v3.exit_code').write_text(str(rc)+'\\n')
