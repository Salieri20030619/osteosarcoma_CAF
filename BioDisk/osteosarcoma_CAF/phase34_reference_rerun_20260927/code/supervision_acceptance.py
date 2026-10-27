
from pathlib import Path
import json,datetime,hashlib,subprocess,traceback,sys
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase34_reference_rerun_20260927';B=P/'phase26_client_reanalysis_20260922'
def now():return datetime.datetime.now().astimezone().isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
kind=sys.argv[1]
try:
 if kind=='rerun':
  assert (O/'FINISHED').exists() and (O/'exit_code').read_text().strip()=='0','Controller failed or completion record absent'
  protocol=json.loads((O/'resources/RERUN_PROTOCOL.json').read_text());s=json.loads((O/'STATUS.json').read_text())
  assert len(s['outcomes'])==4 and all(x['exit_code']==0 and x['accepted'] for x in s['outcomes'])
  assert sha(O/'code/reference_sensitivity_frozen.R')==protocol['script_sha256']
  assert sha(O/'code/run_frozen.R')==protocol['wrapper_sha256']
  assert all(sha(x['path'])==x['sha256'] for x in protocol['assets'])
  labels={k:pd.read_csv(P/'phase30_two_layer_classification_20260926'/f'tables/{k}_classified_cells.tsv.gz',sep='\t').set_index('cell_id') for k in ['supported','expanded']}
  controls=[];summaries=[];details=[];records=[]
  for job in protocol['jobs']:
   sm=job['sample'];v=job['variant'];d=O/'cnv_reference'/sm/v
   assert (d/'exit_code').read_text().strip()=='0' and (d/'COMPLETE').exists() and not (d/'FAILED').exists()
   cells=pd.read_csv(B/'inputs'/sm/'cells.tsv',sep='\t');x=pd.read_csv(d/'cell_results.tsv',sep='\t')
   assert x.cell_id.is_unique and x.cell_id.tolist()==cells.cell_id.tolist()
   gf=list(d.glob('17_HMM*cell_groupings'));rf=list(d.glob('HMM_CNV_predictions*Pnorm_*pred_cnv_regions.dat'));assert len(gf)==len(rf)==1
   g=pd.read_csv(gf[0],sep='\t');r=pd.read_csv(rf[0],sep='\t');assert g.cell.is_unique and set(r.state)<=set(range(1,7))
   assigned=cells.tool_id.map(g.set_index('cell').cell_group_name);altered=set(r.loc[r.state.ne(3),'cell_group_name'])
   expected=assigned.map(lambda z:'not_returned_by_software' if pd.isna(z) else 'alteration_called' if z in altered else 'no_alteration_called')
   assert expected.tolist()==x.call.tolist()
   m=pd.read_csv(O/'inputs'/sm/'reference_membership.tsv.gz',sep='\t').set_index('cell_id');calls=x.set_index('cell_id').call
   for name,col in [('original_broad_reference','reference'),('immune_endothelial','immune_endothelial_reference'),('lymphoid','lymphoid_reference')]:
    z=calls.loc[m.index[m[col]]];controls.append(dict(sample=sm,variant=v,group=name,cells=len(z),subcluster_non_neutral_region=int(z.eq('alteration_called').sum()),no_non_neutral_region=int(z.eq('no_alteration_called').sum()),not_returned=int(z.eq('not_returned_by_software').sum())))
   for layer,lab in labels.items():
    for target in ['all_candidates','F1','F2','F3','F4','F5']:
     ids=lab.index[lab['sample'].eq(sm)&(True if target=='all_candidates' else lab.class_id.eq(target))];z=calls.reindex(ids)
     summaries.append(dict(sample=sm,variant=v,layer=layer,target=target,cells=len(ids),subcluster_non_neutral_region=int(z.eq('alteration_called').sum()),no_non_neutral_region=int(z.eq('no_alteration_called').sum()),not_returned=int(z.eq('not_returned_by_software').sum())))
   x['sample']=sm;x['variant']=v;x['provider_annotation']=m.provider_annotation.reindex(x.cell_id).to_numpy();details.append(x)
   previous=P/'phase33_client_report_20260926/cnv_reference'/sm/v/'cell_results.tsv'
   same=None
   if previous.exists():
    old=pd.read_csv(previous,sep='\t').set_index('cell_id').call;assert set(old.index)==set(calls.index);same=int(calls.eq(old.reindex(calls.index)).sum())
   records.append(dict(sample=sm,variant=v,cells=len(cells),exact_reconstruction=True,sha256=sha(d/'cell_results.tsv'),same_legacy_calls=same))
  # A separate immutable R script checks final serialized objects and actual reference membership.
  rs=O/'code/supervision_validate_final_RDS.R'
  assert not rs.exists()
  rs.write_text("suppressPackageStartupMessages(library(infercnv))\nO<-"+json.dumps(str(O))+"\nfor(sm in c('SCPCS000614','SCPCS000610'))for(v in c('immune_endothelial_reference','lymphoid_reference')){x<-readRDS(file.path(O,'cnv_reference',sm,v,'run.final.infercnv_obj'));a<-read.delim(file.path(O,'inputs',sm,paste0(v,'_annotations.tsv')),header=FALSE);stopifnot(inherits(x,'infercnv'),all(is.finite(x@expr.data)),setequal(colnames(x@expr.data),a[[1]]),setequal(colnames(x@expr.data)[unlist(x@reference_grouped_cell_indices,use.names=FALSE)],a[a[[2]]=='reference',1]));cat(sm,v,'FINAL_RDS_AND_REFERENCE_OK\\n');rm(x);gc()}\n")
  rs.chmod(0o444)
  p=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript','--vanilla',str(rs)],text=True,capture_output=True)
  (O/'logs/supervision_final_RDS.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
  pd.DataFrame(controls).to_csv(O/'tables/rerun_reference_control_behavior.tsv',sep='\t',index=False)
  pd.DataFrame(summaries).to_csv(O/'tables/rerun_candidate_evidence_summary.tsv',sep='\t',index=False)
  pd.concat(details,ignore_index=True).to_csv(O/'tables/rerun_all_cell_evidence.tsv.gz',sep='\t',index=False)
  result=dict(at=now(),status='COMPUTATION_ACCEPTED_INTERPRETATION_REVIEW_REQUIRED',jobs=records,all_native_exit_codes_zero=True,inputs_and_frozen_scripts_unchanged=True,final_objects_readable_and_reference_exact=True,biological_identity_confirmed=False,interpretation='Assigned HMM subcluster non-neutral region evidence is not a binary malignancy diagnosis. Normal reference behavior must be considered. Frozen research cell sets unchanged.')
 else:
  S=B/'resource_recovery_SCPCS000525_20260925';d=S/'cnv/copykat/SCPCS000525'
  assert (S/'RUN_FINISHED').exists() and (S/'exit_code').read_text().strip()=='0'
  q=json.loads((S/'delivery/QA.json').read_text());assert q['native_success'] and q['native_exit_code']==0,q
  assert (d/'exit_code').read_text().strip()=='0' and (d/'COMPLETE').exists() and not (d/'FAILED').exists()
  x=pd.read_csv(d/'cell_results.tsv',sep='\t');c=pd.read_csv(B/'inputs/SCPCS000525/cells.tsv',sep='\t')
  assert x.cell_id.is_unique and set(x.cell_id)==set(c.cell_id)
  result=dict(at=now(),status='RECOVERY_OUTPUT_ACCEPTED_REVIEW_REQUIRED',cells=len(x),native_exit_code=0,original_controller_QA=q,sha256=sha(d/'cell_results.tsv'),frozen_layers_not_automatically_replaced=True)
except Exception:
 result=dict(at=now(),status='REVIEW_REQUIRED_FAILURE',error=traceback.format_exc())
f=O/'delivery'/f'SUPERVISION_{kind.upper()}_ACCEPTANCE.json';f.write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False));sys.exit(1 if result['status']=='REVIEW_REQUIRED_FAILURE' else 0)
