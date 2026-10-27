from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,datetime,fcntl,sys
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase90_patient_validation_20261015'
review_lock=(O/'native_evidence_review.lock').open('a');fcntl.flock(review_lock,fcntl.LOCK_EX)
completed_now=len(list((O/'cnv/infercnv').glob('NCC__*/COMPLETE')))
if (O/'resources/NATIVE_REVIEW_PROGRESS.json').exists() and '--force' not in sys.argv:
 saved=json.loads((O/'resources/NATIVE_REVIEW_PROGRESS.json').read_text(encoding='utf-8'))
 if saved['completed_samples_reviewed']==completed_now:
  print(json.dumps(saved,indent=2));raise SystemExit(0)
old=pd.read_csv(P/'audit_expanded_admission_20261015/tables/final_cell_review_decisions.tsv.gz',sep='\t').set_index('cell_id')
manifest=pd.read_csv(O/'resources/manifest.tsv',sep='\t')
rows=[];allc=[];cal=[]
for r in manifest.to_dict('records'):
 sm=r['sample'];cnv=O/'cnv/infercnv'/sm
 if not (cnv/'COMPLETE').exists():continue
 exitfile=O/'logs'/(sm+'_infercnv.exit_code')
 if not exitfile.exists() or exitfile.read_text(encoding='utf-8').strip()!='0':continue
 d=O/'inputs'/sm
 cells=pd.read_csv(d/'cells.tsv',sep='\t').set_index('cell_id')
 calls=pd.read_csv(cnv/'cell_results.tsv',sep='\t').set_index('cell_id')
 dbl=pd.read_csv(d/'full_scDblFinder_calls.tsv',sep='\t').set_index('cell_id')
 assert calls.index.is_unique and dbl.index.is_unique
 assert set(cells.index)==set(calls.index)==set(dbl.index)
 cells['native_HMM_call']=calls['call'];cells['legacy_main_HMM_call']=calls['call'];cells['full_scDblFinder_class']=dbl.scDblFinder_class
 cells['native_reference_used']=cells.reference
 primary_path=cnv/'cell_results.tsv';variant='main_native_run'
 membership=cnv/'reference_membership.tsv'
 if membership.exists():
  rm=pd.read_csv(membership,sep='\t',encoding='utf-8').set_index('cell_id');assert set(rm.index)==set(cells.index)
  cells['native_reference_used']=rm.reference_used.astype(str).str.lower().eq('true')
  variant='prospective_clean_reference_native_run'
 clean_receipt=O/'resources/BS7_CLEAN_REFERENCE_EXIT.json'
 if sm=='NCC__BS7' and clean_receipt.exists() and json.loads(clean_receipt.read_text(encoding='utf-8')).get('exit_code')==0 and (O/'cnv_clean/infercnv'/sm/'COMPLETE').exists():
  primary_path=O/'cnv_clean/infercnv'/sm/'cell_results.tsv'
  clean_calls=pd.read_csv(primary_path,sep='\t',encoding='utf-8').set_index('cell_id');assert clean_calls.index.is_unique and set(clean_calls.index)==set(cells.index)
  cells['native_HMM_call']=clean_calls['call']
  rm=pd.read_csv(O/'inputs_clean'/sm/'cells.tsv',sep='\t',encoding='utf-8').set_index('cell_id')
  cells['native_reference_used']=rm.reference
  variant='clean_reference_primary_with_legacy_sensitivity'
 prior=pd.read_csv(P/'phase88_expansion_identity_20261014/copykat'/sm/'cell_results.tsv',sep='\t').set_index('cell_id')
 cells['prior_subset_native_CopyKAT_call']=prior['call'].reindex(cells.index).fillna('not_in_prior_subset')
 cells['prior_expanded_retained']=old.expanded_review_retained.reindex(cells.index).fillna(False).astype(bool)
 cells['prior_supported_CAF']=old.supported_CAF_candidate.reindex(cells.index).fillna(False).astype(bool)
 cells['prior_reference_pass']=old.reference_QC_pass.reindex(cells.index).fillna(False).astype(bool)
 for col in ['manual_cluster_decision','manual_reason','cluster','infer_call','copy_call','decontx_contamination','corrected_lineage_supported','CTA_context','mixed_lineage_flag']:
  cells['prior_'+col]=old[col].reindex(cells.index)
 cells['known_doublet_any']=cells.full_scDblFinder_class.eq('doublet')|old.audit_doublet_positive.reindex(cells.index).fillna(False).astype(bool)
 normal=cells.native_HMM_call.eq('no_alteration_called')|cells.prior_subset_native_CopyKAT_call.eq('diploid')
 altered=cells.native_HMM_call.eq('alteration_called')|cells.prior_subset_native_CopyKAT_call.eq('aneuploid')
 if variant=='clean_reference_primary_with_legacy_sensitivity':altered=altered|cells.legacy_main_HMM_call.eq('alteration_called')
 cells['native_entry_no_opposing_alteration']=normal&~altered
 cells['full_input_mechanical_candidate']=cells.provisional_lineage.eq('Fibroblast')&normal&~altered&~cells.known_doublet_any
 cells['reviewed_native_candidate']=cells.full_input_mechanical_candidate&cells.prior_expanded_retained
 cells['full_native_HMM_supported_candidate']=cells.reviewed_native_candidate&cells.native_HMM_call.eq('no_alteration_called')
 cells['supported_native_candidate']=cells.full_native_HMM_supported_candidate&cells.prior_supported_CAF
 # This is an evidence review, not automatic patient admission.
 cells['formal_admission']=False
 hold_original=cells[cells.heldout];hold=hold_original[~hold_original.known_doublet_any]
 h_eval=hold.native_HMM_call.isin(['no_alteration_called','alteration_called'])
 h_n=int(h_eval.sum());h_normal=int(hold.native_HMM_call.eq('no_alteration_called').sum())
 for group,g0 in hold_original.groupby('phase88_cnv_role'):
  g=g0[~g0.known_doublet_any]
  e=g.native_HMM_call.isin(['no_alteration_called','alteration_called']);nn=int(e.sum());normal_n=int(g.native_HMM_call.eq('no_alteration_called').sum())
  cal.append({'sample':sm,'reference_holdout_group':group,'original_cells':len(g0),'new_doublets_excluded':int(g0.known_doublet_any.sum()),'clean_cells':len(g),'native_evaluable':nn,'native_no_alteration':normal_n,'native_normal_fraction':normal_n/nn if nn else None})
 # Explicitly assess new-positive references, without silently discarding them.
 ref=cells[cells.native_reference_used]
 z={'sample':sm,'patient_id':r['patient_id'],'full_cells':len(cells),'prior_subset_cells':int(cells.phase88_retained.sum()),'full_doublet_cells':int(cells.known_doublet_any.sum()),'native_HMM_no_alteration':int(cells.native_HMM_call.eq('no_alteration_called').sum()),'native_HMM_alteration':int(cells.native_HMM_call.eq('alteration_called').sum()),'native_HMM_unassessed':int((~cells.native_HMM_call.isin(['no_alteration_called','alteration_called'])).sum()),'full_mechanical_candidates':int(cells.full_input_mechanical_candidate.sum()),'previously_reviewed_native_candidates':int(cells.reviewed_native_candidate.sum()),'previously_supported_native_candidates':int(cells.supported_native_candidate.sum()),'new_candidate_cells_pending_cluster_review':int((cells.full_input_mechanical_candidate&~cells.prior_expanded_retained).sum()),'heldout_evaluable':h_n,'heldout_normal':h_normal,'heldout_normal_fraction':h_normal/h_n if h_n else None,'reference_cells':len(ref),'reference_new_doublet_positive':int(ref.known_doublet_any.sum()),'formal_admission':False}
 z.update(primary_native_HMM_file=str(primary_path),reference_variant=variant,full_native_HMM_supported_candidates=int(cells.full_native_HMM_supported_candidate.sum()),original_reference_new_doublet_positive=int((cells.reference&cells.known_doublet_any).sum()),heldout_original_cells=len(hold_original),heldout_new_doublets_excluded=int(hold_original.known_doublet_any.sum()))
 rows.append(z);allc.append(cells.reset_index())
 cells.to_csv(O/'tables'/(sm+'_full_native_review.tsv.gz'),sep='\t',index_label='cell_id',encoding='utf-8')
pd.DataFrame(rows).to_csv(O/'tables/native_sample_review.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(cal).to_csv(O/'tables/native_holdout_calibration.tsv',sep='\t',index=False,encoding='utf-8')
summary={'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'completed_samples_reviewed':len(rows),'samples':rows,'formal_admission':False,'note':'Counts are review evidence only. Patient provenance, reference calibration, new doublet flags and cluster evidence need adjudication before admission.'}
(O/'resources/NATIVE_REVIEW_PROGRESS.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
