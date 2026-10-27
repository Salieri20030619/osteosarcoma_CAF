from pathlib import Path
import json,datetime,hashlib
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase95_NCI_expanded_reanalysis_20261016')
run=json.loads((P/'COMPUTATION_FINISHED.json').read_text(encoding='utf-8'))
assert len(run['samples'])==19
allcounts=[];parts=[];missing_manual=[]
for result in run['samples']:
 sm=result['sample'];f=P/'tables'/(sm+'_all_cell_evidence.tsv.gz')
 if not f.exists():allcounts.append({'sample':sm,'status':result['status'],'cells':None});continue
 c=pd.read_csv(f,sep='\t',encoding='utf-8');assert c.cell_id.is_unique
 m=P/'tables'/(sm+'_manual_cluster_review.tsv')
 c['cluster_review']='pending'
 if m.exists():
  review=pd.read_csv(m,sep='\t',encoding='utf-8');review['cluster']=review.cluster.astype(str)
  assert review.cluster.is_unique
  c['cluster_review']=c.cluster.astype(str).map(review.set_index('cluster').decision).fillna('pending')
 if c.cluster_review.eq('pending').any():missing_manual.append(sm)
 c['reviewed_expanded_candidate']=c.mechanical_expanded_candidate & c.cluster_review.eq('retain_expanded_identity_unresolved')
 c['candidate_ambient_supported']=False
 if 'decontX_lineage_winner' in c:c['candidate_ambient_supported']=c.reviewed_expanded_candidate & c.decontX_lineage_winner.eq('Fibroblast') & c.decontX_detected_Fibroblast.ge(2)
 parts.append(c)
 allcounts.append({'sample':sm,'status':result['status'],'cells':len(c),'putative_patient_group':c.putative_patient_group.iloc[0],'author_CAF':int(c.author_CAF.sum()),'reference_pass':result['reference_pass'],'mechanical_candidates':int(c.mechanical_expanded_candidate.sum()),'reviewed_expanded_candidates':int(c.reviewed_expanded_candidate.sum()),'ambient_supported_candidates':int(c.candidate_ambient_supported.sum()),'manual_complete':not c.cluster_review.eq('pending').any()})
pd.DataFrame(allcounts).to_csv(P/'tables/ALL_LIBRARY_PROGRESS_COUNTS.tsv',sep='\t',index=False,encoding='utf-8')
if parts:
 c=pd.concat(parts,ignore_index=True);assert c.cell_id.is_unique
 c.to_csv(P/'tables/ALL_CELL_REVIEW_EVIDENCE.tsv.gz',sep='\t',index=False,encoding='utf-8')
 g=c.groupby('putative_patient_group').agg(libraries=('sample',lambda x:';'.join(sorted(set(x)))),cells=('cell_id','size'),mechanical_candidates=('mechanical_expanded_candidate','sum'),reviewed_expanded_candidates=('reviewed_expanded_candidate','sum'),ambient_supported_candidates=('candidate_ambient_supported','sum')).reset_index()
 g['patient_identity_verified']=False;g['grouping_basis']='sample suffix hypothesis; not verified independent biological N'
 g.to_csv(P/'tables/PUTATIVE_PATIENT_CANDIDATE_COUNTS.tsv',sep='\t',index=False,encoding='utf-8')
q={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'computational_evidence_summarized_manual_review_pending','computation_exit_code':run['exit_code'],'libraries_expected':19,'libraries_with_evidence':len(parts),'manual_review_pending_samples':sorted(missing_manual),'verified_independent_patient_N':None,'scientific_admission_complete':False,'no_strict_layer_gate':True}
(P/'REVIEW_PROGRESS.json').write_text(json.dumps(q,indent=2)+'\n',encoding='utf-8');print(json.dumps(q),flush=True)
