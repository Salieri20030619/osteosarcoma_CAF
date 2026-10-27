from pathlib import Path
import json,datetime
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase33_client_report_20260926';M=P/'phase30_two_layer_classification_20260926'
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
ev=read(M/'tables/expanded_classified_cells.tsv.gz').set_index('cell_id');su=read(M/'tables/supported_classified_cells.tsv.gz').set_index('cell_id')
summary=[];details=[];qa=[]
for sm in ['SCPCS000614','SCPCS000610']:
 paths={v:O/'cnv_reference'/sm/v for v in ['source_reference','lymphoid_reference']}
 accepted=json.loads((O/'delivery/REFERENCE_OUTPUT_INDEPENDENT_ACCEPTANCE.json').read_text()); valid={v:any(r['sample']==sm and r['variant']==v and r['independent_call_reconstruction_exact'] and r['source_cell_order_exact'] for r in accepted['records']) for v in paths}
 qa.append(dict(sample=sm,**{v+'_complete':ok for v,ok in valid.items()}))
 if not all(valid.values()):continue
 a=read(paths['source_reference']/'cell_results.tsv').set_index('cell_id');b=read(paths['lymphoid_reference']/'cell_results.tsv').set_index('cell_id');refs=read(O/'inputs'/sm/'reference_membership.tsv.gz').set_index('cell_id')
 assert a.index.equals(b.index) and a.index.equals(refs.index) and a.index.is_unique
 tab=pd.DataFrame({'source_reference_call':a['call'],'lymphoid_reference_call':b['call'],'source_reference_member':refs.reference,'lymphoid_reference_member':refs.lymphoid_reference,'provider_annotation':refs.provider_annotation})
 tab['expanded_class']=ev.class_id.reindex(tab.index);tab['supported_class']=su.class_id.reindex(tab.index)
 tab['original_copykat']=ev.copykat_call.reindex(tab.index);tab['original_provider_or_native_infercnv']=ev.infercnv_call.reindex(tab.index)
 tab['sample']=sm
 for layer,labels in [('expanded',ev),('supported',su)]:
  ids=labels.index[labels['sample'].eq(sm)]
  for target in ['all_candidates','F1','F2','F3','F4','F5']:
   subids=ids if target=='all_candidates' else labels.loc[ids].index[labels.loc[ids,'class_id'].eq(target)]
   s=tab.loc[subids];joint=s.source_reference_call.ne('not_returned_by_software')&s.lymphoid_reference_call.ne('not_returned_by_software')
   z=s.loc[joint]
   summary.append(dict(sample=sm,layer=layer,target=target,cells=len(s),baseline_reference_members=int(s.source_reference_member.sum()),alternate_reference_members=int(s.lymphoid_reference_member.sum()),jointly_evaluable=len(z),same_call=int(z.source_reference_call.eq(z.lymphoid_reference_call).sum()),source_no_alteration=int(z.source_reference_call.eq('no_alteration_called').sum()),alternate_no_alteration=int(z.lymphoid_reference_call.eq('no_alteration_called').sum()),changed_to_alteration=int((z.source_reference_call.eq('no_alteration_called')&z.lymphoid_reference_call.eq('alteration_called')).sum()),changed_to_no_alteration=int((z.source_reference_call.eq('alteration_called')&z.lymphoid_reference_call.eq('no_alteration_called')).sum()),source_not_independently_called=int(s.source_reference_call.eq('not_returned_by_software').sum()),alternate_not_returned=int(s.lymphoid_reference_call.eq('not_returned_by_software').sum())))
  # All cells and reference-control behavior remain available for review.
 details.append(tab.reset_index(names='cell_id'))
if details:pd.concat(details,ignore_index=True).to_csv(O/'tables/reference_sensitivity_all_cell_comparison_accepted.tsv.gz',sep='\t',index=False)
pd.DataFrame(summary).to_csv(O/'tables/reference_sensitivity_candidate_summary_accepted.tsv',sep='\t',index=False)
out={'at':datetime.datetime.now().astimezone().isoformat(),'samples':qa,'complete_comparisons':sum(all(v for k,v in r.items() if k.endswith('_complete')) for r in qa),'primary_labels_unchanged':True,'no_hard_CNV_signal_threshold_added':True,'reference_members_not_automatically_counted_as_independent_no_CNV_predictions':True}
(O/'delivery/REFERENCE_COMPARISON_ACCEPTED_QA.json').write_text(json.dumps(out,indent=2))
(O/'REFERENCE_COMPARISON_ACCEPTED_FINISHED').write_text(out['at']);print(json.dumps(out));print(pd.DataFrame(summary).to_string(index=False))
