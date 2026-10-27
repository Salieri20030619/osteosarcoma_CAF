from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_identity_20260919');T=O/'tables'
assert (O/'logs/identity.exit_code').read_text().strip()=='0'
a=pd.read_csv(T/'all_cell_identity_evidence.tsv.gz',sep='\t',dtype={'cluster':str})
original=a[['cell_id','raw_strict_candidate','caf_supported_candidate','caf_supported_positive_sensitivity']].copy()
a['count_threshold_sensitive_candidate']=a.raw_strict_candidate&~a.caf_supported_candidate&a.caf_supported_positive_sensitivity
a['trace_competing_sensitive_primary']=a.caf_supported_candidate&~a.caf_supported_positive_sensitivity
a['supported_under_both_detection_rules']=a.caf_supported_candidate&a.caf_supported_positive_sensitivity
a.loc[a.count_threshold_sensitive_candidate,'identity_decision']='corrected_detection_threshold_sensitive'
a.loc[a.raw_strict_candidate&~a.caf_supported_candidate&~a.caf_supported_positive_sensitivity,'identity_decision']='corrected_lineage_evidence_unresolved'
assert a[original.columns].equals(original)
a.to_csv(T/'all_cell_identity_evidence.tsv.gz',sep='\t',index=False)
for sm,g in a.groupby('sample',sort=False):g.to_csv(T/(sm+'_identity.tsv.gz'),sep='\t',index=False)
cols=['cell_id','sample','patient_id','dataset','site','technology','identity_decision','raw_det_Fib','corrected_det_Fib','corrected_pos_Fib','raw_pct_Fib','corrected_pct_Fib','decontx_contamination','exclusion_reasons']
a.loc[a.count_threshold_sensitive_candidate,cols].to_csv(T/'sensitivity_only_CAF_candidates.tsv',sep='\t',index=False)
a.loc[a.trace_competing_sensitive_primary,cols].to_csv(T/'primary_trace_competitor_sensitivity.tsv',sep='\t',index=False)
a[a.lineage_phase3=='FibLike'].groupby(['sample','identity_decision']).size().rename('cells').reset_index().rename(columns={'identity_decision':'decision'}).to_csv(T/'identity_transitions.tsv',sep='\t',index=False)
s=a.groupby('sample')[['caf_supported_candidate','count_threshold_sensitive_candidate','trace_competing_sensitive_primary','supported_under_both_detection_rules']].sum().reset_index()
s.to_csv(T/'candidate_sensitivity_tiers.tsv',sep='\t',index=False)
for sm,g in a[a.raw_strict_candidate].groupby('sample'):
 z=g.groupby(['caf_supported_candidate','caf_supported_positive_sensitivity']).agg(cells=('cell_id','size'),n_genes_median=('n_genes','median'),contamination_median=('decontx_contamination','median'),raw_Fib_pct_median=('raw_pct_Fib','median'),corrected_Fib_pct_median=('corrected_pct_Fib','median')).reset_index()
 z.insert(0,'sample',sm)
 z.to_csv(T/(sm+'_sensitivity_evidence.tsv'),sep='\t',index=False)
state=json.loads((O/'STATUS.json').read_text())
state.update(sensitivity_only_candidates=int(a.count_threshold_sensitive_candidate.sum()),primary_trace_competitor_sensitive=int(a.trace_competing_sensitive_primary.sum()),supported_under_both_detection_rules=int(a.supported_under_both_detection_rules.sum()),sensitivity_annotation='Distinguish fractional-count detection sensitivity from loss of fibroblast expression; primary candidate flags unchanged')
(O/'STATUS.json').write_text(json.dumps(state,indent=2))
print(s.to_string(index=False))
