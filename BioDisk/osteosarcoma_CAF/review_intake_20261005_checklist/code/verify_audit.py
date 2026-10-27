from pathlib import Path
import pandas as pd,json,hashlib
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'review_intake_20261005_checklist'
a=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False)
r=pd.read_csv(O/'tables/doublet_source_reconciliation.tsv.gz',sep='\t')
j=r.merge(a[['cell_id','doublet_positive_any']],on='cell_id',validate='one_to_one');known=j.doublet_positive_any.eq(True)
assert known.sum()==11 and j.loc[known,'source_positive'].all()
extra=j.source_positive&~known
print('EXTRA_FLAGS',int(extra.sum()),'previously_unassessed',int((extra&~j.exported_assessed).sum()),'previously_assessed',int((extra&j.exported_assessed).sum()))
j[j.source_positive].to_csv(O/'tables/model_positive_cells_for_review.tsv',sep='\t',index=False)
for d in json.loads((B/'phase65_spatial_structure_20261005/resources/INHERITED_DEPENDENCIES.json').read_text()):
 if 'phase63' in d['path']:assert hashlib.sha256(Path(d['path']).read_bytes()).hexdigest()==d['sha256']
q={'passed':True,'cells':len(r),'upstream_assessment_coverage':int(r.source_assessed.sum()),'downstream_assessment_records':int(r.exported_assessed.sum()),'upstream_model_positive':int(r.source_positive.sum()),'downstream_model_positive':int(known.sum()),'additional_positive_calls':int(extra.sum()),'additional_positive_previously_unassessed':int((extra&~j.exported_assessed).sum()),'counts_metadata_roster_match_pre_audit_phase65_hashes':True,'active_phase':json.loads((B/'CAF_NEXT_STAGE_CURRENT.json').read_text())['current_phase'],'active_pipeline_modified':False,'new_doublet_detection_performed':False,'negative_control':'Model-negative is not proof of singlet; no classifier gold-standard claim.'}
(O/'resources/AUDIT_ACCEPTANCE.json').write_text(json.dumps(q,indent=2));print(json.dumps(q),flush=True)
