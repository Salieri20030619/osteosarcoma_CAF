from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005'
a=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False)
r=pd.read_csv(O/'resources/doublet_source_reconciliation.tsv.gz',sep='\t').set_index('cell_id')
assert len(a)==8174 and a.cell_id.is_unique and r.index.is_unique and set(a.cell_id)==set(r.index)
r=r.loc[a.cell_id];assert r.source_assessed.all() and r.source_positive.sum()==26
for col in ['doublet_assessed_any','doublet_positive_any']:a['legacy_'+col]=a[col]
a['doublet_assessed_any']=r.source_assessed.to_numpy(bool);a['doublet_positive_any']=r.source_positive.to_numpy(bool)
a['doublet_evidence_methods']=r.source_methods.to_numpy()
a['doublet_evidence_source']='phase26 per-sample raw obs plus successful phase31 scDblFinder; reconciled and verified'
a['doublet_evidence_version']=O.name
a['doublet_review_status']=a.doublet_positive_any.map({True:'model_positive_requires_review',False:'assessed_no_model_positive'})
a['include_exclude_doublet_sensitivity']=~a.doublet_positive_any
a.to_csv(O/'tables/cells_repaired_with_confidence.tsv.gz',sep='\t',index=False)
a[['cell_id','doublet_assessed_any','doublet_positive_any','doublet_evidence_methods','doublet_evidence_source','doublet_evidence_version','doublet_review_status','include_exclude_doublet_sensitivity']].to_csv(O/'tables/current_doublet_evidence.tsv.gz',sep='\t',index=False)
changes=a[a.legacy_doublet_assessed_any.isna() | a.legacy_doublet_positive_any.isna()]
changes[['cell_id','class_id','subgroup_id','legacy_doublet_assessed_any','legacy_doublet_positive_any','doublet_assessed_any','doublet_positive_any','doublet_evidence_methods']].to_csv(O/'tables/repair_change_log.tsv.gz',sep='\t',index=False)
summary=a.groupby(['class_id','subgroup_id']).agg(cells=('cell_id','size'),assessed=('doublet_assessed_any','sum'),model_positive=('doublet_positive_any','sum'),retained_sensitivity=('include_exclude_doublet_sensitivity','sum')).reset_index()
summary.to_csv(O/'tables/subgroup_doublet_summary.tsv',sep='\t',index=False)
(O/'resources/STATISTICAL_IMPLEMENTATION_CLARIFICATION.json').write_text(json.dumps({'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rule':'Use same phase58 common.tests as inherited phase63: exact for N<=15, MonteCarlo19999+observed otherwise. This implementation clarification precedes any phase66 sensitivity results and preserves comparability. No single-cell-as-N tests.'},indent=2))
(O/'REPAIR_COMPLETE.json').write_text(json.dumps({'cells':8174,'assessed':8174,'positive':26,'repaired_missing_records':len(changes),'sensitivity_retained':8148,'original_assignments_and_confidence_preserved':True,'legacy_files_not_overwritten':True},indent=2))
shutil.copy2(B/'tmp_phase66_freeze.py',O/'code/freeze.py')
print('REPAIRED',len(changes),summary[summary.subgroup_id.isin(['F2_C1','F2_C3'])].to_string(index=False))
