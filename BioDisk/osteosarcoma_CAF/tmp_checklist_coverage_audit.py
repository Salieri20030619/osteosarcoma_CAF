from pathlib import Path
import pandas as pd,numpy as np,json
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'review_intake_20261005_checklist';O.mkdir(exist_ok=True)
for d in ['sources','tables','resources','delivery','code']:(O/d).mkdir(exist_ok=True)
a=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False)
z=np.load(B/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz');assert list(z['cells'])==list(a.cell_id)
a['mechanism_count_changed']=np.any(np.abs(z['raw']-z['corrected'])>1e-8,axis=1)
a['library_changed']=abs(z['raw_library']-z['corrected_library'])>1e-8
a['doublet_record_available']=a.doublet_assessed_any.notna()&a.doublet_assessed_any.eq(True)
a['ambient_fraction_available']=a.decontx_contamination.notna()
rows=[]
for typ in ['class_id','subgroup_id','dataset']:
 for k,t in a.groupby(typ):
  rows.append({'group_type':typ,'group':k,'cells':len(t),'doublet_assessed':int(t.doublet_record_available.sum()),'doublet_unassessed':int((~t.doublet_record_available).sum()),'doublet_positive_among_assessed':int(t.doublet_positive_any.eq(True).sum()),'ambient_fraction_available':int(t.ambient_fraction_available.sum()),'mechanism_counts_changed':int(t.mechanism_count_changed.sum()),'library_changed':int(t.library_changed.sum()),'two_CNV_methods_no_alteration':int(t.CNV_evidence_current.eq('concordant_no_alteration').sum())})
d=pd.DataFrame(rows);d.to_csv(O/'tables/existing_QC_coverage_audit.tsv',sep='\t',index=False)
print(d[(d.group_type=='class_id')|d['group'].isin(['F2_C1','F2_C3','F4_C8','F4_C10'])].to_string(index=False))
for rel in ['phase54_raw_review_20261004/RAW_REVIEW_QA.json','phase54_raw_review_20261004/CNV_REVIEW_QA.json','phase54_raw_review_20261004/FINAL_DATA_QA.json','phase32_lung_sensitivity_20260926/FINAL_STATUS.json','phase32_lung_sensitivity_20260926/resources/SENSITIVITY_PROTOCOL.json']:
 p=B/rel;d=json.loads(p.read_text());print('JSON',rel,json.dumps(d,ensure_ascii=False)[:5500])
for sub in ['phase54_raw_review_20261004/tables','phase54_raw_review_20261004/code','phase41_full_transcriptome_setup_20261001/resources']:
 p=B/sub;print(sub,[x.name for x in p.iterdir()][:45])
