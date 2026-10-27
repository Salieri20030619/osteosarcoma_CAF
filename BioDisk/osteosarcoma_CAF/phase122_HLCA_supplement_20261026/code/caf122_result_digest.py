from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026')
for n in ['ANALYSIS_COMPLETE.json','SUPPLEMENT_QA.json','resources/MAPPING_EXTENDED_REVIEW.json']:
 if (O/n).exists():print(n,(O/n).read_text())
b=pd.read_csv(O/'tables/normal_IPF_background_comparison.tsv',sep='\t');programs=['mCAF_matrix','collagen_biosynthesis','collagen_crosslinking','collagen_processing_nonCOL','HLCA_IPF_reference','Osteogenic_customer']
for cond in ['Healthy','Healthy (tumor adjacent)','IPF']:
 for typ in ['Alveolar fibroblasts','Adventitial fibroblasts','Myofibroblasts']:
  g=b[(b.query_min_cells==5)&(b.reference_min_cells==20)&(b.condition==cond)&(b.reference_type==typ)&b.program.isin(programs)];print(cond,typ);print(g.pivot(index='node_id',columns='program',values='median_reference_percentile').round(3).to_string())
t=pd.read_csv(O/'tables/IPF_within_study_patient_tests.tsv',sep='\t');print('WITHIN_STUDY_IPF6');print(t[(t.reference_min_cells==20)&t.program.eq('HLCA_IPF_reference')&t.p.notna()].drop(columns='family').to_string(index=False))
i=pd.read_csv(O/'tables/identity_patient_threshold_summary.tsv',sep='\t');print('IDENTITY');print(i[(i.minimum_candidate_cells==5)&i.identity.isin(['Fibroblast','Smooth_muscle'])].to_string(index=False))
print('SENSITIVITY_BG');print(b[(b.query_min_cells==10)&(b.reference_min_cells==50)&(b.condition=='Healthy')&(b.reference_type=='Alveolar fibroblasts')&b.program.isin(programs)][['node_id','program','query_patients','median_reference_percentile']].to_string(index=False))
