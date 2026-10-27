from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase65_spatial_structure_20261005')
s=pd.read_csv(O/'tables/discovery_structure_statistics.tsv',sep='\t');r=pd.read_csv(O/'tables/discovery_patient_leave_one_out.tsv',sep='\t')
mask=s.dataset.eq('SCPCP000017')&s.matrix.eq('corrected')&s.scope.eq('all')&s.minimum_cells.eq(1)&s.subgroup_id.isin(['F2_C3','F2_C1','F3_C1','F4_C8','F1_originalC10'])
cols=['subgroup_id','program','version','N','mean_effect','positive_patients','P','FDR']
print(s.loc[mask,cols].to_string(index=False))
print('LOPO',r[r.dataset.eq('SCPCP000017')&r.matrix.eq('corrected')&r.scope.eq('all')&r.minimum_cells.eq(1)&r.subgroup_id.isin(['F2_C1','F2_C3'])][['subgroup_id','program','version','N','loo_min','loo_max','all_loo_same_direction']].to_string(index=False))
