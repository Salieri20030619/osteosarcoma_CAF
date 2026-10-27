from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase66_doublet_repair_DSC1_20261005')
d=pd.read_csv(O/'tables/cells_repaired_with_confidence.tsv.gz',sep='\t',low_memory=False)
print('POSITIVES',d[d.doublet_positive_any].groupby(['dataset','subgroup_id']).size().to_string())
w=pd.read_csv(O/'tables/program_before_after_comparison.tsv',sep='\t')
w['changed']=w.effect_change.abs()>1e-10
print('COMPARE',len(w),int(w.changed.sum()),'nominalchanges',int(w.nominal_change.sum()),'lost',int(w.lost_evaluability.sum()),'flips',int(w.direction_flip.sum()))
print(w[w.nominal_change | w.lost_evaluability][['dataset','subgroup_id','matrix','scope','minimum_cells','program','version','N_all','N_excluded','mean_effect_all','mean_effect_excluded','P_all','P_excluded']].to_string(index=False)[:12000])
print('NCHANGE',w[w.N_all!=w.N_excluded][['dataset','subgroup_id','matrix','scope','minimum_cells','N_all','N_excluded']].drop_duplicates().to_string(index=False))
m=pd.read_csv(O/'tables/desmosome_cell_counts.tsv.gz',sep='\t',low_memory=False,float_precision='round_trip')
print('DSC1POS',m[m.DSC1>0][['cell_id','dataset','subgroup_id','matrix','DSC1','doublet_positive_any']].to_string(index=False))
