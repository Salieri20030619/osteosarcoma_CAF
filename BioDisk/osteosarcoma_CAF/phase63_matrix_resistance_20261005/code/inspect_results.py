from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase63_matrix_resistance_20261005')
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
s=read(O/'tables/patient_program_statistics.tsv.gz');base=read(O/'tables/subgroup_baseline.tsv')
main=s[(s.matrix=='corrected')&(s.scope=='all')&(s.minimum_cells==1)]
sel=['collagen_biosynthesis','collagen_assembly','collagen_crosslinking','collagen_degradation','LOX_family']
for r in base.itertuples():
 z=main[(main.subgroup_id==r.subgroup_id)&(main.dataset==r.primary_source)&main.program.isin(sel)]
 print(r.subgroup_id,r.primary_source,' | '.join(f'{t.program}:N{t.N},d{t.mean_effect:.2f},P{t.P:.4g},q{t.FDR:.4g}' for t in z.itertuples()))
print('ALL_NOMINAL_POSITIVE')
z=main[main.program.isin(sel)&main.P.lt(.05)&main.mean_effect.gt(0)].sort_values(['program','P'])
print(z[['dataset','subgroup_id','program','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
if (O/'SPATIAL_COMPLETE').exists():
 f=read(O/'tables/spatial_patient_statistics.tsv');print('SPATIAL');print(f[(f.version=='full_program')&f.program.isin(sel)][['program','endpoint','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
