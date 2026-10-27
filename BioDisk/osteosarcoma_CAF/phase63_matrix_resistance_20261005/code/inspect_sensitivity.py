from pathlib import Path
import json,pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase63_matrix_resistance_20261005')
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
sets=json.loads((O/'resources/PROGRAMS.json').read_text())['sets'];print('CROSSLINK_GENES',sets['collagen_crosslinking'])
s=read(O/'tables/patient_program_statistics.tsv.gz');print('CROSSLINK_SENSITIVITY');print(s[(s.dataset=='SCPCP000017')&s.subgroup_id.isin(['F2_C1','F2_C3'])&s.program.isin(['collagen_crosslinking','LOX_family'])][['subgroup_id','matrix','scope','minimum_cells','program','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
g=read(O/'tables/targeted_mechanism_gene_statistics.tsv.gz');print('GENE_EVIDENCE');print(g[(g.dataset=='SCPCP000017')&(g.matrix=='corrected')&g.subgroup_id.isin(['F2_C1','F2_C3'])&g.gene.isin(sets['collagen_crosslinking'])][['subgroup_id','gene','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
f=read(O/'tables/spatial_patient_statistics.tsv');print('SPATIAL_SENSITIVITY');print(f[f.program.isin(['collagen_biosynthesis','collagen_assembly','collagen_crosslinking','LOX_family'])&f.endpoint.isin(['tumor_local_adjusted','tumor_neighbor_adjusted','CD8_local_adjusted'])][['version','program','endpoint','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
f=read(O/'tables/spatial_site_statistics.tsv');print('SITES',f.site_group.unique());print(f[f.program.eq('collagen_biosynthesis')&f.version.eq('full_program')&f.endpoint.isin(['tumor_local_adjusted','CD8_local_adjusted'])][['site_group','N','endpoint','mean_effect','positive_patients','P','FDR']].to_string(index=False))
