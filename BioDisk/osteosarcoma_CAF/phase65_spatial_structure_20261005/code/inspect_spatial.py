from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase65_spatial_structure_20261005')
d=pd.read_csv(O/'tables/spatial_gene_and_program_scores.tsv.gz',sep='\t',low_memory=False);print('ZERO',d[d.library_UMI<=0][['sample','spot_id']].to_string(index=False));d[d.library_UMI<=0][['sample','spot_id','patient_id']].to_csv(O/'QA/zero_library_spots.tsv',sep='\t',index=False)
s=pd.read_csv(O/'tables/spatial_patient_statistics.tsv',sep='\t')
cols=['cohort_scope','version','program','context','model','N','mean_effect','positive_patients','P','FDR']
print('MAIN',s[(s.version=='full')&(s.model=='depth_xy')][cols].to_string(index=False))
print('MUTUAL',s[(s.version=='full')&(s.model=='mutual_contexts')][cols].to_string(index=False))
print('DROP',s[(s.version=='drop_dominant_pair')&(s.model=='depth_xy')][cols].to_string(index=False))
c=pd.read_csv(O/'tables/source_preference_statistics.tsv',sep='\t')
print('CONTRAST',c[(c.version=='full')&(c.model=='depth_xy')][['cohort_scope','program','contrast','N','mean_effect','P','FDR']].to_string(index=False))
p=O/'code/verify.py';v=p.read_text().replace("sep='\\t')","sep='\\t',low_memory=False)")
v=v.replace("invalid=s[s.library_UMI<=0];assert len(invalid)==1 and invalid.spot_id.iloc[0]=='MR4__GGAGGGTCAAGTAAGA-1'","invalid=s[s.library_UMI<=0]\nbase=pd.read_csv(B/'phase46_spatial_functional_20261001/tables/all_spot_RNA_and_components.tsv.gz',sep='\\t',low_memory=False)\nassert set(invalid.spot_id)==set(base.loc[base.library_UMI<=0,'spot_id'])")
p.write_text(v)
