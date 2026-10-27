from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase86_functional_bone_lung_20261012')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t')
f=read('all_functional_bone_lung_statistics.tsv');a=read('all_marker_functional_associations.tsv.gz');p=read('patient_site_programs.tsv');r=read('ALL14_retention_registry.tsv')
print('MINIMUM_FUNCTIONAL_P',f.P.min(),'valid',f.P.notna().sum())
print('EXPANDED_PRIMARY',f[(f.scope=='expanded_candidates')&(f.matrix=='raw')&(f.minimum_cells==1)&(f.design=='independent')&(f.endpoint=='mean_score')][['dataset','program','bone_N','lung_N','mean_lung_minus_bone','P']].to_string(index=False))
print('MARKER_SUMMARY',a.groupby(['scope','dataset']).agg(rows=('gene','size'),valid=('P','count'),nominal=('nominal','sum')).to_string())
names=['CD248','FAP','LRRC15','PDPN','THY1','PDGFRA','PDGFRB','MCAM','ANTXR1','CD36','GPC1','GPC3','SDC2','SDC3']
t=a[(a.scope=='Fibro_reviewed')&(a.dataset=='SCPCP000017')&a.gene.isin(names)&a.program.isin(['collagen_biosynthesis','collagen_assembly','collagen_crosslinking','ECM_organization','inflammation','TGF_beta'])]
print('REFERENCE_PROXIES\n'+t[['program','gene','N','mean_high_minus_low','P','FDR','component_of_score']].sort_values(['program','P']).to_string(index=False))
top=a[(a.scope=='Fibro_reviewed')&a.surface_annotation&a.P.lt(.05)&a.mean_high_minus_low.gt(0)].sort_values(['program','P']).groupby('program',group_keys=False).head(5)
print('TOP_FUNCTION_MARKERS\n'+top[['dataset','program','gene','N','mean_high_minus_low','P','FDR']].to_string(index=False))
print('ALL14\n'+r.to_string(index=False))
