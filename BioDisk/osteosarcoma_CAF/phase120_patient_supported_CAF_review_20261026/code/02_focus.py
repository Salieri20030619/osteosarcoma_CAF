from pathlib import Path
import pandas as pd,json
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase120_patient_supported_CAF_review_20261026';I=R/'phase119_Bone_Lung_functional_subgroups_20261024'
d=pd.read_csv(O/'tables/functional_audit.tsv',sep='\t');sp=pd.read_csv(O/'tables/node_donor_support.tsv',sep='\t');m=pd.read_csv(O/'tables/cell_identity_audit.tsv.gz',sep='\t');l=pd.read_csv(I/'tables/hierarchical_cell_labels.tsv.gz',sep='\t');g=pd.read_csv(O/'tables/source_candidate_gate.tsv',sep='\t')
focus=['Bone_07','Bone_07_s1','Lung_04','Lung_04_s2'];ends=['program__mCAF_matrix','program__collagen_crosslinking','program__fibrillar_I_III_V','audit__collagen_processing_nonCOL','audit__FNDC1_paper','audit__MYH11_paper','audit__CD36_paper','audit__C7_paper','audit__fibro_identity','audit__immune_identity','audit__osteogenic_identity']
z=d[d.node_id.isin(focus)&d.endpoint.isin(ends)&d['mode'].isin(['main_min5','identity_min5'])&d.N.gt(0)];z.to_csv(O/'tables/focus_functional_audit.tsv',sep='\t',index=False)
print(z[['node_id','source','mode','endpoint','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
print('LUNG04_DONORS',sp[sp.node_id.eq('Lung_04')].to_string(index=False))
for node in ['Lung_04','Bone_07']:
 a=m[~m.doublet&m.site.eq('Lung' if node.startswith('Lung') else 'Bone')].copy();a['target']=l.loc[a.index,'parent'].eq(node);c=a.groupby(['source','donor','target']).agg(cells=('cell_id','size'),identity=('identity_sensitivity','sum'),immune=('immune_coexpression','sum'),endo=('endothelial_coexpression','sum'));c.to_csv(O/'tables'/(node+'_identity_attrition.tsv'),sep='\t');print(node,c.to_string())
print('LITERATURE_POSITIVE',d[d['mode'].eq('main_min5')&d.endpoint.str.startswith('audit__')&d.N.ge(6)&d.mean_effect.gt(0)&d.P.lt(.05)][['node_id','source','endpoint','N','mean_effect','positive_patients','P','FDR']].to_string(index=False))
print('COUNTS',m[~m.doublet].groupby(['source','site']).agg(cells=('cell_id','size'),donors=('donor','nunique'),samples=('sample','nunique')).to_string())
