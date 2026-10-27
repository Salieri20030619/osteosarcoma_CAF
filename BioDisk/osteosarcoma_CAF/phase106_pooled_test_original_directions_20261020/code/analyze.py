import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pathlib,json,hashlib,importlib.util,time
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase106_pooled_test_original_directions_20261020'
for d in ['code','tables','resources','logs']:(P/d).mkdir(parents=True,exist_ok=True)
a=np.load(R/'phase104_full_gene_four_models_20261020/objects/tumor_common_patient_expression.npz');Y=a['expression'];features=a['features'];site=a['sites'];src=a['sources']
base=pd.read_csv(R/'phase104_full_gene_four_models_20261020/tables/pooled_tumor_allgenes.tsv.gz',sep='\t').set_index('feature').loc[features]
original=pd.read_csv(R/'phase105_original_rule_NCI_20261020/tables/original_two_source_rule_reproduced.tsv.gz',sep='\t').set_index('feature').loc[features]
target=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text());target_ids={g['feature_ids'] for g in target}
protocol=dict(user_rule='Pool original two cohorts + NCI for statistical testing, and require direction consistency separately in original two cohorts without requiring original source P<0.05.',
primary_test='Original frozen patient-label permutation method:65536 seeded MC permutations + observed; Bone16 Lung23.',
seed=20261012,primary_selection='pooled P<0.05 AND signs(original17 effect, original23 effect, pooled effect) all nonzero and equal',
FDR='BH over every valid pooled gene test before directional filtering; report FDR-supported subset without requiring it for nominal selection',
direction='Frozen original independent patient rosters. Zero within numerical tolerance1e-12 is not a direction; no biological effect-size cutoff.',
gene_universe='20315 common measured tumor genes as phase104; gene-wise nonvariation => invalid test',
Welch='Apply same direction rule to existing phase104 pooled Welch results as separately labelled method sensitivity; do not mix P across methods.',
independence='Original direction cohorts are part of pooled test data; internal consistency not independent replication.',
normal_background='Not part of current rule',strict_subset='Not used',new_NCI_cell_count=1530,main_total_cells=6633)
(P/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
corepath=R/'phase73_F134_bone_lung_20261009/strict_v2/code/analyze_bone_lung_v2.py'
spec=importlib.util.spec_from_file_location('frozen',corepath);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.SEED=20261012
print('START pooled permutation',Y.shape,time.strftime('%H:%M:%S'),flush=True)
z=core.stats(Y[site=='Bone'],Y[site=='Lung'],False);z['feature']=features;z['gene']=base.gene.to_numpy()
valid=z.P.notna();z['FDR']=np.nan;z.loc[valid,'FDR']=multipletests(z.loc[valid,'P'],method='fdr_bh')[1]
sign=lambda x:np.where(np.abs(x)<=1e-12,0,np.sign(x))
def classify(d,method):
 d=d.copy();d['source17_effect']=original.mean_lung_minus_bone_17.to_numpy();d['source23_effect']=original.mean_lung_minus_bone_23.to_numpy()
 s=sign(d.mean_lung_minus_bone.to_numpy());s17=sign(d.source17_effect.to_numpy());s23=sign(d.source23_effect.to_numpy())
 d['direction_consistent']=(s!=0)&(s==s17)&(s==s23);d['selected_nominal']=d.direction_consistent&d.P.lt(.05);d['selected_FDR']=d.direction_consistent&d.FDR.lt(.05)
 d['original14']=d.feature.isin(target_ids);d['method']=method
 d.to_csv(P/'tables'/(method+'_allgenes_with_directions.tsv.gz'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 d[d.selected_nominal].to_csv(P/'tables'/(method+'_selected_genes.tsv'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 d[d.original14].to_csv(P/'tables'/(method+'_original14.tsv'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 return dict(method=method,valid_tests=int(d.P.notna().sum()),pooled_nominal=int(d.P.lt(.05).sum()),pooled_FDR=int(d.FDR.lt(.05).sum()),direction_consistent=int(d.direction_consistent.sum()),selected_nominal=int(d.selected_nominal.sum()),selected_FDR=int(d.selected_FDR.sum()),selected_lung_higher=int((d.selected_nominal&d.mean_lung_minus_bone.gt(0)).sum()),selected_bone_higher=int((d.selected_nominal&d.mean_lung_minus_bone.lt(0)).sum()),original14_retained=d.loc[d.selected_nominal&d.original14,'gene'].tolist(),original14_FDR=d.loc[d.selected_FDR&d.original14,'gene'].tolist())
summ=[classify(z,'patient_permutation')]
w=base.reset_index().rename(columns={'effect':'mean_lung_minus_bone'})
summ.append(classify(w,'Welch_sensitivity'))
(P/'resources/RESULT_SUMMARY.json').write_text(json.dumps(summ,indent=2),encoding='utf-8')
(P/'COMPUTATION_COMPLETE.json').write_text(json.dumps(dict(status='computed_QA_pending',results=summ),indent=2),encoding='utf-8')
print(json.dumps(summ,indent=2),flush=True)
