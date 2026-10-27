import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pathlib,json,importlib.util,time
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase106_pooled_test_original_directions_20261020'
a=np.load(R/'phase104_full_gene_four_models_20261020/objects/tumor_common_patient_expression.npz')
Y=a['expression'];features=a['features'];site=a['sites'];src=a['sources']
names=pd.read_csv(R/'phase104_full_gene_four_models_20261020/tables/pooled_tumor_allgenes.tsv.gz',sep='\t').set_index('feature').gene.to_dict()
amend=dict(user_request='然后可以再加一个NCI和骨比作为第三，如果出来结果是0那就算了',third='NCI Lung17 vs original pooled Bone16',layer1='Pooled all39 P<0.05 + original17/23 concordant direction + third NCI contrast same direction',layer2='Layer1 plus third P<0.05',FDR='Each comparison full valid-gene BH, computed before any directional filter; no threshold tuning if zero',interpretation='Third contrast reuses original Bone patients and NCI patients from pooled test; supplementary consistency, not independent validation')
(P/'resources/PROTOCOL_THIRD_COMPARISON_AMENDMENT.json').write_text(json.dumps(amend,ensure_ascii=False,indent=2),encoding='utf-8')
corepath=R/'phase73_F134_bone_lung_20261009/strict_v2/code/analyze_bone_lung_v2.py'
spec=importlib.util.spec_from_file_location('frozen',corepath);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.SEED=20261012
b=Y[site=='Bone'];l=Y[src=='NCI'];assert len(b)==16 and len(l)==17
print('START third NCI17 vs Bone16',time.strftime('%H:%M:%S'),flush=True)
z=core.stats(b,l,False);z['feature']=features;z['gene']=[names[x] for x in features]
valid=z.P.notna();z['FDR']=np.nan;z.loc[valid,'FDR']=multipletests(z.loc[valid,'P'],method='fdr_bh')[1]
z.to_csv(P/'tables/third_NCI_vs_Bone_allgenes.tsv.gz',sep='\t',index=False,encoding='utf-8',na_rep='NA')
d=dict(status='computed',Bone_N=16,Lung_N=17,valid_tests=int(valid.sum()),nominal=int(z.P.lt(.05).sum()),FDR=int(z.FDR.lt(.05).sum()),test='MC65536_plus_observed_patient_label_permutation')
(P/'resources/THIRD_COMPARISON_COMPLETE.json').write_text(json.dumps(d,indent=2),encoding='utf-8');print(json.dumps(d),flush=True)
