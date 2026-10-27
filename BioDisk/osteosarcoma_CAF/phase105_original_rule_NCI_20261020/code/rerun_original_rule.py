import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pathlib,json,hashlib,importlib.util,time
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase105_original_rule_NCI_20261020'
for d in ['code','resources','tables','logs']:(P/d).mkdir(parents=True,exist_ok=True)
corepath=R/'phase73_F134_bone_lung_20261009/strict_v2/code/analyze_bone_lung_v2.py'
spec=importlib.util.spec_from_file_location('frozen_permutation',corepath);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.SEED=20261012
a=np.load(R/'phase104_full_gene_four_models_20261020/objects/tumor_common_patient_expression.npz')
Y=a['expression'];features=a['features'];src=a['sources'];site=a['sites']
names=pd.read_csv(R/'phase104_full_gene_four_models_20261020/tables/pooled_tumor_allgenes.tsv.gz',sep='\t').set_index('feature').gene.to_dict()
gene=[names[x] for x in features];target=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text());target_ids={g['feature_ids'] for g in target}
old=pd.read_csv(R/'phase81_allCAF_bone_lung_20261012/tables/all_genes_statistics.tsv.gz',sep='\t')
old=old[old.layer.eq('expanded')&old.sensitivity.eq('all_cells')&old.design.eq('independent')&old.dataset.isin(['SCPCP000017','SCPCP000023'])]
cols=['feature','symbol','mean_lung_minus_bone','P','FDR']
o17=old[old.dataset.eq('SCPCP000017')][cols];o23=old[old.dataset.eq('SCPCP000023')][cols]
base=o17.merge(o23,on='feature',suffixes=('_17','_23'))
mask=(base.P_17<.05)&(base.P_23<.05)&(base.mean_lung_minus_bone_17*base.mean_lung_minus_bone_23>0)
assert set(base.loc[mask,'feature'])==target_ids
base['selected']=mask;base['gene']=base.feature.map(names).fillna(base.symbol_17)
base.to_csv(P/'tables/original_two_source_rule_reproduced.tsv.gz',sep='\t',index=False,encoding='utf-8')
protocol=dict(original_rule='Separate unpaired source-specific patient label permutation tests; intersect raw P<0.05 and concordant nonzero signs; no FDR or effect-size gate.',
original_reproduced=14,original_seed=20261012,exact_policy='All assignments if <=65536; otherwise65536 fixed MC permutations plus observed correction; frozen original code and 1e-12 comparison tolerance.',
expanded17='SCPCP000017 Bone12 vs its Lung3 + NCI17 = Lung20',
expanded23='SCPCP000023 Bone4 vs its Lung3 + NCI17 = Lung20',
NCI_only17='SCPCP000017 Bone12 vs NCI Lung17',
NCI_only23='SCPCP000023 Bone4 vs NCI Lung17',
gene_scope='20315 complete common eligible genes from phase104; all14 original candidates present; preserve all-gene P/FDR.',
FDR='BH per comparison across valid genes; intersection uses raw P not FDR.',
shared_controls='Augmented comparisons share17 NCI patients; not independent cohorts. NCI-only comparisons share same17 lungs, and source/assay/annotation differences remain.',
normal_background='Not applied: original14 selection was tumor-site contrasts without normal subtraction.',
strict_subsets='Not used, following user instruction',raw_counts_scope='Original expanded candidate set + all NCI author CAF')
(P/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
(P/'resources/INPUT_HASHES.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [corepath,R/'phase104_full_gene_four_models_20261020/objects/tumor_common_patient_expression.npz']],indent=2),encoding='utf-8')
results={};counts=[]
for strategy in ['augmented','NCI_only']:
 for source in ['SCPCP000017','SCPCP000023']:
  b=Y[(src==source)&(site=='Bone')]
  lungmask=(src=='NCI')|((src==source)&(site=='Lung')) if strategy=='augmented' else src=='NCI'
  l=Y[lungmask]
  print('START',strategy,source,b.shape,l.shape,time.strftime('%H:%M:%S'),flush=True)
  z=core.stats(b,l,False);z['feature']=features;z['gene']=gene;z['FDR']=np.nan
  valid=z.P.notna();z.loc[valid,'FDR']=multipletests(z.loc[valid,'P'],method='fdr_bh')[1]
  z['original14']=z.feature.isin(target_ids);z['strategy']=strategy;z['source']=source
  z.to_csv(P/'tables'/(strategy+'_'+source+'_allgenes.tsv.gz'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
  results[strategy,source]=z
  counts.append(dict(strategy=strategy,source=source,Bone_N=len(b),Lung_N=len(l),valid_tests=int(valid.sum()),nominal=int(z.P.lt(.05).sum()),FDR=int(z.FDR.lt(.05).sum()),test=';'.join(z.loc[valid,'test'].unique())))
  print('DONE',counts[-1],flush=True)
summary=[]
for strategy in ['augmented','NCI_only']:
 z=results[strategy,'SCPCP000017'].merge(results[strategy,'SCPCP000023'],on='feature',suffixes=('_17','_23'))
 z['gene']=z.feature.map(names)
 z['selected']=(z.P_17<.05)&(z.P_23<.05)&(z.mean_lung_minus_bone_17*z.mean_lung_minus_bone_23>0)
 z['same_as_original_direction']=np.where(z.gene.eq('EBF1'),(z.mean_lung_minus_bone_17<0)&(z.mean_lung_minus_bone_23<0),(z.mean_lung_minus_bone_17>0)&(z.mean_lung_minus_bone_23>0))
 z['original14']=z.feature.isin(target_ids)
 z['both_FDR']=z.FDR_17.lt(.05)&z.FDR_23.lt(.05)&z.selected
 z.to_csv(P/'tables'/(strategy+'_intersection_allgenes.tsv.gz'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
 z[z.selected].to_csv(P/'tables'/(strategy+'_selected_genes.tsv'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
 z[z.original14].to_csv(P/'tables'/(strategy+'_original14_followup.tsv'),sep='\t',index=False,encoding='utf-8',na_rep='NA')
 keep=z.selected&z.original14&z.same_as_original_direction
 d=dict(strategy=strategy,intersection=int(z.selected.sum()),lung_higher=int((z.selected&z.mean_lung_minus_bone_17.gt(0)).sum()),bone_higher=int((z.selected&z.mean_lung_minus_bone_17.lt(0)).sum()),both_FDR=int(z.both_FDR.sum()),original14_selected=int((z.selected&z.original14).sum()),original14_same_direction_retained=int(keep.sum()),retained_original14=z.loc[keep,'gene'].tolist())
 summary.append(d);print('INTERSECTION',json.dumps(d),flush=True)
pd.DataFrame(counts).to_csv(P/'tables/comparison_summary.tsv',sep='\t',index=False,encoding='utf-8')
(P/'resources/RESULT_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(P/'COMPUTATION_COMPLETE.json').write_text(json.dumps(dict(status='computed_QA_pending',original_rule_reproduced=14,comparisons=counts,intersections=summary),indent=2),encoding='utf-8')
