import pathlib,json
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase107_third_NCI_normal_background_20261020'
for d in ['code','tables','resources']:(P/d).mkdir(parents=True,exist_ok=True)
Q=R/'phase104_full_gene_four_models_20261020';B=R/'phase106_pooled_test_original_directions_20261020'
a=np.load(Q/'objects/tumor_common_patient_expression.npz');Y=a['expression'];features=a['features'];site=a['sites'];src=a['sources']
maptab=pd.read_csv(Q/'tables/gene_mapping_audit.tsv',sep='\t')
maps={}
for ds in ['normal_bone','normal_lung']:
 z=maptab[maptab.dataset.eq(ds)&maptab.status.eq('unique_match')];assert z.feature.is_unique
 maps[ds]=dict(zip(z.feature,z.row))
bm=pd.read_csv(Q/'tables/normal_bone_donor_totals.tsv',sep='\t');bc=pd.read_csv(Q/'tables/normal_bone_allgene_raw.tsv',sep='\t')
bone=np.log2(bc[bm.donor.tolist()].to_numpy().T/bm.library_UMIs.to_numpy()[:,None]*1e4+.1)
lung=np.load(Q/'objects/normal_lung_allgene.npz')
assert len(bm)==5 and len(lung['donors'])==3 and bm.cells.min()>=20 and lung['cells'].min()>=20
validmap=np.array([k in maps['normal_bone'] and k in maps['normal_lung'] for k in features])
nb=np.full((5,len(features)),np.nan);nl=np.full((3,len(features)),np.nan)
for j,k in enumerate(features):
 if validmap[j]:
  nb[:,j]=bone[:,maps['normal_bone'][k]];nl[:,j]=lung['expression'][:,maps['normal_lung'][k]]
protocol=dict(starting_sets='phase106 pooled significance + original17/23 and NCI direction + NCI P<0.05: permutation3956 primary, Welch3718 method-matched sensitivity',
correction='(tumor Lung - tumor Bone) - (normal Lung - normal Bone)',
comparisons=['all Lung23 vs Bone16','NCI Lung17 vs Bone16'],
normal='Bone5 donors1439 cells; Lung3 donors2136 cells; >=20 cells per normal donor',
inference='Four independent group Welch-Satterthwaite contrast, including normal sampling variance. Not reuse original permutation P.',
screen='Start from prior selected set; require corrected pooled and corrected NCI P<0.05 and both corrected signs agree with original two source tumor directions. Original source directions remain uncorrected. No original-source P gate.',
FDR='BH separately over all valid genes per corrected comparison BEFORE filtering prior3956; no FDR selection required for nominal endpoint',
missing='Normal coverage missing remains not evaluable, never zero or tested-negative',
same_method_control='Reconstruct raw two-group Welch before correction; compare to phase106 Welch baseline, plus apply raw Welch to prior permutation3956 to separate test-method effects',
strict_subset=False)
(P/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
names=pd.read_csv(Q/'tables/pooled_tumor_allgenes.tsv.gz',sep='\t').set_index('feature').gene.to_dict()
outputs=[];checks=[]
for label,lmask in [('pooled',site=='Lung'),('NCI',src=='NCI')]:
 b=Y[site=='Bone'];l=Y[lmask]
 raw_effect=l.mean(0)-b.mean(0)
 arrays=[l,b,nl,nb];signs=[1,-1,-1,1]
 corrected=raw_effect-(nl.mean(0)-nb.mean(0))
 v=np.vstack([v.var(0,ddof=1)/len(v) for v in arrays])
 vv=v.sum(0);df=vv**2/sum(v[i]**2/(len(arrays[i])-1) for i in range(4))
 rawv=v[0]+v[1];rawdf=rawv**2/(v[0]**2/(len(l)-1)+v[1]**2/(len(b)-1))
 valid=validmap&np.isfinite(df)&(vv>1e-24)&(rawv>1e-24)
 pv=np.full(len(features),np.nan);q=pv.copy();lo=pv.copy();hi=pv.copy()
 pv[valid]=2*stats.t.sf(abs(corrected[valid]/np.sqrt(vv[valid])),df[valid])
 q[valid]=multipletests(pv[valid],method='fdr_bh')[1]
 ci=stats.t.ppf(.975,df[valid])*np.sqrt(vv[valid]);lo[valid]=corrected[valid]-ci;hi[valid]=corrected[valid]+ci
 rawp=2*stats.t.sf(abs(raw_effect/np.sqrt(rawv)),rawdf);rawp[rawv<=1e-24]=np.nan
 data=pd.DataFrame(dict(feature=features,gene=[names[k] for k in features],comparison=label,normal_available=validmap,normal_effect=nl.mean(0)-nb.mean(0),raw_effect=raw_effect,raw_Welch_P=rawp,effect=corrected,SE=np.sqrt(vv),df=df,P=pv,FDR=q,CI_low=lo,CI_high=hi,valid_test=valid,family_valid_tests=int(valid.sum())))
 data.to_csv(P/'tables'/(label+'_normal_allgenes.tsv.gz'),sep='\t',index=False,encoding='utf-8',na_rep='NA');outputs.append(data)
 # Independent saturated group OLS + HC2 verification of effects and variance.
 inds=np.flatnonzero(valid)[np.linspace(0,int(valid.sum())-1,25,dtype=int)]
 for j in inds:
  yy=np.concatenate([x[:,j] for x in arrays]);X=np.zeros((len(yy),4));start=0
  for i,arr in enumerate(arrays):X[start:start+len(arr),i]=1;start+=len(arr)
  inv=np.linalg.inv(X.T@X);beta=inv@X.T@yy;h=np.einsum('ij,jk,ik->i',X,inv,X);res=yy-X@beta
  cov=inv@(X.T@((res*res/(1-h))[:,None]*X))@inv;c=np.array(signs)
  checks.extend([abs(c@beta-corrected[j]),abs(c@cov@c-vv[j])])
 if label=='pooled':
  old=pd.read_csv(Q/'tables/pooled_normal_min20_allgenes.tsv.gz',sep='\t').set_index('feature').loc[features]
  for col in ['effect','P','FDR']:
   assert np.allclose(data[col],old[col],equal_nan=True,atol=1e-10),col
print('CORRECTIONS_DONE',flush=True)
merged=outputs[0].merge(outputs[1],on='feature',suffixes=('_pooled','_NCI'))
summ=[];countrows=[]
sign=lambda v:np.where(np.abs(v)<=1e-12,0,np.sign(v))
for method in ['patient_permutation','Welch_sensitivity']:
 prev=pd.read_csv(B/'tables'/(method+'_three_direction_allgenes.tsv.gz'),sep='\t')
 d=prev.merge(merged,on='feature',validate='one_to_one')
 d['baseline_selected']=d.select_third_also_nominal
 s17=sign(d.source17_effect.to_numpy());s23=sign(d.source23_effect.to_numpy())
 def consistent(e):
  s=sign(e.to_numpy());return (s!=0)&(s==s17)&(s==s23)
 d['corrected_pooled_same_direction']=consistent(d.effect_pooled)
 d['corrected_NCI_same_direction']=consistent(d.effect_NCI)
 d['normal_evaluable']=d.valid_test_pooled&d.valid_test_NCI
 d['remain_after_pooled_correction']=d.baseline_selected&d.normal_evaluable&d.corrected_pooled_same_direction&d.P_pooled.lt(.05)
 d['remain_after_both_corrections']=d.remain_after_pooled_correction&d.corrected_NCI_same_direction&d.P_NCI.lt(.05)
 d['remain_both_full_FDR']=d.remain_after_both_corrections&d.FDR_pooled.lt(.05)&d.FDR_NCI.lt(.05)
 d['raw_Welch_method_control']=d.baseline_selected&d.raw_Welch_P_pooled.lt(.05)&d.raw_Welch_P_NCI.lt(.05)
 d.to_csv(P/'tables'/(method+'_allgenes_followup.tsv.gz'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 d[d.remain_after_both_corrections].to_csv(P/'tables'/(method+'_normal_retained_genes.tsv'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 d[d.original14].to_csv(P/'tables'/(method+'_original14_normal.tsv'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 for layer,mask in [('baseline',d.baseline_selected),('normal_evaluable',d.baseline_selected&d.normal_evaluable),('normal_missing',d.baseline_selected&~d.normal_evaluable),('raw_Welch_method_control',d.raw_Welch_method_control),('pooled_corrected',d.remain_after_pooled_correction),('both_corrected',d.remain_after_both_corrections),('both_corrected_FDR',d.remain_both_full_FDR)]:
  row=dict(method=method,layer=layer,genes=int(mask.sum()),original14=int((mask&d.original14).sum()),original14_genes=d.loc[mask&d.original14,'gene'].tolist())
  summ.append(row)
(P/'resources/RESULT_SUMMARY.json').write_text(json.dumps(summ,indent=2),encoding='utf-8')
pd.DataFrame([{**r,'original14_genes':';'.join(r['original14_genes'])} for r in summ]).to_csv(P/'tables/filter_summary.tsv',sep='\t',encoding='utf-8',index=False)
assert max(checks)<1e-10
qa=dict(status='PASS',independent_OLS_HC2_checks=len(checks)//2,max_error=float(max(checks)),phase104_pooled_normal_full_results_reproduced=True,normal_donors={'Bone':5,'Lung':3},normal_min_cells=20)
(P/'resources/QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(summ,indent=2));print(json.dumps(qa))
