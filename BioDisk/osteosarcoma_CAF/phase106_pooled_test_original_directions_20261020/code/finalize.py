import pathlib,json,hashlib,shutil
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase106_pooled_test_original_directions_20261020'
a=np.load(R/'phase104_full_gene_four_models_20261020/objects/tumor_common_patient_expression.npz');features=a['features'];Y=a['expression'];site=a['sites'];src=a['sources']
third=pd.read_csv(P/'tables/third_NCI_vs_Bone_allgenes.tsv.gz',sep='\t').set_index('feature').loc[features]
sign=lambda v:np.where(np.abs(v)<=1e-12,0,np.sign(v))
b=Y[site=='Bone'];l=Y[src=='NCI'];vp=b.var(0,ddof=1)/len(b)+l.var(0,ddof=1)/len(l)
wp=stats.ttest_ind(l,b,axis=0,equal_var=False).pvalue;wp[vp<=1e-24]=np.nan
wf=np.full(len(wp),np.nan);ok=np.isfinite(wp);wf[ok]=multipletests(wp[ok],method='fdr_bh')[1]
results=[];master=[]
for method in ['patient_permutation','Welch_sensitivity']:
 d=pd.read_csv(P/'tables'/(method+'_allgenes_with_directions.tsv.gz'),sep='\t').set_index('feature').loc[features].reset_index()
 d['third_effect']=third.mean_lung_minus_bone.to_numpy();d['third_P']=third.P.to_numpy() if method=='patient_permutation' else wp
 d['third_FDR']=third.FDR.to_numpy() if method=='patient_permutation' else wf
 d['third_direction_agrees']=(sign(d.third_effect.to_numpy())==sign(d.mean_lung_minus_bone.to_numpy()))&(sign(d.third_effect.to_numpy())!=0)
 d['select_two_directions']=d.selected_nominal
 d['select_three_directions']=d.select_two_directions&d.third_direction_agrees
 d['select_third_also_nominal']=d.select_three_directions&d.third_P.lt(.05)
 d['select_both_FDR']=d.select_third_also_nominal&d.FDR.lt(.05)&d.third_FDR.lt(.05)
 d.to_csv(P/'tables'/(method+'_three_direction_allgenes.tsv.gz'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 d[d.select_third_also_nominal].to_csv(P/'tables'/(method+'_third_P_selected_genes.tsv'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 d[d.original14].to_csv(P/'tables'/(method+'_original14_final.tsv'),sep='\t',encoding='utf-8',index=False,na_rep='NA')
 for layer in ['select_two_directions','select_three_directions','select_third_also_nominal','select_both_FDR']:
  m=d[layer]
  results.append(dict(method=method,layer=layer,genes=int(m.sum()),lung_higher=int((m&d.mean_lung_minus_bone.gt(0)).sum()),bone_higher=int((m&d.mean_lung_minus_bone.lt(0)).sum()),original14_genes=d.loc[m&d.original14,'gene'].tolist(),original14_count=int((m&d.original14).sum())))
 master.append(d)
(P/'resources/FINAL_RESULTS.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
pd.DataFrame([{**r,'original14_genes':';'.join(r['original14_genes'])} for r in results]).to_csv(P/'tables/final_filter_counts.tsv',sep='\t',encoding='utf-8',index=False)
# Independent verification of both MC calculations using equivalent direct group means.
qa=[];fi={v:i for i,v in enumerate(features)}
for name,dd,lmask in [('pooled',master[0],site=='Lung'),('third',third.reset_index(),src=='NCI')]:
 valid=dd[dd.P.notna()].sort_values('P')
 chosen=valid.iloc[np.unique(np.r_[np.arange(3),np.linspace(0,len(valid)-1,6,dtype=int)])].feature.tolist();ix=[fi[k] for k in chosen]
 v=np.vstack([Y[site=='Bone'][:,ix],Y[lmask][:,ix]]);nb=16;nl=int(lmask.sum());obs=abs(v[nb:].mean(0)-v[:nb].mean(0));hits=np.zeros(len(ix),int);total=v.sum(0)
 for j in range(65536):
  selected=np.random.default_rng(20261012+j).choice(nb+nl,nl,replace=False)
  sl=v[selected].sum(0);null=abs(sl/nl-(total-sl)/nb);hits+=null>=obs-1e-12
 pv=(hits+1)/65537;reported=dd.set_index('feature').loc[chosen,'P'].to_numpy();err=float(np.max(abs(pv-reported)));assert err<1e-12
 valid=dd[dd.P.notna()];errq=float(np.max(abs(multipletests(valid.P,method='fdr_bh')[1]-valid.FDR)));assert errq<1e-12
 qa.append(dict(comparison=name,genes=len(ix),max_P_error=err,max_BH_error=errq))
for d in master:
 assert ((d.select_third_also_nominal & ~d.select_three_directions).sum()==0)
 assert ((d.select_three_directions & ~d.select_two_directions).sum()==0)
 # Original directions must match direct source means.
 for source,col in [('SCPCP000017','source17_effect'),('SCPCP000023','source23_effect')]:
  e=Y[(src==source)&(site=='Lung')].mean(0)-Y[(src==source)&(site=='Bone')].mean(0)
  assert np.allclose(e,d[col].to_numpy(),atol=1e-12)
(P/'resources/FINAL_QA.json').write_text(json.dumps(dict(status='PASS',checks=qa,all_three_directions_verified=True,filter_nesting_verified=True),indent=2),encoding='utf-8')
print(json.dumps(results,indent=2));print('QA_PASS',qa)
