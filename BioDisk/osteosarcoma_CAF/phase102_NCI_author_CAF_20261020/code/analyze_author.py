import pathlib,json,math
import numpy as np,pandas as pd
from scipy import stats
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase102_NCI_author_CAF_20261020'
nci=pd.read_csv(P/'tables/author_CAF_patient_expression.tsv',sep='\t')
old=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/original_tumor_patient_expression.tsv',sep='\t')
bone=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/bone_reference_with_coverage_sensitivity.tsv',sep='\t')
lung=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/lung_reference_pseudobulk.tsv',sep='\t')
bone=bone[bone.specification.eq('min20')];lung=lung[lung.specification.eq('min20')]
registry=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text())
protocol=dict(main='Original source-specific exact permutation results preserved. NCI author CAF vs each original source Bone: cross-source exploratory comparison.',scale='patient pseudobulk log2(CP10k+0.1); sum raw counts and totals across paired libraries',strata=['author_all','author_singlets'],normal='5 bone reference donors and 3 lung reference donors, each >=20 cells; four-group Welch contrast (NCI Lung - original Bone) - (normal Lung - normal Bone)',family='BH separately per source x stratum x contrast, 14 prespecified genes, 12 normal-evaluable',combined='Common site effect controlling source using within-source contrast weights nBone*nLung/n; heteroscedastic variance and Welch-Satterthwaite df. NCI only Lung provides zero within-source site information.',boundaries='NCI vs original Bone confounds cohort, assay and annotation; cannot independently validate site effect.')
(P/'resources/PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
def contrast(arrays,signs):
 if any(len(x)<2 or not np.isfinite(x).all() for x in arrays):return dict(status='not_evaluable',effect=np.nan,SE=np.nan,df=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 effect=sum(s*x.mean() for s,x in zip(signs,arrays))
 v=np.array([s*s*x.var(ddof=1)/len(x) for s,x in zip(signs,arrays)])
 if v.sum()<=0:return dict(status='zero_variance_no_test',effect=effect,SE=0,df=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 se=np.sqrt(v.sum());df=v.sum()**2/sum(vi**2/(len(x)-1) for vi,x in zip(v,arrays));t=effect/se
 ci=stats.t.ppf(.975,df)*se
 return dict(status='evaluated',effect=effect,SE=se,df=df,P=2*stats.t.sf(abs(t),df),CI_low=effect-ci,CI_high=effect+ci)
def bh(v):
 v=np.asarray(v);ans=np.full(len(v),np.nan);ix=np.flatnonzero(np.isfinite(v));order=ix[np.argsort(v[ix])];m=len(order)
 if m:ans[order]=np.minimum(1,np.minimum.accumulate((v[order]*m/np.arange(1,m+1))[::-1])[::-1])
 return ans
rows=[];combined=[];validation=[]
for st,nn in nci.groupby('stratum'):
 for g in registry:
  gene=g['gene'];x=nn[nn.gene.eq(gene)].expression.to_numpy()
  nb=bone[bone.gene.eq(gene)].expression.to_numpy();nl=lung[lung.gene.eq(gene)].expression.to_numpy()
  o=old[old.gene.eq(gene)]
  arrays=[];weights=[];dfs=[]
  for src,z in o.groupby('source'):
   b=z[z.site.eq('Bone')].expression.to_numpy();l=z[z.site.eq('Lung')].expression.to_numpy()
   for name,aa,ss in [('tumor_only',[x,b],[1,-1]),('normal_min20',[x,b,nl,nb],[1,-1,-1,1])]:
    res=contrast(aa,ss)
    res.update(gene=gene,stratum=st,bone_source=src,contrast=name,NCI_N=len(x),Bone_N=len(b),normal_bone_N=len(nb),normal_lung_N=len(nl),original_effect=g[src+'__mean_lung_minus_bone'],original_P=g[src+'__P'],original_FDR=g[src+'__FDR'])
    res['same_direction']=bool(res['effect']*res['original_effect']>0) if np.isfinite(res['effect']) else False
    rows.append(res)
    if name=='tumor_only':
     ref=stats.ttest_ind(x,b,equal_var=False)
     validation.append(abs(ref.pvalue-res['P']))
   w=len(l)*len(b)/(len(l)+len(b));weights.append(w);arrays.extend([l,b]);dfs.extend([w,-w])
  scale=sum(weights);signs=np.array(dfs)/scale
  for name,aa,ss in [('tumor_only',arrays,signs),('normal_min20',arrays+[nl,nb],list(signs)+[-1,1])]:
   rr=contrast(aa,ss);rr.update(gene=gene,stratum=st,contrast=name,NCI_N=len(x),site_information_weight_NCI=0)
   combined.append(rr)
  # Independent numerical OLS check of site coefficient before and after adding NCI.
  tmp=o.copy();tmp['donor']=tmp.donor.astype(str)
  nc=pd.DataFrame({'source':'NCI','site':'Lung','expression':x})
  coefs=[]
  for dd in [tmp,pd.concat([tmp,nc],ignore_index=True)]:
   X=np.column_stack([np.ones(len(dd)),dd.site.eq('Lung').to_numpy().astype(float),pd.get_dummies(dd.source,drop_first=True,dtype=float).to_numpy()])
   coefs.append(np.linalg.lstsq(X,dd.expression.to_numpy(),rcond=None)[0][1])
  expected=contrast(arrays,signs)['effect']
  assert max(abs(v-expected) for v in coefs)<1e-10
  validation.extend([abs(coefs[0]-coefs[1]),abs(coefs[0]-expected)])
res=pd.DataFrame(rows)
for _,ix in res.groupby(['stratum','bone_source','contrast']).groups.items():
 res.loc[ix,'FDR']=bh(res.loc[ix,'P']);res.loc[ix,'family_valid_tests']=res.loc[ix,'P'].notna().sum()
res['nominal']=res.P.lt(.05);res['FDR_significant']=res.FDR.lt(.05)
res.to_csv(P/'tables/NCI_author_CAF_external_contrasts.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
comb=pd.DataFrame(combined)
for _,ix in comb.groupby(['stratum','contrast']).groups.items():comb.loc[ix,'FDR']=bh(comb.loc[ix,'P'])
comb.to_csv(P/'tables/source_adjusted_site_contrasts.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
summ=[]
for (st,ct),z in res.groupby(['stratum','contrast']):
 pivot=z.pivot(index='gene',columns='bone_source',values='same_direction')
 bothdir=pivot.all(axis=1)
 nom=z.assign(ok=z.same_direction & z.nominal).pivot(index='gene',columns='bone_source',values='ok').all(axis=1)
 fdr=z.assign(ok=z.same_direction & z.FDR_significant).pivot(index='gene',columns='bone_source',values='ok').all(axis=1)
 d=dict(stratum=st,contrast=ct,same_direction_both=bothdir[bothdir].index.tolist(),same_direction_nominal_both=nom[nom].index.tolist(),same_direction_FDR_both=fdr[fdr].index.tolist())
 summ.append(d)
qa=dict(status='PASS',scipy_Welch_and_OLS_checks=len(validation),max_absolute_error=float(max(validation)),normal_bone_donors=int(bone.donor.nunique()),normal_lung_donors=int(lung.donor.nunique()),source_adjustment_NCI_site_weight_zero=True)
(P/'resources/STATISTICAL_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
(P/'resources/RESULT_SUMMARY.json').write_text(json.dumps(summ,indent=2),encoding='utf-8')
print(json.dumps(summ,indent=2));print(json.dumps(qa))
print(res[res.stratum.eq('author_all')][['gene','bone_source','contrast','effect','P','FDR','same_direction']].to_string(index=False))
