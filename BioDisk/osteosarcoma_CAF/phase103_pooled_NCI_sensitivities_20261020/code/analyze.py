import pathlib,json
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase103_pooled_NCI_sensitivities_20261020'
dat=pd.read_csv(P/'tables/all_patient_expression.tsv',sep='\t')
registry=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text())
nb=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/bone_reference_with_coverage_sensitivity.tsv',sep='\t');nb=nb[nb.specification.eq('min20')]
nl=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/lung_reference_pseudobulk.tsv',sep='\t');nl=nl[nl.specification.eq('min20')]
assert nb.donor.nunique()==5 and nl.donor.nunique()==3 and nb.cells.min()>=20 and nl.cells.min()>=20
protocol=dict(primary='expanded + NCI, pooled unpaired patient-level Welch (16 Bone vs 23 Lung), exploratory',
cellsets=dict(expanded='Original two main cohorts expanded candidates + NCI author CAF',
singlets='Original doublet_positive_any excluded + NCI author singlets',
conservative='Original frozen conservative_fibro and no doublet + NCI author singlets',
expanded_matched_strict_donors='Expanded expression restricted to patients retained in conservative set, NCI author all; controls patient-roster shift'),
models=['pooled_Welch','source_adjusted_OLS_HC3'],
cohorts=['with_NCI','without_NCI'],
normal='Separate raw site contrast and subtract normal Lung - normal Bone; retain >=20 cells per normal donor',
source_adjusted='OLS intercept + Lung indicator + source dummy variables; target Lung coefficient, HC3 covariance, residual-df t inference',
normal_variance='Independent normal donor variance/n added to tumor contrast variance; Satterthwaite df combining tumor variance with its df and each reference variance component',
FDR='BH separately within each cellset x cohort x model x contrast across valid prespecified14 genes (12 for normal)',
sample_scope='Two original main sources SCPCP000017 and SCPCP000023 only; original paired-donor exclusions preserved; source-aware donor ID matching; no shared public IDs, no genotype cross-study identity proof',
no_new_patient_cell_threshold=True,cell_count_is_not_biological_N=True,all_zero_or_no_variance='effect only, P=NA',
source_adjustment_limit='NCI source has only Lung; source-adjusted site coefficient unchanged by adding NCI, though HC3 df and variance may change.',
old_genomewide_FDR_preserved=True)
(P/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
def stat(effect,vs,dfs):
 vs=np.asarray(vs,dtype=float);dfs=np.asarray(dfs,dtype=float)
 if not np.isfinite(effect) or not np.isfinite(vs).all() or not np.isfinite(dfs).all() or (dfs<=0).any():
  return dict(status='not_evaluable',effect=effect,SE=np.nan,df=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 v=vs.sum()
 if v<1e-24:return dict(status='zero_variance_effect_only',effect=effect,SE=0,df=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 df=v*v/np.sum(vs*vs/dfs);se=np.sqrt(v);ci=stats.t.ppf(.975,df)*se
 return dict(status='evaluated',effect=float(effect),SE=float(se),df=float(df),P=float(2*stats.t.sf(abs(effect/se),df)),CI_low=float(effect-ci),CI_high=float(effect+ci))
rows=[];qa=[];coeff=[]
for cellset,dd in dat.groupby('cellset',sort=False):
 for cohort in ['with_NCI','without_NCI']:
  ds=dd if cohort=='with_NCI' else dd[dd.source.ne('NCI')]
  for gene in registry:
   g=gene['gene'];d=ds[ds.gene.eq(g)].copy();d=d[np.isfinite(d.expression)]
   y=d.expression.to_numpy();b=d[d.site.eq('Bone')].expression.to_numpy();l=d[d.site.eq('Lung')].expression.to_numpy()
   assert len(b)>=2 and len(l)>=2
   X=np.column_stack([np.ones(len(d)),d.site.eq('Lung').to_numpy().astype(float),pd.get_dummies(d.source,drop_first=True,dtype=float).to_numpy()])
   fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True)
   # Independent sandwich check.
   inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;h=np.einsum('ij,jk,ik->i',X,inv,X)
   assert np.max(h)<1-1e-10
   cov=inv@(X.T@(((y-X@beta)**2/(1-h)**2)[:,None]*X))@inv
   qa.extend([abs(beta[1]-fit.params[1]),abs(cov[1,1]-fit.cov_params()[1,1])])
   coeff.append(dict(cellset=cellset,cohort=cohort,gene=g,effect=fit.params[1]))
   v0=[l.var(ddof=1)/len(l),b.var(ddof=1)/len(b)];df0=[len(l)-1,len(b)-1]
   models=[('pooled_Welch',l.mean()-b.mean(),v0,df0),('source_adjusted_OLS_HC3',fit.params[1],[fit.cov_params()[1,1]],[fit.df_resid])]
   base=stat(l.mean()-b.mean(),v0,df0);ref=stats.ttest_ind(l,b,equal_var=False)
   if np.isfinite(ref.pvalue):qa.append(abs(base['P']-ref.pvalue))
   bn=nb[nb.gene.eq(g)].expression.to_numpy();ln=nl[nl.gene.eq(g)].expression.to_numpy()
   for model,effect,vs,dfs in models:
    for contrast in ['tumor_only','normal_min20']:
     if contrast=='tumor_only':out=stat(effect,vs,dfs)
     elif np.isfinite(bn).all() and np.isfinite(ln).all():
      out=stat(effect-(ln.mean()-bn.mean()),list(vs)+[ln.var(ddof=1)/len(ln),bn.var(ddof=1)/len(bn)],list(dfs)+[len(ln)-1,len(bn)-1])
     else:out=dict(status='reference_feature_missing',effect=np.nan,SE=np.nan,df=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
     out.update(cellset=cellset,cohort=cohort,model=model,contrast=contrast,gene=g,Bone_N=len(b),Lung_N=len(l),Bone_cells=int(d[d.site.eq('Bone')].cells.sum()),Lung_cells=int(d[d.site.eq('Lung')].cells.sum()),original_direction='Bone_higher' if g=='EBF1' else 'Lung_higher')
     out['same_direction']=bool(out['effect']<0 if g=='EBF1' else out['effect']>0) if np.isfinite(out['effect']) else False
     out['family']='|'.join([cellset,cohort,model,contrast]);rows.append(out)
out=pd.DataFrame(rows)
for fam,ix in out.groupby('family').groups.items():
 valid=out.loc[ix,'P'].notna();sel=out.loc[ix].index[valid]
 out.loc[ix,'family_valid_tests']=len(sel);out.loc[sel,'FDR']=multipletests(out.loc[sel,'P'],method='fdr_bh')[1]
out['nominal']=out.P.lt(.05);out['FDR_significant']=out.FDR.lt(.05)
out.to_csv(P/'tables/all_14gene_results.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
summ=[]
for fam,z in out.groupby('family',sort=False):
 r=z.iloc[0]
 d={k:r[k] for k in ['cellset','cohort','model','contrast','Bone_N','Lung_N','Bone_cells','Lung_cells']}
 for key,mask in [('same_direction',z.same_direction),('same_direction_nominal',z.same_direction&z.nominal),('same_direction_FDR',z.same_direction&z.FDR_significant),('opposite_nominal',~z.same_direction&z.nominal)]:
  d[key+'_n']=int(mask.sum());d[key+'_genes']=';'.join(z.loc[mask,'gene'])
 summ.append(d)
summary=pd.DataFrame(summ);summary.to_csv(P/'tables/scenario_summary.tsv',sep='\t',encoding='utf-8',index=False)
co=pd.DataFrame(coeff).pivot(index=['cellset','gene'],columns='cohort',values='effect');diff=(co.with_NCI-co.without_NCI).abs().max()
assert diff<1e-10 and max(qa)<1e-10
qaout=dict(status='PASS',results=len(out),families=out.family.nunique(),independent_checks=len(qa),max_error=float(max(qa)),source_adjusted_coefficient_change_max=float(diff),statsmodels_version=__import__('statsmodels').__version__)
(P/'resources/STATISTICAL_QA.json').write_text(json.dumps(qaout,indent=2),encoding='utf-8')
print(summary[['cellset','cohort','model','contrast','Bone_N','Lung_N','same_direction_nominal_genes','same_direction_FDR_genes','opposite_nominal_genes']].to_string(index=False));print(json.dumps(qaout))
