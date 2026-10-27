"""Rebuild source-adjusted informative-source-df results from preserved trial."""
import pathlib,json
import pandas as pd,numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
P=pathlib.Path(__file__).resolve().parents[1];R=P.parent
x=pd.read_csv(P/'tables/all_14gene_results_global_df_trial.tsv',sep='\t')
extra=x[x.model.eq('source_adjusted_OLS_HC3')].copy();extra['model']='source_adjusted_OLS_HC3_global_df';extra['family']=extra.family.str.replace('source_adjusted_OLS_HC3','source_adjusted_OLS_HC3_global_df',regex=False)
refs=[]
for name in ['bone_reference_with_coverage_sensitivity.tsv','lung_reference_pseudobulk.tsv']:
 z=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables'/name,sep='\t');refs.append(z[z.specification.eq('min20')])
for i,r in x[x.model.eq('source_adjusted_OLS_HC3') & x.status.eq('evaluated')].iterrows():
 base=x[x.cellset.eq(r.cellset)&x.cohort.eq('without_NCI')&x.model.eq(r.model)&x.contrast.eq('tumor_only')&x.gene.eq(r.gene)].iloc[0]
 df=base.df
 if r.contrast=='normal_min20':
  b,l=[z[z.gene.eq(r.gene)].expression.to_numpy() for z in refs]
  v=[base.SE**2,l.var(ddof=1)/len(l),b.var(ddof=1)/len(b)]
  df=sum(v)**2/(v[0]**2/base.df+v[1]**2/(len(l)-1)+v[2]**2/(len(b)-1))
 ci=stats.t.ppf(.975,df)*r.SE
 x.loc[i,['df','P','CI_low','CI_high']]=[df,2*stats.t.sf(abs(r.effect/r.SE),df),r.effect-ci,r.effect+ci]
x=pd.concat([x,extra],ignore_index=True)
for _,idx in x.groupby('family').groups.items():
 ix=x.loc[idx].index[x.loc[idx,'P'].notna()];x.loc[ix,'FDR']=multipletests(x.loc[ix,'P'],method='fdr_bh')[1]
x['nominal']=x.P.lt(.05);x['FDR_significant']=x.FDR.lt(.05)
x.to_csv(P/'tables/all_14gene_results.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')

rows=[]
for _,z in x.groupby('family',sort=False):
 r=z.iloc[0];d={k:r[k] for k in ['cellset','cohort','model','contrast','Bone_N','Lung_N','Bone_cells','Lung_cells']}
 for key,mask in [('same_direction',z.same_direction),('same_direction_nominal',z.same_direction&z.nominal),('same_direction_FDR',z.same_direction&z.FDR_significant),('opposite_nominal',~z.same_direction&z.nominal)]:
  d[key+'_n']=int(mask.sum());d[key+'_genes']=';'.join(z.loc[mask,'gene'])
 rows.append(d)
pd.DataFrame(rows).to_csv(P/'tables/scenario_summary.tsv',sep='\t',encoding='utf-8',index=False)
