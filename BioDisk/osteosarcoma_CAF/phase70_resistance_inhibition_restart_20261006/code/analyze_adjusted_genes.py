from pathlib import Path
import pandas as pd, numpy as np, scipy, scipy.stats as st, statsmodels, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import json, hashlib, datetime
B=Path(__file__).resolve().parents[1]; I=B/'inputs'; T=B/'tables'
proto=json.loads((B/'resources/PHASE70_PROTOCOL_FROZEN.json').read_text(encoding='utf-8'))
for name,sha in proto['input_sha256'].items():assert hashlib.sha256((I/name).read_bytes()).hexdigest()==sha,name
old=pd.read_csv(I/'fixed_program_gene_response_statistics.tsv',sep='\t')
meta=pd.read_csv(I/'patient_metadata.tsv',sep='\t').set_index('gsm',verify_integrity=True)
cand=pd.read_csv(I/'nominal_upregulated_all_sources_with_sensitivity.tsv',sep='\t')
assert len(old)==1303 and not old.duplicated(['cohort','gene']).any()
assert len(cand)==54 and cand.gene.nunique()==39
out=[];qc=[];members=[]; max_p=0.; max_e=0.
for cohort, legacy in old.groupby('cohort',sort=True):
 expr=pd.read_csv(I/f'{cohort}_gene_VST.tsv.gz',sep='\t',index_col=0)
 assert expr.index.is_unique and expr.columns.is_unique
 m=meta.loc[(meta.cohort==cohort)&meta.included_primary.eq(True)].copy()
 assert m.poor.isin([0,1]).all() and set(m.index)<=set(expr.columns)
 cov=['poor','age_years','male']+(['metastases_at_diagnosis'] if cohort=='GSE21257' else [])
 d=m.dropna(subset=cov)
 X=sm.add_constant(d[cov].astype(float),has_constant='add')
 assert np.linalg.matrix_rank(X)==X.shape[1]
 lev=np.diag(X.to_numpy()@np.linalg.pinv(X.to_numpy()))
 qc.append({'cohort':cohort,'n_primary':len(m),'n_adjusted':len(d),'n_poor':int(d.poor.sum()),'n_good':int((d.poor==0).sum()),'covariates':','.join(cov[1:]),'family_size':len(legacy),'model_rank':int(np.linalg.matrix_rank(X)),'condition_number':float(np.linalg.cond(X)),'max_leverage':float(lev.max()),'excluded_for_missing_covariates':','.join(m.index.difference(d.index)),'chemotherapy_counts':json.dumps(m.chemotherapy.fillna('not recorded').value_counts().to_dict())})
 for gsm in m.index:members.append({'cohort':cohort,'gsm':gsm,'poor':m.loc[gsm,'poor'],'included_adjusted':gsm in d.index})
 for row in legacy.itertuples():
  assert row.gene in expr.index
  a=expr.loc[row.gene,m.index[m.poor==1]].to_numpy(float);b=expr.loc[row.gene,m.index[m.poor==0]].to_numpy(float)
  assert np.isfinite(a).all() and np.isfinite(b).all()
  assert len(a)==row.n_poor and len(b)==row.n_good
  p=st.ttest_ind(a,b,equal_var=False).pvalue;effect=a.mean()-b.mean()
  max_p=max(max_p,abs(p-row.p));max_e=max(max_e,abs(effect-row.effect))
  y=expr.loc[row.gene,d.index].astype(float)
  assert np.isfinite(y).all() and y.var()>0
  fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['poor']
  out.append({'cohort':cohort,'gene':row.gene,'model':'adjusted_OLS_HC3_t','units':'submitted_gene_VST','n_poor':int(d.poor.sum()),'n_good':int((d.poor==0).sum()),'effect':fit.params['poor'],'ci_low':ci.iloc[0],'ci_high':ci.iloc[1],'p':fit.pvalues['poor'],'standard_error_HC3':fit.bse['poor'],'df':fit.df_resid,'covariates':','.join(cov[1:]),'family':f'{cohort}:all_measured_fixed_program_genes:adjusted_HC3'})
assert max_p<1e-8 and max_e<1e-8,(max_p,max_e)
adj=pd.DataFrame(out)
assert len(adj)==1303 and np.isfinite(adj[['effect','p','ci_low','ci_high']]).all().all()
for cohort,idx in adj.groupby('cohort').groups.items():adj.loc[idx,'fdr_gene_family']=multipletests(adj.loc[idx,'p'],method='fdr_bh')[1]
adj['nominal']=adj.p<.05;adj['FDR_significant']=adj.fdr_gene_family<.05
adj.to_csv(T/'all_fixed_family_adjusted_gene_response.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(T/'clinical_model_QC.tsv',sep='\t',index=False)
pd.DataFrame(members).to_csv(T/'clinical_model_membership.tsv',sep='\t',index=False)
paired=old.merge(adj,on=['cohort','gene'],how='outer',suffixes=('_unadjusted','_adjusted'),validate='one_to_one',indicator=True)
assert paired['_merge'].eq('both').all();paired=paired.drop(columns='_merge')
paired['same_direction']=np.sign(paired.effect_unadjusted)==np.sign(paired.effect_adjusted)
paired['nominal_both_models']=paired.p_unadjusted.lt(.05)&paired.p_adjusted.lt(.05)
paired['positive_nominal_both_models']=paired.nominal_both_models&paired.effect_unadjusted.gt(0)&paired.effect_adjusted.gt(0)
paired.to_csv(T/'all_fixed_family_unadjusted_adjusted_comparison.tsv',sep='\t',index=False)
grid=pd.MultiIndex.from_product([sorted(cand.gene.unique()),sorted(old.cohort.unique())],names=['gene','cohort']).to_frame(index=False)
gene=grid.merge(paired,on=['gene','cohort'],how='left',validate='one_to_one')
gene['clinical_assayed']=gene.p_adjusted.notna()
gene.to_csv(T/'CAF_39_genes_clinical_evidence.tsv',sep='\t',index=False)
record=cand.merge(gene,on='gene',how='left',validate='many_to_many')
assert len(record)==108
record.to_csv(T/'CAF_54_records_clinical_evidence.tsv',sep='\t',index=False)
cross=[]
for g,v in gene.groupby('gene'):
 cross.append({'gene':g,'cohorts_assayed':int(v.clinical_assayed.sum()),'both_cohorts_positive_adjusted_nominal':bool(len(v)==2 and v.p_adjusted.lt(.05).all() and v.effect_adjusted.gt(0).all()),'both_cohorts_any_same_direction_adjusted_nominal':bool(v.p_adjusted.lt(.05).all() and ((v.effect_adjusted>0).all() or (v.effect_adjusted<0).all())),'cohorts_positive_both_models':','.join(v.loc[v.positive_nominal_both_models.eq(True),'cohort'])})
cross=pd.DataFrame(cross);cross.to_csv(T/'CAF_cross_cohort_summary.tsv',sep='\t',index=False)
result={'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'completed','protocol_sha256':hashlib.sha256((B/'resources/PHASE70_PROTOCOL_FROZEN.json').read_bytes()).hexdigest(),'N_patients_primary':len(members),'N_patients_adjusted':sum(q['n_adjusted'] for q in qc),'gene_cohort_tests':len(adj),'candidate_genes':39,'candidate_records':54,'record_cohort_rows':len(record),'reproduction_max_abs_p_error':max_p,'reproduction_max_abs_effect_error':max_e,'all_gene_adjusted_nominal_by_cohort':adj[adj.nominal].groupby('cohort').size().to_dict(),'all_gene_adjusted_FDR_significant':int(adj.FDR_significant.sum()),'candidate_missing':gene.loc[~gene.clinical_assayed,['cohort','gene']].to_dict('records'),'candidate_adjusted_nominal':gene.loc[gene.p_adjusted.lt(.05),['cohort','gene','effect_adjusted','ci_low_adjusted','ci_high_adjusted','p_adjusted','fdr_gene_family_adjusted','p_unadjusted','effect_unadjusted']].to_dict('records'),'candidate_both_cohort_nominal':cross.loc[cross.both_cohorts_any_same_direction_adjusted_nominal].to_dict('records'),'versions':{'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__}}
(B/'resources/ANALYSIS_COMPLETE.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
