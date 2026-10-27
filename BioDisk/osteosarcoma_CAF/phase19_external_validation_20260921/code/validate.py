from pathlib import Path
import json,gzip,csv,io,hashlib,datetime,os
import pandas as pd,numpy as np
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(os.environ.get('CAF_VALIDATION_ROOT','C:/codex_test/CAF_external_delivery_20260921'));I=O/'inputs';T=O/'tables'
protocol=json.loads((O/'FROZEN_PROTOCOL.json').read_text());sigs=protocol['signatures'];rng=np.random.default_rng(20260921)
def save(d,n):d.to_csv(T/n,sep='\t',index=False)
def zscore(v):return (v-v.mean())/v.std(ddof=1)
X1=pd.read_csv(I/'GSE21257_gene_median_expression.tsv.gz',sep='\t',index_col=0)
C1=pd.read_csv(I/'GSE21257_clinical_mapping.tsv',sep='\t').set_index('sample');assert X1.columns.is_unique and C1.index.is_unique
assert set(X1.columns)==set(C1.index) and len(C1)==53
X1=X1.loc[:,C1.index];C1['event']=C1.metastasis.astype(int)
raw=gzip.open(I/'GSE39055_series_matrix.local.txt.gz','rt').read();ls=raw.splitlines();a=ls.index('!series_matrix_table_begin');b=ls.index('!series_matrix_table_end')
PX=pd.read_csv(io.StringIO('\n'.join(ls[a+1:b])),sep='\t',index_col=0)
samples=next(csv.reader([next(x for x in ls if x.startswith('!Sample_geo_accession'))],delimiter='\t'))[1:]
C2=pd.DataFrame(index=samples)
for line in ls:
 if line.startswith('!Sample_characteristics_ch1'):
  vals=next(csv.reader([line],delimiter='\t'))[1:];keys=[x.split(': ',1)[0] for x in vals];assert len(set(keys))==1
  C2[keys[0]]=[x.split(': ',1)[1] for x in vals]
C2['title']=next(csv.reader([next(x for x in ls if x.startswith('!Sample_title'))],delimiter='\t'))[1:]
C2['time_months']=pd.to_numeric(C2['time until first recurrence or latest follow-up (months)']);assert set(C2.recurrence)=={'Y','N'}
C2['event']=C2.recurrence.map({'Y':1,'N':0});C2['age_years']=C2.age.str.extract(r'(\d+)')[0].astype(float);C2['male']=C2.gender.map({'male':1,'female':0});C2['eligible_RFS']=C2.time_months.gt(0)&C2.event.notna()
ann=pd.read_csv(I/'GPL14951_platform_table.tsv',sep='\t',dtype=str).set_index('ID');assert len(ann)==29377 and ann.index.is_unique
PX.index=PX.index.astype(str);symbols=ann.Symbol.reindex(PX.index);assert symbols.notna().mean()>.8
valid=symbols.notna()&~symbols.str.contains('///|;',na=False)&symbols.ne('')
X2=PX.loc[valid].groupby(symbols[valid]).median();assert set(X2.columns)==set(C2.index) and len(C2)==37
X2=X2.loc[:,C2.index];X2.to_csv(I/'GSE39055_gene_median_expression.tsv.gz',sep='\t')
assert np.isfinite(X1.to_numpy()).all() and np.isfinite(X2.to_numpy()).all()
assert not set(C1.index)&set(C2.index)
C1.to_csv(T/'GSE21257_clinical_used.tsv',sep='\t',index_label='sample');C2.to_csv(T/'GSE39055_clinical_used.tsv',sep='\t',index_label='sample')
coverage=[];scoretables=[];coherence=[];fits=[];sens=[];cohorts=[('GSE21257',X1,C1),('GSE39055',X2,C2)]
def scores(X,cohort):
 out=pd.DataFrame(index=X.columns)
 Z=X.sub(X.mean(axis=1),axis=0).div(X.std(axis=1,ddof=1).replace(0,np.nan),axis=0)
 ranks=X.rank(axis=0,pct=True)
 for name,genes in sigs.items():
  measured=[g for g in genes if g in X.index and X.loc[g].std(ddof=1)>0];cov=len(measured)/len(genes);assert cov>=.8
  coverage.append({'cohort':cohort,'signature':name,'requested_genes':len(genes),'measured_variable_genes':len(measured),'coverage':cov,'missing_genes':';'.join(set(genes)-set(measured))})
  out[name]=zscore(Z.loc[measured].mean());out[name+'_rank']=zscore(ranks.loc[measured].mean())
  for gene in measured:
   rest=[x for x in measured if x!=gene];corr=np.corrcoef(Z.loc[gene],Z.loc[rest].mean())[0,1]
   coherence.append({'cohort':cohort,'signature':name,'gene':gene,'gene_vs_rest_signature_Pearson':corr,'evidence':'descriptive tissue coexpression; not cell identity'})
 return out
def logistic(score,covars,y,cohort,name,model):
 xx=pd.DataFrame({'score':score});
 if covars is not None:xx=xx.join(covars)
 xx=sm.add_constant(xx).astype(float);yy=y.reindex(xx.index).astype(float);ok=xx.notna().all(axis=1)&yy.notna();xx=xx[ok];yy=yy[ok]
 res=sm.Logit(yy,xx).fit(disp=False,maxiter=200);assert res.mle_retvals['converged'] and np.isfinite(res.params).all()
 ci=res.conf_int().loc['score'];beta=float(res.params['score']);se=float(res.bse['score'])
 return {'cohort':cohort,'signature':name,'endpoint':'any_metastasis_original_group','model':model,'n':len(yy),'events':int(yy.sum()),'non_events':int((1-yy).sum()),'effect_type':'OR_per_score_SD','effect':np.exp(beta),'CI_low':np.exp(ci.iloc[0]),'CI_high':np.exp(ci.iloc[1]),'beta':beta,'SE':se,'P':float(res.pvalues['score']),'converged':True}
for cohort,X,C in cohorts:
 sc=scores(X,cohort);sc['cohort']=cohort;sc.index.name='sample';scoretables.append(sc.reset_index());sc=sc.drop(columns='cohort')
 if cohort=='GSE21257':
  for name in sigs:
   row=logistic(sc[name],None,C.event,cohort,name,'primary_unadjusted');yy=C.event;row['AUC_high_score_for_event']=roc_auc_score(yy,sc[name]);row['mean_difference_SD']=sc.loc[yy==1,name].mean()-sc.loc[yy==0,name].mean()
   event=sc.loc[yy==1,name].to_numpy();non=sc.loc[yy==0,name].to_numpy();boot=rng.choice(event,(2000,len(event)),replace=True).mean(1)-rng.choice(non,(2000,len(non)),replace=True).mean(1);row['mean_diff_CI_low'],row['mean_diff_CI_high']=np.quantile(boot,[.025,.975]);fits.append(row)
   subset=~C.baseline_metastasis
   sens.append(logistic(sc.loc[subset,name],None,C.loc[subset,'event'],cohort,name,'exclude_baseline_metastasis'))
   sens.append(logistic(sc[name],C[['age_years','male']],C.event,cohort,name,'adjust_age_sex'))
   sens.append(logistic(sc[name+'_rank'],None,C.event,cohort,name,'rank_score_sensitivity'))
 else:
  cc=C[C.eligible_RFS]
  for name in sigs:
   ss=sc.loc[cc.index,name];res=PHReg(cc.time_months.to_numpy(),ss.to_numpy()[:,None],status=cc.event.to_numpy(),ties='efron').fit()
   assert np.isfinite(res.params).all() and np.isfinite(res.bse).all()
   beta=float(res.params[0]);ci=res.conf_int()[0]
   fits.append({'cohort':cohort,'signature':name,'endpoint':'recurrence_free_followup','model':'primary_Cox_Efron','n':len(cc),'events':int(cc.event.sum()),'non_events':int((1-cc.event).sum()),'effect_type':'HR_per_score_SD','effect':np.exp(beta),'CI_low':np.exp(ci[0]),'CI_high':np.exp(ci[1]),'beta':beta,'SE':float(res.bse[0]),'P':float(res.pvalues[0]),'converged':True})
   # No OS model: the available time field is RFS, not death/censoring time.
F=pd.DataFrame(fits);F['FDR']=np.nan
for _,ix in F.groupby(['cohort','endpoint']).groups.items():F.loc[ix,'FDR']=multipletests(F.loc[ix,'P'],method='fdr_bh')[1]
S=pd.DataFrame(sens);S['FDR_sensitivity_family']=np.nan
for _,ix in S.groupby('model').groups.items():S.loc[ix,'FDR_sensitivity_family']=multipletests(S.loc[ix,'P'],method='fdr_bh')[1]
save(F,'external_primary_associations.tsv');save(S,'external_sensitivity_associations.tsv');save(pd.DataFrame(coverage),'external_signature_coverage.tsv');save(pd.DataFrame(coherence),'external_gene_coherence.tsv');save(pd.concat(scoretables),'external_patient_scores.tsv')
coh=pd.DataFrame(coherence).groupby(['cohort','signature']).agg(genes=('gene','size'),median_gene_rest_correlation=('gene_vs_rest_signature_Pearson','median'),positive_gene_correlations=('gene_vs_rest_signature_Pearson',lambda v:int((v>0).sum()))).reset_index();save(coh,'external_coherence_summary.tsv')
fig,axs=plt.subplots(1,2,figsize=(10,4.2))
for ax,cohort in zip(axs,['GSE21257','GSE39055']):
 rr=F[F.cohort==cohort].reset_index(drop=True)
 for i,r in rr.iterrows():ax.errorbar(r.effect,i,xerr=[[r.effect-r.CI_low],[r.CI_high-r.effect]],fmt='o',color='#175d82',capsize=4);ax.text(.02,.08+.09*i,f'{r.signature}: P={r.P:.3g}, FDR={r.FDR:.3g}',transform=ax.transAxes,fontsize=9)
 ax.axvline(1,color='gray',linestyle='--');ax.set_xscale('log');ax.set_yticks(range(len(rr)),rr.signature);ax.set_ylim(-.6,1.7);ax.set_title(cohort+'  n='+str(int(rr.n.iloc[0])));ax.set_xlabel('Odds ratio / 1 SD score' if cohort=='GSE21257' else 'Hazard ratio / 1 SD score')
fig.tight_layout();fig.savefig(O/'figures/external_validation_forest.png',dpi=180);plt.close(fig)
accepted=[I/'GSE21257_gene_median_expression.tsv.gz',I/'GSE21257_clinical_mapping.tsv',I/'GSE39055_series_matrix.local.txt.gz',I/'GPL14951_platform_table.tsv']
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in accepted};(I/'accepted_input_sha256.json').write_text(json.dumps(manifest,indent=2))
qa={'GSE21257_patients':len(C1),'GSE21257_metastasis':int(C1.event.sum()),'GSE21257_baseline_metastasis':int(C1.baseline_metastasis.sum()),'GSE39055_biopsies':len(C2),'GSE39055_RFS_evaluable':int(C2.eligible_RFS.sum()),'GSE39055_RFS_events':int(C2.loc[C2.eligible_RFS,'event'].sum()),'GSE39055_RFS_zero_time_excluded':C2.index[~C2.eligible_RFS].tolist(),'OS_model_not_run':'death status exists but a valid OS time field was not provided','signatures_frozen_before_testing':True,'outcome_driven_gene_selection':False,'both_signature_coverage_at_least_80_percent':True,'FDR_families':'2 signatures separately per cohort-specific endpoint; sensitivities separately by model','validation_scope':'external tissue outcome associations, not single-cell identity or protein sorting','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'delivery/VALIDATION_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2));(O/'VALIDATION_COMPLETE').write_text('EXTERNAL_TISSUE_VALIDATION_COMPLETE\n')
print(F.to_string(index=False),flush=True);print(S.to_string(index=False),flush=True);print(json.dumps(qa),flush=True)
