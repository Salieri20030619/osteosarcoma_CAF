from pathlib import Path
import json,csv,gzip,io,re,hashlib,datetime
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005/interim_GSE21257'
P=json.loads((B/'phase63_matrix_resistance_20261005/resources/PROGRAMS.json').read_text())['sets']
P.update(json.loads((B/'phase63_matrix_resistance_20261005/resources/POSTHOC_COLLAGEN_STRUCTURE_PROTOCOL.json').read_text())['panels'])
assert len(P)==14
(O/'resources/FIXED_PROGRAMS.json').write_text(json.dumps(P,indent=2))
h=pd.read_csv(B/'phase41_full_transcriptome_setup_20261001/resources/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
approved=set(h.symbol);aliases={}
for _,r in h.iterrows():
 for col in ['prev_symbol','alias_symbol']:
  for a in r[col].split('|'):
   if a:aliases.setdefault(a,set()).add(r.symbol)
def canon(s):
 if s in approved:return s
 a=aliases.get(s,set())
 return next(iter(a)) if len(a)==1 else None
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def welch(y,label):
 a=y[label==1];b=y[label==0];na=len(a);nb=len(b)
 if min(na,nb)<2 or np.var(y)==0:return None
 effect=a.mean()-b.mean();se=np.sqrt(a.var(ddof=1)/na+b.var(ddof=1)/nb)
 df=(a.var(ddof=1)/na+b.var(ddof=1)/nb)**2/((a.var(ddof=1)/na)**2/(na-1)+(b.var(ddof=1)/nb)**2/(nb-1))
 pool=np.sqrt(((na-1)*a.var(ddof=1)+(nb-1)*b.var(ddof=1))/(na+nb-2))
 return dict(n_poor=na,n_good=nb,effect=effect,ci_low=effect-stats.t.ppf(.975,df)*se,ci_high=effect+stats.t.ppf(.975,df)*se,p=stats.ttest_ind(a,b,equal_var=False).pvalue,hedges_g=effect/pool*(1-3/(4*(na+nb)-9)),df=df)
allstats=[];allgenes=[];coverage=[];qcs=[];allmeta=[];allscores=[]
for cohort,gpl in [('GSE21257','GPL10295')]:
 print('START',cohort,flush=True)
 raw=pd.DataFrame(json.loads((O/'resources'/f'{cohort}_clinical_raw.json').read_text())).set_index('gsm')
 raw['age_years']=raw.age.str.extract(r'(\d+(?:\.\d+)?)').astype(float)
 if cohort=='GSE21257':raw['age_years']/=12
 raw['male']=raw.gender.str.lower().isin(['m','male']).astype(int)
 if cohort=='GSE39055':
  def nec(s):
   nums=list(map(float,re.findall(r'\d+(?:\.\d+)?',s)))
   lo=nums[0];hi=100 if '>' in s else nums[-1]
   return 0. if lo>=90 else (1. if hi<90 else np.nan)
  raw['poor']=raw['percent necrosis'].map(nec);raw['response_raw']=raw['percent necrosis']
 else:
  raw['poor']=raw['huvos grade'].map({'1':1.,'2':1.,'3':0.,'4':0.});raw['response_raw']=raw['huvos grade']
  raw['metastases_at_diagnosis']=(raw['group']=='Metastases present at diagnosis').astype(int)
 raw['cohort']=cohort;raw['included_primary']=raw.poor.notna()
 assert raw.index.is_unique and raw.title.is_unique
 allmeta.append(raw.reset_index())
 text=gzip.open(O/'sources'/f'{cohort}_series_matrix.txt.gz','rt').read()
 table=text.split('!series_matrix_table_begin')[1].split('!series_matrix_table_end')[0].strip()
 x=pd.read_csv(io.StringIO(table),sep='\t',index_col=0);x.index=x.index.astype(str)
 assert list(x.columns)==list(raw.index) and x.index.is_unique and np.isfinite(x.to_numpy()).all()
 atext=(O/'sources'/f'{gpl}_platform_only.soft').read_text()
 atable=atext.split('!platform_table_begin')[1].split('!platform_table_end')[0].strip()
 ann=pd.read_csv(io.StringIO(atable),sep='\t',dtype=str,keep_default_na=False)
 assert ann.ID.is_unique
 ann['canonical_symbol']=ann.Symbol.map(canon)
 mapping=ann[['ID','Symbol','canonical_symbol']].copy()
 mapping['in_matrix']=mapping.ID.isin(x.index)
 save(mapping,f'{cohort}_probe_mapping.tsv')
 mapped=x.join(ann.set_index('ID').canonical_symbol).dropna(subset=['canonical_symbol'])
 genes=mapped.groupby('canonical_symbol',sort=True)[x.columns].median()
 genes.to_csv(O/'objects'/f'{cohort}_gene_VST.tsv.gz',sep='\t',compression='gzip')
 sd=genes.std(axis=1,ddof=1)
 z=genes.sub(genes.mean(axis=1),axis=0).div(sd,axis=0)
 assert not z.loc[sd>0].isna().any().any()
 scores=pd.DataFrame(index=x.columns)
 for name,gs in P.items():
  use=sorted(set(gs)&set(genes.index[sd>0]));frac=len(use)/len(gs);ok=frac>=.8 and len(use)>=3
  coverage.append(dict(cohort=cohort,program=name,n_expected=len(gs),n_measured=len(use),coverage=frac,evaluable=ok,measured=';'.join(use),missing=';'.join(sorted(set(gs)-set(use)))))
  if ok:scores[name]=z.loc[use].mean(axis=0)
 scores.index.name='gsm'
 alls=scores.reset_index();alls['cohort']=cohort;allscores.append(alls)
 keep=raw.poor.notna();labels=raw.loc[keep,'poor'].to_numpy()
 covars=['poor','age_years','male']+(['metastases_at_diagnosis'] if cohort=='GSE21257' else [])
 for name in scores:
  y=scores.loc[keep,name].to_numpy();w=welch(y,labels)
  if w:allstats.append(dict(cohort=cohort,program=name,model='unadjusted_Welch',family='primary' if name=='fibrillar_I_III_V' else 'secondary13',**w))
  design=sm.add_constant(raw.loc[keep,covars].astype(float))
  assert np.linalg.matrix_rank(design)==design.shape[1]
  fit=sm.OLS(y,design).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['poor']
  allstats.append(dict(cohort=cohort,program=name,model='adjusted_HC3',family='primary' if name=='fibrillar_I_III_V' else 'secondary13',n_poor=int(sum(labels==1)),n_good=int(sum(labels==0)),effect=fit.params['poor'],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues['poor'],df=fit.df_resid))
  if name=='fibrillar_I_III_V':
   a=y[labels==1];b=y[labels==0]
   perm=stats.permutation_test((a,b),lambda a,b,axis=0:np.mean(a,axis=axis)-np.mean(b,axis=axis),vectorized=True,n_resamples=99999,batch=1000,random_state=20261005,alternative='two-sided')
   allstats.append(dict(cohort=cohort,program=name,model='permutation_99999',family='primary_sensitivity',n_poor=len(a),n_good=len(b),effect=a.mean()-b.mean(),p=perm.pvalue))
   if cohort=='GSE39055':
    for lab in [0.,1.]:
     ss=raw.poor.fillna(lab).to_numpy();ww=welch(scores[name].to_numpy(),ss)
     allstats.append(dict(cohort=cohort,program=name,model=f'ambiguous_as_{"poor" if lab else "good"}',family='label_sensitivity',**ww))
 union=sorted(set().union(*map(set,P.values()))&set(genes.index[sd>0]))
 for gene in union:
  w=welch(genes.loc[gene,keep.to_numpy()].to_numpy(),labels)
  if w:allgenes.append(dict(cohort=cohort,gene=gene,model='unadjusted_Welch',units='submitted_VST',**w))
 qcs.append(dict(cohort=cohort,n_samples=len(raw),poor=int(sum(labels==1)),good=int(sum(labels==0)),unclassified=int(sum(~keep)),probes=len(x),matched_probes=len(mapped),mapped_genes=len(genes),programs_evaluable=len(scores),raw_min=float(x.min().min()),raw_max=float(x.max().max()),age_range=[float(raw.age_years.min()),float(raw.age_years.max())]))
 print('READY',qcs[-1],flush=True)
# Correct within complete predefined families, never mix genes and pathways.
st=pd.DataFrame(allstats);st['fdr_secondary13']=np.nan;st['fdr_all14']=np.nan
for (c,m),ix in st.groupby(['cohort','model']).groups.items():
 if m not in ['unadjusted_Welch','adjusted_HC3']:continue
 st.loc[ix,'fdr_all14']=multipletests(st.loc[ix,'p'],method='fdr_bh')[1]
 sec=st.index.isin(ix)&st.family.eq('secondary13')
 st.loc[sec,'fdr_secondary13']=multipletests(st.loc[sec,'p'],method='fdr_bh')[1]
ge=pd.DataFrame(allgenes);ge['fdr_gene_family']=np.nan
for c,ix in ge.groupby('cohort').groups.items():ge.loc[ix,'fdr_gene_family']=multipletests(ge.loc[ix,'p'],method='fdr_bh')[1]
save(st,'program_response_statistics.tsv');save(ge,'fixed_program_gene_response_statistics.tsv');save(pd.DataFrame(coverage),'program_coverage.tsv');save(pd.concat(allmeta,ignore_index=True),'patient_metadata.tsv');scoresall=pd.concat(allscores,ignore_index=True);save(scoresall,'patient_program_scores.tsv')
(O/'QA/COMPUTATION_QC.json').write_text(json.dumps(qcs,indent=2))
meta=pd.concat(allmeta).set_index(['cohort','gsm']);fig,axs=plt.subplots(1,2,figsize=(9,4.2),layout='constrained');rng=np.random.default_rng(20261005)
for ax,c in zip(axs,['GSE21257']):
 ss=scoresall[scoresall.cohort==c].set_index('gsm');mm=meta.loc[c].loc[ss.index];keep=mm.poor.notna();ss=ss.loc[keep];mm=mm.loc[keep]
 for label,color in [(0,'#2c7fb8'),(1,'#d95f0e')]:
  yy=ss.loc[mm.poor==label,'fibrillar_I_III_V'];ax.scatter(label+rng.uniform(-.09,.09,len(yy)),yy,color=color,alpha=.75,s=25);ax.plot([label-.15,label+.15],[yy.mean()]*2,color='black',lw=2)
 row=st[(st.cohort==c)&(st.program=='fibrillar_I_III_V')&(st.model=='unadjusted_Welch')].iloc[0]
 ax.set_title(f'{c}: poor - good = {row.effect:.2f}\nWelch P = {row.p:.3g}');ax.set_xticks([0,1],[f'Good (n={int(row.n_good)})',f'Poor (n={int(row.n_poor)})']);ax.set_ylabel('Fibrillar I/III/V mean gene z score')
fig.savefig(O/'figures/fibrillar_response.png',dpi=180);fig.savefig(O/'figures/fibrillar_response.pdf');plt.close(fig)
print(st[st.program.isin(['fibrillar_I_III_V','collagen_IV','collagen_crosslinking'])].to_string(index=False),flush=True)
print(ge[ge.gene.isin(['PCOLCE','LOX','LOXL2','PXDN'])].to_string(index=False),flush=True)
(O/'COMPUTATION_COMPLETE.json').write_text(json.dumps({'exit_status':0,'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'program_rows':len(st),'gene_rows':len(ge)}))
