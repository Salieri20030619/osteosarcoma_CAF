import sys
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase46_spatial_functional_20261001/code')
from common import *
import ast,math,itertools,anndata as ad,statsmodels.api as sm
from scipy import sparse
Q=P/'phase48_clinical_programs_delivery_20261001';D=P/'phase44_patient_mechanism_20261001';pr=json.loads((O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json').read_text());programs=pr['programs'];source=D/'code/stageD_clinical.py'
# Recover the already accepted independent-patient permutation function only.
source=D/'code/stageD_clinical.py';nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='independent'];assert len(nodes)==1;rng=np.random.default_rng(20261001);exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'))
m=pd.read_csv(D/'clinical_site_review/tables/clinical_manifest_with_strict_site_group.tsv',sep='\t',keep_default_na=False).set_index('sample');a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad');gm=feature_map(a);presence=sparse.load_npz(A/'objects/sample_feature_presence.npz');si=pd.read_csv(A/'tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample').presence_row.to_dict();gi=pd.read_csv(A/'tables/feature_presence_gene_index.tsv.gz',sep='\t').set_index('feature').presence_column.to_dict();assert np.allclose(np.asarray(a.X.sum(1)).ravel(),a.obs.full_library_counts.to_numpy(float))
comp=[('strict_Bone_vs_Lung','strict_site_group','Bone','Lung'),('reported_nonmetastatic_vs_metastatic','reported_status_group','Reported_nonmetastatic','Reported_metastatic_bone'),('history_B1_vs_B2','history_group','B1_history_timing_unresolved','B2_strict_lung')]
m['reported_status_group']=m.clinical_group.map({'A_candidate_reported_no_metastasis_timing_unknown':'Reported_nonmetastatic','A_reported_nonmetastasis':'Reported_nonmetastatic','metastatic_bone_primary_lung_and_timing_unconfirmed':'Reported_metastatic_bone','Metastatic_primary_lung_unconfirmed':'Reported_metastatic_bone'}).fillna('Unknown');m['history_group']='Unknown';m.loc[m.clinical_group.eq('B1_history_supported_timing_unresolved')&m.strict_site_group.eq('Bone'),'history_group']='B1_history_timing_unresolved';m.loc[m.clinical_group.eq('B2')&m.strict_site_group.eq('Lung'),'history_group']='B2_strict_lung'
dump(Q/'resources/CLINICAL_PROGRAM_PROTOCOL_FROZEN.json',{'at':datetime.datetime.now().astimezone().isoformat(),'signatures_SHA':sha(O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json'),'programs':list(programs),'classes':['F2','F3','F4'],'unit':'within source patient clinical event group/class counts summed across samples; full-source library pseudobulk log2(CP10k+.1), mean over recorded program genes','comparisons':comp,'paired_removed_from_independent':True,'models':'independent raw label permutation plus optional target group+treatment+technology HC3 with t reference; complete recorded covariates, full rank residual df>=2','family':'layer/source/class/contrast/design across all15 programs; adjusted separate','Unknown_not_No':True,'unavailable_class_not_zero':True,'allzero_program_P_NA':True,'no_source_pooling':True})
values=[];coverage=[]
for layer in ['supported','expanded']:
 labs=a.obs[layer+'_class_id'].astype(str)
 for contrast,field,g0,g1 in comp:
  for group in [g0,g1]:
   samples=m.index[m[field].eq(group)];idx=np.flatnonzero(a.obs['sample'].isin(samples));sub=a.obs.iloc[idx]
   for (ds,donor),t in sub.groupby(['dataset','donor_record'],observed=True):
    for cl in ['F2','F3','F4']:
     ids=idx[sub.index.get_indexer(t.index[labs.loc[t.index].eq(cl)])]
     if not len(ids):continue
     sms=sorted(set(a.obs.iloc[ids]['sample'].astype(str)));den=float(a.obs.iloc[ids].full_library_counts.sum());raw=np.asarray(a.X[ids].sum(0)).ravel();lg=np.log2(raw/den*1e4+.1) if den>0 else np.full(a.n_vars,np.nan);treat=sorted(set(m.loc[sms,'treatment']));tech=sorted(set(m.loc[sms,'technology']));treat=';'.join(treat) if len(treat)==1 and all(x.lower() not in ['','field_absent','na','unknown','not_reported','unavailable'] and ';' not in x for x in treat) else '';tech=';'.join(tech) if len(tech)==1 and all(x.lower() not in ['','field_absent','na','unknown','not_reported','unavailable'] and ';' not in x for x in tech) else ''
     for name,genes in programs.items():
      have=[g for g in genes if g in gm and all(bool(presence[si[sm],gi[str(a.var_names[gm[g]])]]) for sm in sms)];v=np.mean([lg[gm[g]] for g in have]) if have else np.nan
      values.append(dict(layer=layer,dataset=ds,patient_id=donor,class_id=cl,contrast=contrast,group=group,endpoint=name,value=v,cells=len(ids),library_UMI=den,treatment=treat,technology=tech,available_genes=len(have)))
      coverage.append(dict(layer=layer,dataset=ds,patient_id=donor,class_id=cl,contrast=contrast,group=group,endpoint=name,requested=len(genes),present=len(have),missing=';'.join(sorted(set(genes)-set(have)))))
f=pd.DataFrame(values);f.to_csv(Q/'tables/clinical_patient_program_values.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(Q/'tables/clinical_program_coverage.tsv',sep='\t',index=False);tests=[];roster=[];modelrows=[];paired=[]
def pairedtest(d):
 d=np.asarray(d,float);n=len(d)
 if n<2 or np.all(d==0):return np.nan,'unassessable'
 assert n<=20;s=np.asarray(list(itertools.product([-1.,1.],repeat=n)));return float(np.mean(np.abs(s@d/n)>=abs(d.mean())-1e-12)),'exact_patient_signflip'
for key,t in f.groupby(['layer','dataset','class_id','contrast','endpoint'],observed=True):
 layer,ds,cl,contrast,endpoint=key;g0,g1=next((c[2],c[3]) for c in comp if c[0]==contrast);base=dict(zip(['layer','dataset','class_id','contrast','endpoint'],key));v=t[np.isfinite(t.value)];both=set(v.loc[v.group.eq(g0),'patient_id'])&set(v.loc[v.group.eq(g1),'patient_id']);pairs=v[v.patient_id.isin(both)].pivot(index='patient_id',columns='group',values='value');d=pairs[g1]-pairs[g0] if len(pairs) else pd.Series(dtype=float)
 for donor,value in d.items():paired.append({**base,'patient_id':donor,'effect':value})
 p,method=pairedtest(d.to_numpy());tests.append({**base,'design':'paired','N0':len(d),'N1':len(d),'N':len(d),'mean_effect':d.mean() if len(d) else np.nan,'P':p,'test':method,'status':'evaluable' if np.isfinite(p) else 'insufficient_or_constant'})
 single=v[~v.patient_id.isin(both)];a0=single.loc[single.group.eq(g0),'value'].to_numpy();b1=single.loc[single.group.eq(g1),'value'].to_numpy();p,method=independent(a0,b1);tests.append({**base,'design':'independent','N0':len(a0),'N1':len(b1),'N':len(a0)+len(b1),'mean_effect':np.mean(b1)-np.mean(a0) if len(a0) and len(b1) else np.nan,'P':p,'test':method,'status':'evaluable' if np.isfinite(p) else 'insufficient_or_constant'})
 for r in single.itertuples(index=False):roster.append({**base,'patient_id':r.patient_id,'group0_or_1':int(r.group==g1),'value':r.value,'treatment':r.treatment,'technology':r.technology})
 modelp=np.nan;coef=np.nan;reason='insufficient_or_unrecorded_covariates';nmodel=0;dfres=np.nan;se=np.nan;mt=single[single.treatment.ne('')&single.technology.ne('')].copy();mt['target_group']=(mt.group==g1).astype(int)
 if len(a0)>=2 and len(b1)>=2 and mt.target_group.nunique()==2 and mt.value.nunique()>1:
  design=pd.concat([pd.Series(1.,index=mt.index,name='intercept'),mt.target_group.astype(float),pd.get_dummies(mt[['treatment','technology']],drop_first=True,dtype=float)],axis=1);xx=design.to_numpy(float);nmodel=len(mt);rank=np.linalg.matrix_rank(xx);dfres=nmodel-rank
  leverage=np.diag(xx@np.linalg.pinv(xx));valid_hc3=bool(np.all(leverage<1-1e-12))
  if rank==xx.shape[1] and dfres>=2 and valid_hc3:
   fit=sm.OLS(mt.value.to_numpy(float),xx).fit(cov_type='HC3',use_t=True);coef=float(fit.params[1]);modelp=float(fit.pvalues[1]);se=float(fit.bse[1]);reason='evaluable'
   mid='|'.join(map(str,key))
   for pos,(ix,r) in enumerate(mt.iterrows()):modelrows.append({'model_id':mid,'patient_id':r.patient_id,'value':r.value,**{f'X{j}':xx[pos,j] for j in range(xx.shape[1])}})
  else:reason='HC3_undefined_unit_leverage' if not valid_hc3 else 'rank_or_residual_df'
 tests.append({**base,'design':'adjusted_HC3_t','N0':int(mt.target_group.eq(0).sum()),'N1':int(mt.target_group.eq(1).sum()),'N':nmodel,'mean_effect':coef,'P':modelp,'test':'target_group_HC3_t','status':reason,'residual_df':dfres,'HC3_SE':se})
s=pd.DataFrame(tests);s['family']=s.layer+'|'+s.dataset+'|'+s.class_id+'|'+s.contrast+'|'+s.design+'|clinical_programs';s['FDR']=s.groupby('family').P.transform(bh);s['nominal']=s.P.lt(.05);s.to_csv(Q/'tables/clinical_program_statistics.tsv',sep='\t',index=False);pd.DataFrame(roster).to_csv(Q/'tables/independent_patient_roster.tsv',sep='\t',index=False);pd.DataFrame(paired,columns=['layer','dataset','class_id','contrast','endpoint','patient_id','effect']).to_csv(Q/'tables/paired_patient_program_deltas.tsv',sep='\t',index=False);(pd.DataFrame(modelrows) if modelrows else pd.DataFrame(columns=['model_id','patient_id','value','X0','X1'])).to_csv(Q/'tables/HC3_model_designs.tsv',sep='\t',index=False)
qa={'patient_program_records':len(f),'statistics':len(s),'valid':int(s.P.notna().sum()),'nominal':int(s.nominal.sum()),'FDR':int(s.FDR.lt(.05).sum()),'adjusted_fitted':int(s.loc[s.design.eq('adjusted_HC3_t'),'P'].notna().sum()),'no_cross_source_pooling':True,'class_missing_not_zero':True};dump(Q/'delivery/CLINICAL_PROGRAM_QA.json',qa);(Q/'CLINICAL_PROGRAM_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
