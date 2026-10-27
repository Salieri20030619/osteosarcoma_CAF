from pathlib import Path
import json,datetime,numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005';A=B/'phase64_response_followup_20261005'
cfg=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text())
sets={k:[g for g in gs if g not in cfg['drop_dominant'][k]] for k,gs in cfg['structure'].items()}
freeze={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sets':sets,'role':'Exploratory clinical sensitivity prompted by phase65 source/chain-composition results; phase64 full-program tests remain primary historical results. No new independent validation cohort, no replacement of original negative primary results.','methods':'Same phase64 cohort and patient response labels; submitted VST, median probes already mapped, gene z across all cohort samples; poor minus good Welch and age/sex (+baseline metastases GSE21257) OLS HC3. Three reduced panels form BH family within each cohort/model. No outcome-directed gene selection.'}
p=O/'resources/CLINICAL_REDUCED_PANEL_FOLLOWUP_FROZEN.json';assert not p.exists();p.write_text(json.dumps(freeze,indent=2))
meta=pd.read_csv(A/'tables/patient_metadata.tsv',sep='\t');records=[];scores=[];coverage=[]
for c in ['GSE39055','GSE21257']:
 m=meta[meta.cohort==c].set_index('gsm');e=pd.read_csv(A/'objects'/f'{c}_gene_VST.tsv.gz',sep='\t',index_col=0)
 z=e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1,ddof=1),axis=0)
 for name,gs in sets.items():
  have=[g for g in gs if g in z.index];assert len(have)>=3 and len(have)/len(gs)>=.8
  coverage.append({'cohort':c,'program':name,'measured':len(have),'expected':len(gs),'genes':';'.join(have)})
  v=z.loc[have].mean(axis=0).reindex(m.index);keep=m.poor.notna();a=v[keep&m.poor.eq(1)];b=v[keep&m.poor.eq(0)]
  scores.extend([{'cohort':c,'gsm':g,'program':name,'score':float(value)} for g,value in v.items()])
  va=a.var(ddof=1)/len(a);vb=b.var(ddof=1)/len(b);se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1));delta=a.mean()-b.mean();ci=stats.t.ppf(.975,df)*se
  records.append({'cohort':c,'program':name,'model':'unadjusted_Welch','N':int(keep.sum()),'n_poor':len(a),'n_good':len(b),'effect':delta,'ci_low':delta-ci,'ci_high':delta+ci,'P':stats.ttest_ind(a,b,equal_var=False).pvalue})
  cs=['poor','age_years','male']+(['metastases_at_diagnosis'] if c=='GSE21257' else [])
  x=sm.add_constant(m.loc[keep,cs].astype(float));assert np.linalg.matrix_rank(x)==x.shape[1]
  fit=sm.OLS(v[keep],x).fit(cov_type='HC3',use_t=True);ci=fit.conf_int().loc['poor']
  records.append({'cohort':c,'program':name,'model':'adjusted_HC3','N':int(keep.sum()),'n_poor':len(a),'n_good':len(b),'effect':fit.params['poor'],'ci_low':ci.iloc[0],'ci_high':ci.iloc[1],'P':fit.pvalues['poor']})
d=pd.DataFrame(records);d['family']=d.cohort+'|'+d.model+'|three_reduced_structure_panels';d['FDR']=d.groupby('family').P.transform(lambda x:multipletests(x,method='fdr_bh')[1])
d.to_csv(O/'tables/clinical_reduced_panel_statistics.tsv',sep='\t',index=False);pd.DataFrame(scores).to_csv(O/'tables/clinical_reduced_panel_scores.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/clinical_reduced_panel_coverage.tsv',sep='\t',index=False)
print(d.to_string(index=False),flush=True)
# Verify the serialized scores with an independent Welch summary computed from sample variance and scipy t survival.
ss=pd.read_csv(O/'tables/clinical_reduced_panel_scores.tsv',sep='\t',float_precision='round_trip')
for r in d[d.model=='unadjusted_Welch'].itertuples():
 q=ss[(ss.cohort==r.cohort)&(ss.program==r.program)].merge(meta[['gsm','poor']],on='gsm',validate='one_to_one').dropna(subset=['poor'])
 a=q.loc[q.poor==1,'score'].to_numpy();b=q.loc[q.poor==0,'score'].to_numpy();sa=np.var(a,ddof=1)/len(a);sb=np.var(b,ddof=1)/len(b)
 df=(sa+sb)**2/(sa*sa/(len(a)-1)+sb*sb/(len(b)-1));pval=2*stats.t.sf(abs(a.mean()-b.mean())/np.sqrt(sa+sb),df)
 assert np.isclose(pval,r.P) and len(q)==r.N
(O/'QA/CLINICAL_REDUCED_QA.json').write_text(json.dumps({'passed':True,'patient_level_Welch_recalculated':True,'models':len(d),'exploratory_followup':True},indent=2))
report=O/'delivery/CAF胶原结构与空间来源验证报告_20261005.txt';txt=report.read_text()
insert='\n3.6 来源分析后的缩减面板未解决真实应答关联缺失（探索性随访）\n\n'
for r in d[d.program=='fibrillar_I_III_V'].itertuples():insert+=f"剔除COL1A1/COL1A2，{r.cohort}，{r.model}：N={r.N}，低应答减高应答差值{r.effect:+.3f}，95% CI [{r.ci_low:+.3f}, {r.ci_high:+.3f}]，P={r.P:.6f}，本轮三面板家族FDR={r.FDR:.6f}。\n"
insert+='该补充在空间结果之后提出、检查临床缩减面板结果之前冻结，使用相同患者与既有标准化表达、相同调整协变量，对3个缩减面板在各队列/模型内BH校正；属于探索性敏感性分析，不是新独立临床验证。完整数据见clinical_reduced_panel_statistics.tsv。\n'
n=int((d.P<.05).sum());insert+=f"全部12项队列/面板/模型检验中，原始P<0.05共{n}项。原phase64主检验结果完整保留，不用新面板替换。\n"
assert '\n四、Discussion' in txt
report.write_text(txt.replace('\n四、Discussion',insert+'\n四、Discussion'),encoding='utf-8')
