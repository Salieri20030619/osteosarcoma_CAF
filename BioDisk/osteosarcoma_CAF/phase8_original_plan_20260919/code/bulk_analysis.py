from pathlib import Path
import gzip,csv,io,re,json
import numpy as np,pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase8_original_plan_20260919')
assert (O/'bulk/DOWNLOAD_COMPLETE').exists()
with gzip.open(O/'bulk/GSE21257_series_matrix.txt.gz','rt') as f:lines=f.readlines()
meta=[];start=end=None
for i,line in enumerate(lines):
 if line.startswith('!series_matrix_table_begin'):start=i+1
 elif line.startswith('!series_matrix_table_end'):end=i
 elif line.startswith('!Sample_'):meta.append(next(csv.reader([line],delimiter='\t')))
assert start and end
expr=pd.read_csv(io.StringIO(''.join(lines[start:end])),sep='\t',index_col=0)
assert expr.shape[1]==53 and expr.columns.is_unique and expr.index.is_unique
expr=expr.astype(float);assert np.isfinite(expr.to_numpy()).all()
cols=next(x[1:] for x in meta if x[0]=='!Sample_geo_accession')
assert cols==expr.columns.tolist()
c=pd.DataFrame(index=cols);c.index.name='sample'
for row in meta:
 if row[0] in ['!Sample_title','!Sample_source_name_ch1']:c[row[0][8:]]=row[1:]
 if row[0].startswith('!Sample_characteristics'):
  for sample,value in zip(cols,row[1:]):
   key,val=value.split(':',1);c.loc[sample,key.strip()]=val.strip()
assert c['group'].notna().all()
c['metastasis']=~c['group'].str.lower().eq('no metastases')
c['baseline_metastasis']=c['group'].str.contains('present at diagnosis',case=False)
c['age_years']=c['age'].str.extract(r'([0-9.]+)').astype(float)/12
c['male']=c['gender'].map({'M':1,'F':0})
c['followup_months']=c['status'].str.extract(r'([0-9.]+)').astype(float)
c['metastasis_months']=c['group'].str.extract(r'at ([0-9.]+) months').astype(float)
c.loc[c.baseline_metastasis,'metastasis_months']=0
assert int(c.metastasis.sum())==34 and (~c.metastasis).sum()==19
c['short_recorded_followup_nonmetastatic']=(~c.metastasis)&(c.followup_months<60)
c.to_csv(O/'tables/GSE21257_clinical_mapping.tsv',sep='\t')
# GEO states lumi variance-stabilizing transformation and robust spline normalization.
# Use these published normalized values directly; no second log or count normalization.
ann=pd.read_csv(O/'bulk/GPL10295_platform_table.tsv',sep='\t',dtype=str).set_index('ID')
mapping=ann['Symbol'].reindex(expr.index).fillna('').str.strip()
valid=(mapping!='')&~mapping.str.contains(r'///|;|\|',regex=True)
gene=expr.loc[valid].groupby(mapping[valid]).median()
gene.to_csv(O/'bulk/GSE21257_gene_median_expression.tsv.gz',sep='\t')
panels=json.loads((O/'resources/programs_frozen.json').read_text())
family=json.loads((O/'PROTOCOL_FROZEN.json').read_text())['bulk_primary_family']
scores=pd.DataFrame(index=cols);coverage=[];out=[];sens=[]
for name in family:
 gs=panels[name];present=[g for g in gs if g in gene.index and gene.loc[g].std(ddof=1)>0]
 coverage.append({'program':name,'requested':len(gs),'covered':len(present),'missing':';'.join(g for g in gs if g not in present)})
 assert len(present)>=max(3,len(gs)*.5),(name,len(present))
 mat=gene.loc[present]
 zz=mat.sub(mat.mean(1),axis=0).div(mat.std(1,ddof=1),axis=0)
 score=zz.mean(0);score=(score-score.mean())/score.std(ddof=1);scores[name]=score
 def model(mask,adjust=False):
  df=pd.DataFrame({'score':score,'y':c.metastasis.astype(int),'age':c.age_years,'male':c.male}).loc[mask]
  keys=['score']+(['age','male'] if adjust else [])
  df=df.dropna(subset=keys+['y'])
  fit=sm.Logit(df.y,sm.add_constant(df[keys])).fit(disp=False,maxiter=100)
  assert fit.mle_retvals['converged']
  b=float(fit.params['score']);ci=fit.conf_int().loc['score'].to_numpy()
  return {'n':len(df),'events':int(df.y.sum()),'log_OR':b,'OR_per_SD':float(np.exp(b)),'CI_low':float(np.exp(ci[0])),'CI_high':float(np.exp(ci[1])),'P_value':float(fit.pvalues['score'])}
 main=model(np.ones(len(c),bool));main['program']=name;out.append(main)
 for variant,mask,adjust in [('age_sex_adjusted',np.ones(len(c),bool),True),('exclude_baseline_metastasis',~c.baseline_metastasis,False),('exclude_short_nonmet_followup',~c.short_recorded_followup_nonmetastatic,False)]:
  r=model(mask,adjust);r.pop('P_value');r.update(program=name,variant=variant);sens.append(r)
 loo=[]
 for sample in cols:
  r=model(c.index!=sample);loo.append(r['log_OR'])
 sens.append({'program':name,'variant':'leave_one_patient_out','n':52,'log_OR_min':min(loo),'log_OR_max':max(loo),'same_direction_fraction':float((np.sign(loo)==np.sign(main['log_OR'])).mean())})
out=pd.DataFrame(out);out['FDR_BH']=multipletests(out.P_value,method='fdr_bh')[1]
out['FDR_family']='8_frozen_bulk_program_endpoints';out.to_csv(O/'tables/GSE21257_primary_associations.tsv',sep='\t',index=False)
pd.DataFrame(sens).to_csv(O/'tables/GSE21257_sensitivity_effects.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'tables/GSE21257_signature_coverage.tsv',sep='\t',index=False)
scores.to_csv(O/'tables/GSE21257_program_scores.tsv',sep='\t')
summary={'samples':53,'metastasis':34,'no_metastasis':19,'baseline_metastasis':int(c.baseline_metastasis.sum()),'nonmetastatic_short_recorded_followup':int(c.short_recorded_followup_nonmetastatic.sum()),'probe_rows':len(expr),'gene_rows':len(gene),'probe_annotation_matched':int(mapping.ne('').sum()),'input_value_quantiles':np.quantile(expr.to_numpy(),[0,.5,.99,1]).tolist(),'normalization':'author lumi VST and robust spline normalized; no additional log transform','endpoint':'original reported metastasis group, not lung-specific; baseline metastasis included in primary and excluded in sensitivity','FDR_significant_programs':out.loc[out.FDR_BH<.05,'program'].tolist(),'interpretation':'bulk tissue association, not CAF abundance or a clinically validated predictor'}
(O/'bulk/analysis_summary.json').write_text(json.dumps(summary,indent=2))
(O/'BULK_COMPLETE').write_text('BULK_ASSOCIATION_AND_SENSITIVITY_COMPLETE\n')
print(json.dumps(summary),flush=True);print(out.to_string(index=False))
