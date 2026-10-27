# SYNTAX_CHECK: python -m py_compile clinical_small_sample_sensitivity.py
# RUN_COMMAND: batch_bench/bin/python -u clinical_small_sample_sensitivity.py
# EXPECTED_SUCCESS_MARKER: CLINICAL_SMALL_SAMPLE_SENSITIVITY_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,ast,itertools
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import norm,t
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');D=P/'phase44_patient_mechanism_20261001';O=D/'clinical_model_sensitivity'
for sub in ['tables','resources','delivery','logs']:(O/sub).mkdir(parents=True,exist_ok=True)
protocol={'at':datetime.datetime.now().astimezone().isoformat(),'trigger':'review of primary HC3 normal-reference signals at N=6-7','scope':'all evaluable primary adjusted models, not only significant ones','primary_results':'retained unchanged','sensitivity':'same patients/covariates/coefficient/HC3 SE, Student t reference with actual residual df','target_effect':'target_group, no substitution by covariate P','multiplicity':'same complete predeclared adjusted family including unavailable rows, suffix HC3_t_reference_sensitivity','documentation':'https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.RegressionResults.get_robustcov_results.html'}
(O/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2));out=[];rc=1
try:
 names={'strict_Bone_vs_Lung':('Bone','Lung'),'reported_nonmetastatic_vs_metastatic_bone':('Reported_nonmetastatic','Reported_metastatic_bone'),'history_supported_B1_vs_strict_B2':('B1_history_timing_unresolved','B2_strict_lung')}
 for root in [D,D/'clinical_reporting_extension']:
  primary=pd.read_csv(root/'tables/adjusted_clinical_fraction_statistics.tsv',sep='\t');patients=pd.read_csv(root/'tables/patient_event_class_state_fractions.tsv.gz',sep='\t',keep_default_na=False);patients['fraction']=pd.to_numeric(patients.fraction,errors='coerce')
  for row in primary.itertuples(index=False):
   record=row._asdict();record['source_scope']=root.name;record['primary_P']=record.pop('P');record['primary_FDR']=record.pop('FDR');record['P']=np.nan;record['df_resid']=np.nan;record['normal_reference_reconstruction_error']=np.nan;record['family']=record['family']+'|HC3_t_reference_sensitivity'
   if not np.isfinite(row.P):record['sensitivity_status']='primary_not_evaluable';out.append(record);continue
   z=patients[(patients.layer==row.layer)&(patients.dataset==row.dataset)&(patients.contrast==row.contrast)&(patients.feature==row.feature)&(patients.denominator==row.denominator)&patients.fraction.notna()];first,second=names[row.contrast];a=z[z.group.eq(first)].set_index('donor_record');b=z[z.group.eq(second)].set_index('donor_record');paired=set(a.index)&set(b.index);a=a.drop(index=list(paired));b=b.drop(index=list(paired));model=pd.concat([a.assign(target_group=0),b.assign(target_group=1)])
   X=pd.concat([model[['target_group']].astype(float),pd.get_dummies(model[['treatment','technology']],drop_first=True,dtype=float)],axis=1);X=sm.add_constant(X,has_constant='add');fit=sm.OLS(model.fraction.astype(float),X).fit(cov_type='HC3',use_t=True)
   coef=float(fit.params['target_group']);se=float(fit.bse['target_group']);pz=2*norm.sf(abs(coef/se));pt=float(fit.pvalues['target_group']);manual=2*t.sf(abs(coef/se),fit.df_resid)
   assert len(model)==row.N and abs(coef-row.target_effect)<1e-10 and abs(se-row.SE_HC3)<1e-10 and abs(pz-row.P)<1e-10 and abs(pt-manual)<1e-12
   record.update(P=pt,df_resid=fit.df_resid,sensitivity_status='HC3_t_reference_reconstructed',normal_reference_reconstruction_error=abs(pz-row.P));out.append(record)
 df=pd.DataFrame(out);df['FDR']=np.nan
 for family,g in df.groupby('family'):
  q=multipletests(g.P.fillna(1),method='fdr_bh')[1];df.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 df['sensitivity_nominal']=df.P.lt(.05);df['sensitivity_FDR_significant']=df.FDR.lt(.05);df.to_csv(O/'tables/all_adjusted_HC3_t_reference_sensitivity.tsv',sep='\t',index=False)
 print('PRIMARY_NOMINAL_FOLLOWUP',df[df.primary_P.lt(.05)][['layer','dataset','contrast','feature','N','df_resid','target_effect','primary_P','primary_FDR','P','FDR']].to_string(index=False))
 (O/'delivery/SENSITIVITY_QA.json').write_text(json.dumps({'rows':len(df),'primary_evaluable':int(df.primary_P.notna().sum()),'primary_nominal':int(df.primary_P.lt(.05).sum()),'sensitivity_nominal':int(df.P.lt(.05).sum()),'sensitivity_FDR':int(df.FDR.lt(.05).sum()),'coefficients_SE_and_normal_P_reconstructed':True,'primary_results_unchanged':True},indent=2));(O/'COMPLETE').write_text('SUCCESS\n');print('CLINICAL_SMALL_SAMPLE_SENSITIVITY_COMPLETE');rc=0
except Exception:traceback.print_exc();(O/'FAILED').write_text(traceback.format_exc())
finally:(O/'sensitivity.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)