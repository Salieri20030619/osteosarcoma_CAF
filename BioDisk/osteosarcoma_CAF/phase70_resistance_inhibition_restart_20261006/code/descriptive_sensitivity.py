from pathlib import Path
import pandas as pd,numpy as np,statsmodels.api as sm,json
b=Path(__file__).resolve().parents[1];i=b/'inputs'
assert (b/'resources/SENSITIVITY_PROTOCOL.json').exists()
meta=pd.read_csv(i/'patient_metadata.tsv',sep='\t').set_index('gsm');records=[]
for co in ['GSE39055','GSE21257']:
 e=pd.read_csv(i/f'{co}_gene_VST.tsv.gz',sep='\t',index_col=0)
 m=meta[(meta.cohort==co)&meta.included_primary.eq(True)].copy()
 cov=['poor','age_years','male']+(['metastases_at_diagnosis'] if co=='GSE21257' else [])
 m=m.dropna(subset=cov)
 groups=[('leave_one_patient_out',gsm,m.drop(index=gsm)) for gsm in m.index]
 if co=='GSE39055':groups.append(('MAP_only','two_AP_patients',m[m.chemotherapy.str.startswith('MAP')]))
 for kind,excluded,sub in groups:
  x=sm.add_constant(sub[cov].astype(float));assert np.linalg.matrix_rank(x)==x.shape[1]
  for g in ['ACAN','COL6A2']:
   f=sm.OLS(e.loc[g,sub.index].astype(float),x).fit()
   records.append({'cohort':co,'gene':g,'sensitivity':kind,'excluded':excluded,'n':len(sub),'n_poor':int(sub.poor.sum()),'n_good':int((sub.poor==0).sum()),'effect':f.params['poor']})
pd.DataFrame(records).to_csv(b/'tables/positive_candidate_descriptive_sensitivity.tsv',sep='\t',index=False)
