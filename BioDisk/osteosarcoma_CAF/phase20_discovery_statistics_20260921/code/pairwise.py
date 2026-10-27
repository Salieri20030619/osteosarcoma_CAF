from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import wilcoxon,binomtest
from statsmodels.stats.multitest import multipletests
B=Path('C:/codex_test');I=B/'CAF_evidence_20260921';O=B/'CAF_discovery_stats_20260921'
p=json.loads((O/'PROTOCOL.json').read_text());p['additional_direct_comparison']='supported C7 minus C9, all candidate cells and myCAF-positive-only; >=10 cells per cluster per patient; same3 recruitment cohorts; separate5-program BH family per population. Exploratory focused comparison added to distinguish the two frozen states.'
(O/'PROTOCOL.json').write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
d=pd.read_csv(I/'tables/supported_analysis_labels.tsv.gz',sep='\t');d=d[d.dataset.isin(p['primary_cohorts'])&d.known_patient]
progs=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer'];rows=[];results=[]
for population in ['all_candidate_cells','within_myCAF']:
 dd=d if population=='all_candidate_cells' else d[d.myCAF_customer_positive]
 for donor,g in dd.groupby('analysis_donor'):
  a=g[g.cluster==7];b=g[g.cluster==9]
  for pr in progs:
   x=a[pr+'_score'].dropna();y=b[pr+'_score'].dropna()
   if len(x)>=10 and len(y)>=10:rows.append({'population':population,'donor_record':donor,'dataset':g.dataset.iloc[0],'program':pr,'C7_cells':len(x),'C9_cells':len(y),'delta_C7_minus_C9':x.mean()-y.mean()})
 e=pd.DataFrame([r for r in rows if r['population']==population])
 for pr in progs:
  g=e[e.program==pr];v=np.round(g.delta_C7_minus_C9.to_numpy(),12);n=len(v)
  assert np.unique(abs(v)).size==n and (v!=0).all()
  results.append({'population':population,'program':pr,'n_patients':n,'C7_higher_patients':int((v>0).sum()),'C9_higher_patients':int((v<0).sum()),'median_C7_minus_C9':np.median(v) if n else np.nan,'P':float(wilcoxon(v,method='exact').pvalue) if n>=2 else np.nan,'P_sign':binomtest(int((v>0).sum()),n).pvalue if n>=2 else np.nan})
r=pd.DataFrame(results);r['FDR']=np.nan
for _,ix in r.groupby('population').groups.items():r.loc[ix,'FDR']=multipletests(r.loc[ix,'P'].fillna(1),method='fdr_bh')[1]
r.to_csv(O/'tables/C7_vs_C9_direct_tests.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(O/'tables/C7_vs_C9_patient_differences.tsv',sep='\t',index=False)
print(r.to_string(index=False))
