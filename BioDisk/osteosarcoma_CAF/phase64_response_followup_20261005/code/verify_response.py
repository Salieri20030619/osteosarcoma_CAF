from pathlib import Path
import json,numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
meta=pd.read_csv(O/'tables/patient_metadata.tsv',sep='\t');scores=pd.read_csv(O/'tables/patient_program_scores.tsv',sep='\t');st=pd.read_csv(O/'tables/program_response_statistics.tsv',sep='\t');ge=pd.read_csv(O/'tables/fixed_program_gene_response_statistics.tsv',sep='\t');cov=pd.read_csv(O/'tables/program_coverage.tsv',sep='\t')
assert len(meta)==90 and meta.gsm.is_unique and scores.gsm.is_unique
assert len(meta.dropna(subset=['poor']))==83
checks=[]
for c in ['GSE39055','GSE21257']:
 m=meta[meta.cohort==c].set_index('gsm');s=scores[scores.cohort==c].set_index('gsm')
 assert set(s.index)==set(m.index)
 if c=='GSE39055':
  assert len(m)==37 and m.poor.isna().sum()==1
  assert m.loc[m.poor.isna(),'response_raw'].iloc[0]=='>80% necrosis'
 else:
  assert len(m)==53 and m.poor.isna().sum()==6 and sum(m.poor==1)==29 and sum(m.poor==0)==18
  assert (m.age_years*12-m.age.str.extract(r'(\d+)')[0].astype(float)).abs().max()<1e-8
 expr=pd.read_csv(O/'objects'/f'{c}_gene_VST.tsv.gz',sep='\t',index_col=0)
 primary=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL5A3']
 observed=[x for x in primary if x in expr.index]
 z=expr.loc[observed].sub(expr.loc[observed].mean(axis=1),axis=0).div(expr.loc[observed].std(axis=1),axis=0).mean(axis=0)
 assert np.allclose(z.loc[s.index],s.fibrillar_I_III_V)
 known=m.poor.notna();ids=m.index[known];labels=m.loc[ids,'poor'];y=s.loc[ids,'fibrillar_I_III_V']
 row=st[(st.cohort==c)&(st.program=='fibrillar_I_III_V')&(st.model=='unadjusted_Welch')].iloc[0]
 a=y[labels==1];b=y[labels==0]
 assert np.isclose(row.effect,a.mean()-b.mean())
 assert np.isclose(row.p,stats.ttest_ind(a,b,equal_var=False).pvalue)
 for model in ['unadjusted_Welch','adjusted_HC3']:
  d=st[(st.cohort==c)&(st.model==model)]
  assert np.allclose(d.fdr_all14,multipletests(d.p,method='fdr_bh')[1])
  e=d[d.family=='secondary13'];assert np.allclose(e.fdr_secondary13,multipletests(e.p,method='fdr_bh')[1])
 d=ge[ge.cohort==c];assert np.allclose(d.fdr_gene_family,multipletests(d.p,method='fdr_bh')[1])
 checks.append({'cohort':c,'sample_join':'passed','response_label':'passed','primary_score_recalculated':'passed','primary_effect_and_P_recalculated':'passed','program_and_gene_FDR':'passed','n':len(ids)})
mapping=pd.read_csv(O/'tables/GSE21257_probe_mapping.tsv',sep='\t').set_index('ID')
assert mapping.loc['9ilipx9ZCOSmZZ9Vhw','canonical_symbol']=='PLPPR3'
assert mapping.loc['NSjRKdq2eSGf0ur4aQ','canonical_symbol']=='PRG2'
assert pd.isna(mapping.loc['iVWUrXuSKJXuieoXkI','canonical_symbol'])
assert mapping.loc['Bbg0g0HnkhSis7iip4','canonical_symbol']=='ACAN'
assert cov[cov.program=='fibrillar_I_III_V'].evaluable.all()
assert np.isfinite(st[['effect','p']]).all().all()
assert ((st.p>=0)&(st.p<=1)).all()
(O/'QA/ACCEPTANCE.json').write_text(json.dumps({'passed':True,'checks':checks,'scope':'patient-level bulk response analysis; spatial branch source audit only'},indent=2))
print('ACCEPTANCE_PASSED',checks)
print('COVERAGE',cov[['cohort','program','n_measured','n_expected','evaluable']].to_string(index=False))
print('ALL_PROGRAMS',st[st.model=='unadjusted_Welch'][['cohort','program','effect','p','fdr_secondary13','fdr_all14']].to_string(index=False))
