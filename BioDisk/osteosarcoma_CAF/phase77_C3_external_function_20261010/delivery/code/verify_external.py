from pathlib import Path
import pandas as pd,numpy as np,json,itertools,shutil
from scipy.stats import spearmanr
B=Path('C:/Users/31937/CAF_C3_external_function_20261010');R=B/'results';T=R/'tables';D=B/'delivery';D.mkdir(exist_ok=True)
cs=pd.read_csv(T/'query_cell_scores.tsv.gz',sep='\t');stats=pd.read_csv(T/'functional_statistics.tsv',sep='\t');ef=pd.read_csv(T/'patient_effects.tsv',sep='\t');cr=pd.read_csv(T/'continuous_patient_correlations.tsv',sep='\t');cstats=pd.read_csv(T/'continuous_statistics.tsv',sep='\t');proto=json.loads((R/'resources/PROTOCOL_FROZEN.json').read_text());coef=np.array(json.loads((R/'resources/frozen_residual_model.json').read_text())['coefficients']);features=pd.read_csv(R/'resources/frozen_mapping_model.tsv',sep='\t');ex=pd.read_csv(R/'resources/mapping_excluded_genes.tsv',sep='\t');assert not set(features.gene)&set(ex.gene);assert np.allclose(cs.processing_residual,cs.processing2-coef[0]-coef[1]*cs.published_COL5_4-coef[2]*np.log1p(cs.library));assert cs.cell_id.is_unique
records=0
for _,r in stats.iterrows():
 q=cs[(cs.dataset==r.dataset)&(cs.tissue==r.tissue)];vals=[]
 for donor,dd in q.groupby('donor_record'):
  vv=[]
  for site,z in dd.groupby('site'):
   l=z[z.C3like];rt=z[~z.C3like]
   if min(len(l),len(rt))>=r.min_cells and l[r.panel].notna().all() and rt[r.panel].notna().all():vv.append(l[r.panel].mean()-rt[r.panel].mean())
  if vv:vals.append(float(np.mean(vv)))
 assert len(vals)==r.N and sum(v>0 for v in vals)==r.positive
 assert np.isclose(np.mean(vals),r.mean_effect) if vals else pd.isna(r.mean_effect)
 if len(vals)>=2 and np.any(vals):
  p=np.mean([abs(np.mean(np.array(v)*vals))>=abs(np.mean(vals))-1e-12 for v in itertools.product([-1,1],repeat=len(vals))]);assert np.isclose(p,r.P)
 else:assert pd.isna(r.P)
 records+=1
for _,r in cr.iterrows():
 q=cs[(cs.dataset==r.dataset)&(cs.tissue==r.tissue)&(cs.donor_record==r.donor)];v=[]
 for site,z in q.groupby('site'):
  if len(z)>=5 and z.mapping_logit.nunique()>1 and z.processing_residual.nunique()>1:v.append(spearmanr(z.mapping_logit,z.processing_residual).statistic)
 assert np.isclose(np.mean(v),r.rho)
for _,r in cstats.iterrows():
 q=cr[(cr.dataset==r.dataset)&(cr.tissue==r.tissue)];v=q.rho.to_numpy();assert len(v)==r.N and sum(v>0)==r.positive
 assert np.isclose(np.mean(v),r.mean_rho) if len(v) else pd.isna(r.mean_rho)
 if len(v)>=2 and np.any(v):assert np.isclose(np.mean([abs(np.mean(np.array(sign)*v))>=abs(np.mean(v))-1e-12 for sign in itertools.product([-1,1],repeat=len(v))]),r.P)
 else:assert pd.isna(r.P)
 records+=1
for d in [stats,cstats]:
 for family,g in d[d.P.notna()].groupby('family'):
  x=g.P.to_numpy();o=np.argsort(x);expected=np.minimum(1,np.minimum.accumulate((x[o]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]);assert np.allclose(expected,g.FDR.to_numpy()[o])
ref=pd.read_csv(T/'reference_heldout_predictions.tsv.gz',sep='\t');assert ref.cell_id.is_unique and len(ref)==1456 and ref.true_C3.sum()==34;tp=sum(ref.true_C3 & ref.predicted_C3like);fp=sum(~ref.true_C3 & ref.predicted_C3like);fn=sum(ref.true_C3 & ~ref.predicted_C3like);tn=sum(~ref.true_C3 & ~ref.predicted_C3like);q=json.loads((R/'resources/MAPPING_QUALITY.json').read_text());assert [tn,fp,fn,tp]==q['confusion_TN_FP_FN_TP'];assert np.isclose(tp/(tp+fp),q['pooled_precision']) and np.isclose(tp/(tp+fn),q['pooled_C3_recall']);assert (R/'logs/analysis.exit_code').read_text().strip()=='0';assert json.loads((R/'COMPLETE.json').read_text())['status']=='completed'
qa={'status':'passed','analysis_exit_code':0,'all_81_statistics_recomputed':records==81,'query_cells':len(cs),'cell_ids_unique':True,'heldout_mapping_confusion_verified':[tn,fp,fn,tp],'mapping_features':len(features),'no_endpoint_or_matrisome_feature_leak':True,'residual_formula_reproduced':True,'all_new_BH_families_recomputed':True,'external_query_patients_primary':int(cs[(cs.dataset=='GSE237070')&(cs.tissue=='primary')].donor_record.nunique()),'external_query_patients_adjacent':int(cs[(cs.dataset=='GSE237070')&(cs.tissue=='adjacent')].donor_record.nunique()),'scope':'Statistics, file integrity and mapping-exclusion QA; not biological validation of transferred cell identity'}
(D/'QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa));print('MIN5');print(stats[(stats.dataset=='GSE237070')&(stats.min_cells==5)].drop(columns='family').to_csv(index=False))
