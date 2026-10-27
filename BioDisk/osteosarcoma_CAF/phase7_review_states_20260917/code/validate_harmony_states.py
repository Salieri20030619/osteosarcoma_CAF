from pathlib import Path
from itertools import product
import json
import numpy as np,pandas as pd,scanpy as sc
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();out=p/'integration_review/validation';out.mkdir(exist_ok=True)
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');key='harmony_theta2_r0.5';labels=a.obs[key].astype(str).to_numpy();patients=a.obs.patient_id.astype(str).to_numpy();counts=a.layers['counts'];states=sorted(set(labels),key=int);programs=[c for c in a.obs if c.startswith('program_')]
markers=[];pe=[];pg=[];stability=[];audit=[]
for state in states:
 target=labels==state;ds=[];eligible=[]
 for patient in sorted(set(patients)):
  t=target&(patients==patient);r=(~target)&(patients==patient)
  if t.sum()<20 or r.sum()<20:continue
  tc=np.asarray(counts[t].sum(0)).ravel();rc=np.asarray(counts[r].sum(0)).ravel();delta=np.log2(tc/tc.sum()*1e6+1)-np.log2(rc/rc.sum()*1e6+1)
  ds.append(delta);eligible.append(patient)
  pe.append(pd.DataFrame({'state':state,'patient':patient,'gene_id':a.var_names,'symbol':a.var.symbol.to_numpy(),'log2CPM_difference':delta}))
  for program in programs:
   pg.append({'state':state,'patient':patient,'program':program,'target_cells':int(t.sum()),'other_cells':int(r.sum()),'score_difference':float(a.obs.loc[t,program].median()-a.obs.loc[r,program].median())})
 arr=np.asarray(ds,dtype=np.float64);n=len(ds);effect=np.median(arr,axis=0) if n else np.full(a.n_vars,np.nan)
 m=pd.DataFrame({'state':state,'gene_id':a.var_names,'symbol':a.var.symbol.to_numpy(),'paired_patients':n,'median_log2CPM_difference':effect,'positive_patients':(arr>0).sum(0) if n else 0,'negative_patients':(arr<0).sum(0) if n else 0,'target_detection_pct':np.asarray((counts[target]>0).mean(0)).ravel()*100,'other_detection_pct':np.asarray((counts[~target]>0).mean(0)).ravel()*100})
 if n>=3:
  signs=np.array(list(product([-1,1],repeat=n)));pv=(np.abs(signs@arr)>=np.abs(arr.sum(0))[None,:]-1e-10).mean(0);assert np.all(pv>=2/(2**n)), 'Exact enumeration P below attainable minimum'
  m['P_value']=pv;m['FDR_BH']=multipletests(pv,method='fdr_bh')[1]
 else:m['P_value']=np.nan;m['FDR_BH']=np.nan
 m['FDR_family']='Harmony_state_'+state+'_vs_remaining_CAF_all_genes';markers.append(m)
 audit.append({'state':state,'cells':int(target.sum()),'eligible_patients':eligible,'paired_n':n,'theoretical_min_P':2/2**n if n>=3 else None,'nominal_P_lt05':int((m.P_value<.05).sum()),'FDR_lt05':int((m.FDR_BH<.05).sum())})
 for altkey in [c for c in a.obs if c.startswith('harmony_theta') and c!=key]:
  alt=a.obs[altkey].astype(str).to_numpy();options=[]
  for ac in set(alt):
   v=alt==ac;inter=int((v&target).sum());options.append((inter/int((v|target).sum()),ac,inter/int(target.sum()),inter/int(v.sum())))
  j,ac,rec,prec=max(options)
  stability.append({'state':state,'comparison':altkey,'best_match':ac,'Jaccard':j,'target_recall':rec,'match_precision':prec})
pd.concat(markers).to_csv(out/'patient_paired_gene_effects_all.tsv.gz',sep='\t',index=False)
pd.concat(pe).to_csv(out/'individual_patient_gene_effects.tsv.gz',sep='\t',index=False)
pg=pd.DataFrame(pg);pg.to_csv(out/'individual_patient_program_effects.tsv',sep='\t',index=False)
pr=[]
for (s,program),d in pg.groupby(['state','program']):
 vals=d.score_difference.to_numpy();n=len(vals);pv=float((np.abs(np.array(list(product([-1,1],repeat=n)))@vals)>=abs(vals.sum())-1e-10).mean()) if n>=3 else np.nan
 pr.append({'state':s,'program':program,'paired_patients':n,'median_score_difference':np.median(vals),'patients_positive':int((vals>0).sum()),'patients_negative':int((vals<0).sum()),'P_value':pv})
pr=pd.DataFrame(pr);pr['FDR_BH']=np.nan
for s,d in pr.groupby('state'):
 valid=d.P_value.notna()
 if valid.any():pr.loc[d.index[valid],'FDR_BH']=multipletests(d.loc[valid,'P_value'],method='fdr_bh')[1]
pr.to_csv(out/'patient_paired_program_effects.tsv',sep='\t',index=False)
pd.DataFrame(stability).to_csv(out/'parameter_state_overlap.tsv',sep='\t',index=False)
m=pd.concat(markers);top=m[(m.paired_patients>=3)&(m.target_detection_pct>=20)&(m.positive_patients>=m.paired_patients-1)].sort_values(['state','median_log2CPM_difference'],ascending=[True,False]).groupby('state').head(30)
top.to_csv(out/'cross_patient_top30_genes.tsv',sep='\t',index=False)
q=a.obs.quality_sensitivity_keep.to_numpy();pd.crosstab(a.obs.loc[q,'patient_id'],a.obs.loc[q,key]).to_csv(out/'quality_retained_patient_state_counts.tsv',sep='\t')
(out/'statistical_audit.json').write_text(json.dumps(audit,indent=2));(out/'COMPLETE').write_text('HARMONY_STATE_VALIDATION_COMPUTE_COMPLETE\n');print(json.dumps(audit));print('HARMONY_STATE_VALIDATION_COMPUTE_COMPLETE',flush=True)
