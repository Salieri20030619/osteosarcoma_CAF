from pathlib import Path
import json
import numpy as np,pandas as pd,scanpy as sc
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();r=p/'remediation_20260917';out=r/'program_stability';out.mkdir(exist_ok=True)
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');parts=[sc.read_h5ad(f) for f in sorted((r/'decontx_comparison').glob('*_strict_corrected.h5ad'))];b=sc.concat(parts,join='inner',merge='same')[a.obs_names,a.var_names].copy();b.var=a.var.copy();sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
labels=pd.read_csv(r/'state_final_sensitivity/corrected_recluster_cells.tsv.gz',sep='\t',dtype={'corrected_state':str}).set_index('cell_id').reindex(a.obs_names);old=a.obs['harmony_theta2_r0.5'].astype(str).to_numpy();new=labels.corrected_state.to_numpy();patients=a.obs.patient_id.astype(str).to_numpy();panels=json.loads((p/'analysis/gene_panels.json').read_text());names=['Contractile','Matrix_remodeling','Inflammatory','Antigen_presentation','Osteogenic','Mineral_calcium'];rng=np.random.default_rng(20260917)
whole=[];effects=[];cell=[]
for version,x in [('original',a),('corrected',b)]:
 values={}
 for name in names:
  ix=np.where(a.var.symbol.isin(panels[name]))[0];values[name]=np.asarray(x.X[:,ix].mean(1)).ravel()
 for pat in sorted(set(patients)):
  within=patients==pat
  for name,v in values.items():whole.append({'patient':pat,'version':version,'program':name,'cells':int(within.sum()),'median_mean_log_expression':float(np.median(v[within]))})
 for design,lab in [('frozen_states',old),('reclustered_states',new)]:
  for original_state,cluster in [('H0','0'),('H2','2' if design=='frozen_states' else '1')]:
   for pat in sorted(set(patients)):
    t=(patients==pat)&(lab==cluster);o=(patients==pat)&(lab!=cluster);nt=int(t.sum());no=int(o.sum())
    if min(nt,no)<5:continue
    for name,v in values.items():
     tv=v[t];ov=v[o];effect=float(tv.mean()-ov.mean())
     # Conditional resampling uncertainty; no patient replication or P values implied.
     boot=tv[rng.integers(nt,size=(500,nt))].mean(1)-ov[rng.integers(no,size=(500,no))].mean(1)
     row={'patient':pat,'version':version,'design':design,'state':original_state,'program':name,'target_cells':nt,'other_cells':no,'mean_log_expression_difference':effect,'resampling_2.5pct':float(np.quantile(boot,.025)),'resampling_97.5pct':float(np.quantile(boot,.975)),'resampling_fraction_same_direction':float((boot*effect>0).mean()),'eligible_primary_20_per_side':min(nt,no)>=20,'interpretation':'within-patient cell-resampling precision only; not biological inference'}
     for size in [10,20]:
      if min(nt,no)>=size:
       ds=np.array([rng.choice(tv,size,replace=False).mean()-rng.choice(ov,size,replace=False).mean() for _ in range(200)])
       row[f'downsample_{size}_fraction_same_direction']=float((ds*effect>0).mean())
      else:row[f'downsample_{size}_fraction_same_direction']=np.nan
     effects.append(row)
 for i,c in enumerate(a.obs_names):cell.append({'cell_id':c,'patient':patients[i],'version':version,**{k:float(v[i]) for k,v in values.items()}})
w=pd.DataFrame(whole);e=pd.DataFrame(effects);w.to_csv(out/'whole_patient_programs.tsv',sep='\t',index=False);e.to_csv(out/'within_patient_program_stability.tsv',sep='\t',index=False);pd.DataFrame(cell).to_csv(out/'cell_program_expression.tsv.gz',sep='\t',index=False)
summary=[]
for (ver,design,state,name),z in e.groupby(['version','design','state','program']):
 for scope,zz in [('primary20',z[z.eligible_primary_20_per_side]),('descriptive_min5',z)]:
  if not len(zz):continue
  v=zz.mean_log_expression_difference.to_numpy();summary.append({'version':ver,'design':design,'state':state,'program':name,'scope':scope,'patients':len(v),'median_patient_effect':float(np.median(v)),'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'leave_one_patient_out_median_min':float(min(np.median(np.delete(v,i)) for i in range(len(v)))) if len(v)>1 else np.nan,'leave_one_patient_out_median_max':float(max(np.median(np.delete(v,i)) for i in range(len(v)))) if len(v)>1 else np.nan})
pd.DataFrame(summary).to_csv(out/'patient_direction_summary.tsv',sep='\t',index=False)
(out/'parameters.json').write_text(json.dumps({'programs':{k:panels[k] for k in names},'value':'mean log1p normalized expression over covered program genes; compare same program only','bootstrap_replicates':500,'balanced_subsampling_replicates':200,'seed':20260917,'identity_and_primary_thresholds_changed':False,'min5':'descriptive measurement precision only; not alternative primary test','pvalues':'none'},indent=2));(out/'COMPLETE').write_text('PROGRAM_STABILITY_COMPLETE\n');print('PROGRAM_STABILITY_COMPLETE',len(e),flush=True)
