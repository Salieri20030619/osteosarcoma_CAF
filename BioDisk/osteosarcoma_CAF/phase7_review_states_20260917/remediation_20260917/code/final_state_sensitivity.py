from pathlib import Path
import json
import numpy as np,pandas as pd,scanpy as sc,harmonypy as hm
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();r=p/'remediation_20260917';out=r/'state_final_sensitivity';out.mkdir(exist_ok=True)
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');key='harmony_theta2_r0.5'
parts=[sc.read_h5ad(f) for f in sorted((r/'decontx_comparison').glob('*_strict_corrected.h5ad'))];b=sc.concat(parts,join='inner',merge='same');b=b[a.obs_names,a.var_names].copy();assert b.n_obs==2047;b.var=a.var.copy();b.layers['corrected_counts']=b.X.copy();sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
markers=['ACTA2','COL3A1','FN1','RUNX2','OMD','IL6','CXCL12','CHI3L1','PLA2G2A','GPC3','SLPI','CD74','HLA-DRA','CIITA','CD80','CD86'];rows=[]
patients=a.obs.patient_id.astype(str).to_numpy();labels=a.obs[key].astype(str).to_numpy()
for state in ['0','2']:
 for pat in sorted(set(patients)):
  t=(patients==pat)&(labels==state);rest=(patients==pat)&(labels!=state)
  if t.sum()<20 or rest.sum()<20:continue
  for version,counts in [('original',a.layers['counts']),('corrected',b.layers['corrected_counts'])]:
   tc=np.asarray(counts[t].sum(0),dtype=np.float64).ravel();rc=np.asarray(counts[rest].sum(0),dtype=np.float64).ravel();delta=np.log2(tc/tc.sum()*1e6+1)-np.log2(rc/rc.sum()*1e6+1)
   for gene in markers:
    js=np.where(a.var.symbol==gene)[0]
    if len(js)!=1:continue
    rows.append({'state':state,'patient':pat,'version':version,'gene':gene,'log2CPM_difference':float(delta[js[0]])})
d=pd.DataFrame(rows);d.to_csv(out/'individual_patient_gene_effects.tsv',sep='\t',index=False)
s=d.groupby(['state','version','gene']).agg(patients=('patient','nunique'),median_effect=('log2CPM_difference','median'),positive=('log2CPM_difference',lambda x:int((x>0).sum())),negative=('log2CPM_difference',lambda x:int((x<0).sum()))).reset_index();s.to_csv(out/'gene_effect_summary.tsv',sep='\t',index=False)
# Recluster corrected expression with frozen original HVG set and Harmony parameters.
hvg=b[:,a.var.highly_variable.to_numpy()].copy();sc.pp.scale(hvg,max_value=10);sc.tl.pca(hvg,n_comps=30,random_state=20260917,svd_solver='arpack');h=hm.run_harmony(hvg.obsm['X_pca'],hvg.obs,'patient_id',theta=2,random_state=20260917,ncores=4,max_iter_harmony=20)
x=np.asarray(h.Z_corr);hvg.obsm['X_harmony']=x if x.shape[0]==hvg.n_obs else x.T;sc.pp.neighbors(hvg,use_rep='X_harmony',n_neighbors=15,random_state=20260917);sc.tl.leiden(hvg,key_added='corrected_state',resolution=.5,flavor='igraph',n_iterations=2,directed=False,random_state=20260917)
hvg.obs[['patient_id',key,'corrected_state']].to_csv(out/'corrected_recluster_cells.tsv.gz',sep='\t',index_label='cell_id');new=hvg.obs.corrected_state.astype(str).to_numpy();ct=pd.crosstab(hvg.obs.patient_id,new);ct.to_csv(out/'corrected_recluster_patient_counts.tsv',sep='\t')
result=[]
for state in ['0','2']:
 t=labels==state;matches=[]
 for c in set(new):
  q=new==c;inter=int((q&t).sum());matches.append((inter/int((q|t).sum()),c,inter/int(t.sum()),inter/int(q.sum())))
 j,c,rec,prec=max(matches);result.append({'original_state':state,'best_corrected_match':c,'Jaccard':j,'recall':rec,'precision':prec,'paired_patients':int(((ct[c]>=20)&((ct.sum(1)-ct[c])>=20)).sum())})
pd.DataFrame(result).to_csv(out/'corrected_state_overlap.tsv',sep='\t',index=False)
changed=[]
for f in (r/'cnv_comparison').glob('*_cell_comparison.tsv.gz'):
 q=pd.read_csv(f,sep='\t');changed.extend(q.loc[q.strict_CAF&(q.call_immune!='low_signal'),'cell_id'].tolist())
assert len(changed)==25;keep=~a.obs_names.isin(changed);obs=a.obs.loc[keep];tab=pd.crosstab(obs.patient_id,obs[key]);tab.to_csv(out/'exclude_reference_sensitive_25_patient_counts.tsv',sep='\t');removed=a.obs.loc[a.obs_names.isin(changed),['patient_id',key]];removed.to_csv(out/'reference_sensitive_25_cells.tsv',sep='\t',index_label='cell_id')
cnv=[]
for state in ['0','2']:cnv.append({'state':state,'retained':int(tab[state].sum()),'paired_patients':int(((tab[state]>=20)&((tab.sum(1)-tab[state])>=20)).sum())})
(out/'summary.json').write_text(json.dumps({'corrected_expression_recluster':result,'CNV_25_exclusion_sensitivity':cnv,'original_cell_identity_updated':False,'corrected_effect_tests':'descriptive only; no additional inferential tests','correction_model':'DecontX broad-lineage z; not independent ground truth'},indent=2));(out/'COMPLETE').write_text('FINAL_STATE_SENSITIVITY_COMPLETE\n');print(json.dumps(result));print(json.dumps(cnv))
