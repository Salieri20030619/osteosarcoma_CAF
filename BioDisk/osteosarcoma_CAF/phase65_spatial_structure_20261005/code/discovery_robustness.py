from pathlib import Path
import os,sys,json,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005';A=B/'phase63_matrix_resistance_20261005'
sys.path.insert(0,str(B/'phase46_spatial_functional_20261001/code'));from common import test_table
cfg=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text());sets=cfg['structure'];defs={'full':sets,'drop_dominant_pair':{k:[g for g in gs if g not in cfg['drop_dominant'][k]] for k,gs in sets.items()}}
ob=pd.read_csv(A/'tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
roster=pd.read_csv(A/'tables/paired_patient_site_roster.tsv.gz',sep='\t');z=np.load(A/'objects/mechanism_counts.npz')
print('INPUT',len(ob),'npz keys',z.files,flush=True)
if 'cell_ids' in z:assert list(z['cell_ids'])==list(ob.index)
if 'cells' in z:assert list(z['cells'])==list(ob.index)
counts={k:z[k] for k in ['raw','corrected']};den={k:z[k+'_library'] for k in counts};assay=z['assayed'];gi={g:i for i,g in enumerate(z['genes'])}
mem={(v,p):np.array([gi[g] for g in gs]) for v,ds in defs.items() for p,gs in ds.items()}
rows=[]
for (ds,st,mode,scope,minimum),rg in roster.groupby(['dataset','subgroup_id','matrix','scope','minimum_cells']):
 for donor,dg in rg.groupby('donor_record'):
  for r in dg.itertuples():
   base=(ob.dataset.eq(ds)&ob.donor_record.eq(donor)&ob.site.eq(r.site)).to_numpy()
   if scope=='Fibro':base &= ob.lineage_stratum.eq('Fibro').to_numpy()
   il=np.flatnonzero(base&ob.subgroup_id.eq(st).to_numpy());ir=np.flatnonzero(base&ob.subgroup_id.ne(st).to_numpy())
   assert len(il)==r.target_cells and len(ir)==r.reference_cells
   x=counts[mode];lib=den[mode];present=assay[np.r_[il,ir]].all(0)
   dv=np.log2(x[il].sum(0)/lib[il].sum()*1e4+.1)-np.log2(x[ir].sum(0)/lib[ir].sum()*1e4+.1)
   for (v,p),jj in mem.items():
    have=jj[present[jj]];ok=len(have)>=3 and len(have)/len(jj)>=.8
    rows.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'donor_record':donor,'site':r.site,'version':v,'program':p,'effect':(0.0 if abs(float(dv[have].mean()))<1e-12 else float(dv[have].mean())) if ok else np.nan,'target_cells':len(il),'reference_cells':len(ir)})
f=pd.DataFrame(rows);f.to_csv(O/'tables/discovery_site_effects.tsv.gz',sep='\t',index=False)
keys=['dataset','subgroup_id','matrix','scope','minimum_cells','version','program']
pat=f.groupby(keys+['donor_record'],dropna=False).effect.mean().reset_index();pat.to_csv(O/'tables/discovery_patient_effects.tsv',sep='\t',index=False)
stat=test_table(pat,keys,'effect',lambda r:f'{r.dataset}|{r.matrix}|{r.scope}|min{r.minimum_cells}|{r.version}|all33_structure3')
stat.to_csv(O/'tables/discovery_structure_statistics.tsv',sep='\t',index=False)
rob=[]
for k,d in pat.groupby(keys,dropna=False):
 v=d.effect.dropna().to_numpy();n=len(v);loo=(v.sum()-v)/(n-1) if n>1 else np.array([np.nan])
 rob.append(dict(zip(keys,k),N=n,mean_effect=v.mean() if n else np.nan,positive_patients=int(sum(v>0)),negative_patients=int(sum(v<0)),loo_min=float(min(loo)) if n>1 else np.nan,loo_max=float(max(loo)) if n>1 else np.nan,all_loo_same_direction=bool(np.all(np.sign(loo)==np.sign(v.mean()))) if n>1 else False,evidence_level='exploratory_N_lt10' if 2<=n<10 else 'descriptive_only' if n<2 else 'patient_inference'))
pd.DataFrame(rob).to_csv(O/'tables/discovery_patient_leave_one_out.tsv',sep='\t',index=False)
cols=['subgroup_id','working_name','expression_confidence','functional_interpretation_confidence','lineage_confidence','F4_major_branch','cell_id']
# Keep cell confidence and assignment unchanged for traceability.
ob.reset_index()[cols].to_csv(O/'tables/original_cell_assignments_and_confidence.tsv.gz',sep='\t',index=False)
(O/'DISCOVERY_COMPLETE.json').write_text(json.dumps({'success':True,'cells':len(ob),'subgroups':int(ob.subgroup_id.nunique()),'patient_effect_rows':len(pat),'statistics_rows':len(stat),'at':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print('DISCOVERY_COMPLETE',len(stat),flush=True)
