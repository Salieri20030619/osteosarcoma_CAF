from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json,hashlib,warnings
warnings.simplefilter('ignore',pd.errors.PerformanceWarning)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919';assert (O/'UNIFIED_GATES_COMPLETE').exists()
e=pd.read_csv(O/'tables/all_cell_unified_evidence.tsv.gz',sep='\t',index_col='cell_id',low_memory=False);flags=pd.Series(False,index=e.index)
for f in [P/'phase3_identity_20260914/identity/all_cell_identity_evidence.tsv.gz',P/'phase7_review_states_20260917/identity/all_cell_identity_evidence.tsv.gz',P/'geo_identity_20260919/tables/all_cell_identity_evidence.tsv.gz']:
 cols=pd.read_csv(f,sep='\t',nrows=0).columns;keys=[c for c in ['doublet_review_flag','doublet_review_exclude','doublet_exclude'] if c in cols];d=pd.read_csv(f,sep='\t',usecols=['cell_id']+keys).set_index('cell_id');v=np.zeros(len(d),bool)
 for k in keys:v|=d[k].fillna(False).astype(str).str.lower().isin(['true','1','1.0']).to_numpy()
 flags.loc[d.index]=v
sets=['unified_strict','osteo_sensitive','competitor_sensitive','mural_reference','trace_count_sensitive'];before={k:int(e[k].sum()) for k in sets};e['independent_doublet_excluded']=flags
for k in sets:e[k]=e[k]&~flags
# Reapplying the guard is idempotent.
e.loc[flags,'exclusion_reasons']=e.loc[flags,'exclusion_reasons'].fillna('').str.replace('independent_doublet_exclusion;','',regex=False)+'independent_doublet_exclusion;'
assert not e.loc[e.unified_strict,'independent_doublet_excluded'].any()
objects=pd.read_csv(O/'resources/object_manifest.tsv',sep='\t');records=[]
for r in objects.to_dict('records'):
 b=ad.read_h5ad(r['path']);bo=e.loc[b.obs_names].copy();keep=bo[sets].any(axis=1).to_numpy();b=b[keep].copy();bo=bo.iloc[np.flatnonzero(keep)].copy();b.obs=bo
 b.obs.index=pd.Index(b.obs.index.to_numpy(dtype=object),dtype=object);b.var.index=pd.Index(b.var.index.to_numpy(dtype=object),dtype=object)
 for frame in [b.obs,b.var]:
  for col in frame:
   if isinstance(frame[col].dtype,pd.StringDtype):frame[col]=frame[col].astype(object)
 b.write_h5ad(r['path'],compression='gzip');r['cells']=b.n_obs;records.append(r)
pd.DataFrame(records).to_csv(O/'resources/object_manifest.tsv',sep='\t',index=False)
summary=[];reasons=pd.read_csv(O/'tables/old_joint_exclusion_reason_counts.tsv',sep='\t');reasonrows=[]
for sm,g in e.groupby('sample'):
 old=g.old_joint_member;summary.append({'sample':sm,'dataset':g.dataset.iloc[0],'source_stage':g.source_stage.iloc[0],'technology':g.technology.iloc[0],'retained_cells':len(g),'ambient_checked':bool(g.ambient_checked.iloc[0]),'old_joint_cells':int(old.sum()),**{k:int(g[k].sum()) for k in sets},'old_retained_in_new':int((old&g.unified_strict).sum()),'old_removed':int((old&~g.unified_strict).sum()),'new_vs_old':int((~old&g.unified_strict).sum()),'independent_doublet_excluded':int(g.independent_doublet_excluded.sum())});reasonrows.append({'sample':sm,'dataset':g.dataset.iloc[0],'reason':'independent_doublet_exclusion','old_joint_excluded_cells':int((old&g.independent_doublet_excluded).sum())});g.to_csv(O/'tables'/(sm+'_unified_evidence.tsv.gz'),sep='\t',index_label='cell_id')
e.to_csv(O/'tables/all_cell_unified_evidence.tsv.gz',sep='\t',index_label='cell_id');pd.DataFrame(summary).to_csv(O/'tables/sample_gate_changes.tsv',sep='\t',index=False);pd.concat([reasons[reasons.reason!='independent_doublet_exclusion'],pd.DataFrame(reasonrows)]).to_csv(O/'tables/old_joint_exclusion_reason_counts.tsv',sep='\t',index=False)
# Recompute molecular negative-control gates, including the independent quality guard.
raw=e[[c for c in e if c.startswith('raw_')]].rename(columns=lambda c:c[4:]);cor=e[[c for c in e if c.startswith('corrected_')]].rename(columns=lambda c:c[10:])
def conflict(d,frac=.2,osteo=.2,mural=True,det='det_'):
 out=np.zeros(len(d),bool)
 for name,n,pct in [('Mural',3,frac),('OC',2,frac/2),('T_NK',3,frac),('Myeloid',4,frac),('Endo',3,frac),('Epi',3,frac),('Myogenic',2,frac)]:
  if name=='Mural' and not mural:continue
  pp='OC_specific' if name=='OC' else name;out|=(d[det+name]>=n)&(d['pct_'+pp]>pct)
 out|=(d[det+'Osteo']>=3)&(d.pct_Osteo>osteo);immune=sum(d['pct_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC']);mes=d.pct_Fib+d.pct_Osteo+d.pct_Mural
 mix=((immune>frac)&(mes>frac))|((d.pct_Endo>frac)&(mes>frac))|((d.pct_Epi>frac)&(mes>frac));return np.asarray(out|mix)
fixed=e.dual_low&e.reference_QC_pass&~e.CTA_context&~flags&(raw.det_CTA<2)&(raw.det_Fib>=2)&(e.n_genes>=500)&(e.cluster_both_high_fraction.fillna(1)<.2)&(e.cluster_any_high_fraction.fillna(1)<.2)&e.ambient_checked&(e.corrected_library_counts>0)
cal=[]
for key,frac,osteo,mural,det in [('unified_strict',.2,.2,True,'det_'),('osteo_sensitive',.2,.5,True,'det_'),('competitor_sensitive',.5,.5,True,'det_'),('mural_reference',.2,.2,False,'det_'),('trace_count_sensitive',.2,.2,True,'pos_')]:
 mask=fixed&~conflict(raw,frac,osteo,mural)&~conflict(cor,frac,osteo,mural,det)&(cor[det+'Fib']>=2)&(cor[det+'CTA']<2)
 if key!='mural_reference':assert np.array_equal((mask&e.core_scope).to_numpy(),e[key].to_numpy())
 for (sm,role),g in e[e.cnv_role.isin(['holdout_immune','holdout_endo'])].groupby(['sample','cnv_role']):cal.append({'sample':sm,'dataset':g.dataset.iloc[0],'rule':key,'role':role,'cells':len(g),'molecular_pass_without_scope_gate':int(mask.loc[g.index].sum()),'note':'descriptive negative-control leakage; not independent stromal truth'})
pd.DataFrame(cal).to_csv(O/'tables/negative_control_gate_checks.tsv',sep='\t',index=False)
protocol=json.loads((O/'PROTOCOL.json').read_text());protocol['mandatory_quality_guard']='Any existing reviewed doublet-exclusion flag remains exclusionary for all candidate sets; detector provenance preserved. Stage code/apply_quality_guard.py is mandatory before batch models.';(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
status=json.loads((O/'STATUS.json').read_text());status.update({k:int(e[k].sum()) for k in sets});status['old_retained_in_new']=int((e.unified_strict&e.old_joint_member).sum());status['quality_guard_passed']=True;status['counts_before_quality_guard']=before;(O/'STATUS.json').write_text(json.dumps(status,indent=2));(O/'QUALITY_GUARD_COMPLETE').write_text('QUALITY_GUARD_COMPLETE\n');print(json.dumps(status),flush=True)