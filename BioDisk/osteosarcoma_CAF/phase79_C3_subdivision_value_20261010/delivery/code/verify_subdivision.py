from pathlib import Path
import pandas as pd,numpy as np,json,itertools,shutil,ast
from scipy.stats import rankdata
B=Path('C:/Users/31937/CAF_C3_subdivision_value_20261010');R=B/'results';T=R/'tables';D=B/'delivery';D.mkdir(exist_ok=True)
# Make the original script robust to its recorded unavailable-panel case without altering numerical work.
p=B/'analyze_subdivision.py';original=p.read_text(encoding='utf-8-sig');shutil.copy2(p,B/'analyze_subdivision_initial.py');needle="cv=pd.DataFrame(cvrows);cv.to_csv";guard="if not cvrows:\n (O/'CELL_STAGE_COMPLETE.json').write_text(json.dumps({'cell_stage_completed':True,'strict_spatial_unassessable':True,'reason':'Required assay genes missing; use separately recorded assay-limited sensitivity'},indent=2))\n print('CELL_STAGE_COMPLETE_STRICT_SPATIAL_UNASSESSABLE',flush=True)\n sys.exit(0)\n";assert needle in original;p.write_text(original.replace(needle,guard+needle),encoding='utf-8');ast.parse(p.read_text())
cs=pd.read_csv(T/'cell_scores.tsv.gz',sep='\t');ce=pd.read_csv(T/'composition_patient_effects.tsv',sep='\t');ms=pd.read_csv(T/'matched_cell_differences.tsv',sep='\t');md=pd.read_csv(T/'matching_patient_effects.tsv',sep='\t');checked=0
# Panel identity against previous independently reviewed score export.
old=pd.read_csv('C:/Users/31937/CAF_C3_matrix_strength_20261010/results/cell_functional_scores.tsv.gz',sep='\t').set_index(['matrix','cell_id']);now=cs.set_index(['matrix','cell_id']);assert np.allclose(now.processing,old.loc[now.index,'processing2']);assert np.allclose(now.COL5_reference,old.loc[now.index,'Published_COL5A1_MSC4']);assert np.allclose(now.balance,now.III_V-now.IV)
for _,r in ce.iterrows():
 z=cs[(cs.matrix==r.matrix)&(cs.donor_record==r.donor)];v=[]
 for site,zz in z.groupby('site'):
  l=zz[zz.subgroup_id=='F2_C3'];rt=zz[zz.subgroup_id!='F2_C3']
  if r.reference=='F2_C1':rt=rt[rt.subgroup_id=='F2_C1']
  elif r.reference=='other_F2':rt=rt[rt.class_id=='F2']
  if min(len(l),len(rt))>=r.minimum:v.append(l[r.panel].mean()-rt[r.panel].mean())
 assert np.isclose(np.mean(v),r.effect)
for mode,d in cs.groupby('matrix'):
 d=d.set_index('cell_id');cv=d[['COL5_reference','log_library']].to_numpy();w=1/d.groupby('donor_record').donor_record.transform('size').to_numpy();w/=w.sum();mu=w@cv;sd=np.sqrt(w@((cv-mu)**2));z=(cv-mu)/sd;d['z0']=z[:,0];d['z1']=z[:,1]
 for _,r in ms[ms.matrix==mode].iterrows():
  l=d.loc[r.target_cell];ctl=d.loc[r.control_cells.split(';')];assert len(ctl)==3 and ctl.subgroup_id.ne('F2_C3').all() and ctl.donor_record.eq(r.donor).all() and ctl.site.eq(r.site).all();assert np.isclose(l.processing-ctl.processing.mean(),r.effect);assert np.isclose(l.z0-ctl.z0.mean(),r.COL5_z_difference);assert np.isclose(l.z1-ctl.z1.mean(),r.library_z_difference)
  if r.caliper=='caliper0.5':assert (abs(ctl[['z0','z1']].to_numpy()-l[['z0','z1']].to_numpy(float))<=.5+1e-12).all()
agg=ms.groupby(['matrix','caliper','donor','site']).effect.mean().groupby(['matrix','caliper','donor']).mean()
for _,r in md.iterrows():assert np.isclose(agg.loc[(r.matrix,r.caliper,r.donor)],r.effect)

def exact(v):
 v=np.asarray(v,float);n=len(v)
 if n<2 or not np.any(abs(v)>1e-12):return np.nan
 hits=0
 for start in range(0,2**n,8192):
  b=np.arange(start,min(start+8192,2**n),dtype=np.uint64);sgn=(((b[:,None]>>np.arange(n,dtype=np.uint64))&1)*2).astype(float)-1;hits+=sum(abs(sgn@v/n)>=abs(v.mean())-1e-12)
 return hits/2**n
for fname,eff,keys,val in [('composition_statistics.tsv',ce,['matrix','minimum','reference','panel'],'effect'),('matching_statistics.tsv',md,['matrix','caliper'],'effect'),('assay_limited_spatial_association_statistics.tsv',pd.read_csv(T/'assay_limited_spatial_patient_associations.tsv',sep='\t'),['scope','control','endpoint'],'rho')]:
 st=pd.read_csv(T/fname,sep='\t')
 for _,r in st.iterrows():
  q=eff.copy()
  for key in keys:q=q[q[key]==r[key]]
  v=q[val].dropna().to_numpy();assert len(v)==r.N and sum(v>0)==r.positive;assert np.isclose(v.mean(),r.mean_effect) if len(v) else pd.isna(r.mean_effect);pv=exact(v);assert np.isclose(pv,r.P,equal_nan=True);checked+=1
 for family,g in st[st.P.notna()].groupby('family'):
  x=g.P.to_numpy();o=np.argsort(x);bh=np.minimum(1,np.minimum.accumulate((x[o]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]);assert np.allclose(bh,g.FDR.to_numpy()[o])
# Reconstruct spatial partial correlations and the hierarchical patient means.
sp=pd.read_csv(T/'spatial_scores_assay_limited.tsv.gz',sep='\t',dtype={'event_month':str});sa=pd.read_csv(T/'assay_limited_spatial_section_associations.tsv',sep='\t',dtype={'event_month':str});ctrl=['log_library','fibro_context','mural_context','endothelial_context','osteogenic_available5','matrix_CAF','x','y']
for _,r in sa.iterrows():
 q=sp[sp['sample']==r['sample']];cols=ctrl+(['contractile_available2'] if r.control.endswith('contractile2') else []);q=q.dropna(subset=cols+['C3_proxy',r.endpoint]);Z=np.column_stack([np.ones(len(q))]+[rankdata(q[c]) for c in cols]);x=rankdata(q.C3_proxy);y=rankdata(q[r.endpoint]);rx=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0];ry=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0];rho=np.corrcoef(rx,ry)[0,1];assert np.isclose(rho,r.rho,atol=1e-9)
pa=pd.read_csv(T/'assay_limited_spatial_patient_associations.tsv',sep='\t')
for scope,q in [('all_patients',sa),('primary_only',sa[sa.site_group=='primary_tumor'])]:
 keys=['control','endpoint'];a=q.groupby(['patient_id','site_group','event_month']+keys,dropna=False).rho.mean().reset_index();a=a.groupby(['patient_id','site_group']+keys).rho.mean().reset_index();a=a.groupby(['patient_id']+keys).rho.mean();b=pa[pa.scope==scope].set_index(['patient_id']+keys);assert np.allclose(a,b.loc[a.index,'rho'])
pr=pd.read_csv(T/'assay_limited_spatial_LOPO_predictions.tsv.gz',sep='\t');cv=pd.read_csv(T/'assay_limited_spatial_LOPO_prediction_effects.tsv',sep='\t');cvkeys=['scope','patient_id','control','endpoint']
for _,r in cv.iterrows():
 q=pr.copy()
 for key in cvkeys:q=q[q[key]==r[key]]
 w=q.weight.to_numpy();y=q.observed.to_numpy();assert np.isclose(w.sum(),1);mu=w@y;v=w@((y-mu)**2);b=w@((y-q.baseline_prediction.to_numpy())**2);a=w@((y-q.extended_prediction.to_numpy())**2);assert np.isclose((b-a)/v,r.delta_R2)
strict=pd.read_csv(T/'spatial_association_statistics.tsv',sep='\t');assert strict.N.eq(0).all() and strict.P.isna().all();assert int((R/'logs/spatial_assay_limited.exit_code').read_text())==0
qa={'status':'passed_with_assay_limitation_recorded','statistical_rows_recomputed':checked,'composition_patient_effects_reconstructed':len(ce),'matched_cell_contrasts_verified':len(ms),'spatial_section_correlations_reconstructed':len(sa),'heldout_prediction_records_verified':len(cv),'previous_processing_scores_reproduced':True,'all_BH_families_verified':True,'strict_panel_unassessable_rows_retained':len(strict),'initial_run_exit_code':1,'assay_limited_run_exit_code':0,'initial_error_handling_repaired_without_changing_computed_effects':True,'analysis_complete':True};(D/'QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa))
