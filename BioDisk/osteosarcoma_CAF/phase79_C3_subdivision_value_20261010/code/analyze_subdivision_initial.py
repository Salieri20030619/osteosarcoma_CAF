from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,sys,itertools,time
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase79_C3_subdivision_value_20261010'
for z in ['code','tables','resources','logs']:(O/z).mkdir(parents=True,exist_ok=True)
III=['COL3A1','COL5A1','COL5A2','COL5A3'];IV=['COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6'];panels={'III_V':III,'IV':IV,'COL5_reference':['COL5A1','COL16A1','EMILIN1','ADAMTS2'],'processing':['PCOLCE','LOXL2'],'C3_proxy':['CD248','POSTN','VCAN'],'matrix_CAF':['FAP','MMP11','COL1A1','COL11A1','POSTN','FN1'],'C1_proxy':['MYH11','ACTA2','TAGLN']};cfg=json.loads((P/'phase65_spatial_structure_20261005/resources/PROTOCOL_FROZEN.json').read_text());contexts=cfg['contexts'];genes=sorted(set(sum(panels.values(),[])+sum(contexts.values(),[])))
protocol={'frozen_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'questions':['C3 versus C1 collagen composition','processing increment beyond comparable known COL5 programme','spatial RNA information beyond broad CAF'], 'cell_population':'SCPCP000017 current Fibro excluding model-positive doublets; raw/corrected separately','scores':'Mean log1p(full library CP10k); III_V minus IV is relative transcript program balance, not collagen mass ratio','panels':panels,'cell_tests':'Within donor/site C3 minus C1, other F2 Fibro and other Fibro; site then donor equal; minimum1 and5; exact signflip at patient level N>=2 nonzero; three composition endpoints x three references one BH family per matrix/threshold','matching':'Within donor/site, other Fibro controls; equal-donor SD of COL5_reference and log1p(library); nearest3 controls per C3 cell, with replacement; primary requires >=3 controls within0.5 SD on EACH feature; no-caliper sensitivity retains all targets with >=3 controls; no relabeling; aggregate matched-cell deltas equally within site then donor. Patient N not pair N. One BH family across matrix/caliper processing comparisons. Matching is descriptive selection adjustment, not causal.','spatial':'GSE293065 existing26 sections19 patients, all and primary-only9 separately. Zero-library spots excluded, no other new filters. Require complete panel genes, absence not zero. Prior spot scores reproduced.','spatial_association':'Within-section partial Spearman C3_proxy versus III_V, IV, balance; rank-adjust depth, x,y, four context panels, known matrix_CAF; with and without C1_proxy. Repeated sections -> events -> sites -> patient, equal at each level.12 related hypotheses constitute one BH family; exact patient signflip.','spatial_prediction':'Leave-one-patient-out Ridge lambda1. Training standardization only, hierarchical equal patient/site/event/section weighting; base depth+four contexts+matrix_CAF; extended adds C3_proxy; sensitivity includes C1_proxy in both. Heldout delta_R2=(MSEbase-MSEextended)/heldout weighted outcome variance. All and primary-only trained separately. No P on overlapping cross-validation folds; report effects/consistency.','spatial_limit':'Proxy RNA score, not localization of confirmed C3 cells or physical collagen geometry; no outcome genes shared with predictors; cross-feature RNA correlation is not function','no_new_labels_or_master_edits':True}
(O/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2))

def exact(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v)
 if n<2 or not np.any(abs(v)>1e-12):return np.nan
 hits=0
 for start in range(0,2**n,16384):
  b=np.arange(start,min(2**n,start+16384),dtype=np.uint64);sgn=(((b[:,None]>>np.arange(n,dtype=np.uint64))&1)*2).astype(float)-1;hits+=np.sum(abs(sgn@v/n)>=abs(v.mean())-1e-12)
 return hits/2**n

def summarize(d,keys,val,family):
 rows=[]
 for k,z in d.groupby(keys,dropna=False):
  k=k if isinstance(k,tuple) else (k,);v=z[val].dropna().to_numpy();r=dict(zip(keys,k));r.update(N=len(v),mean_effect=np.mean(v) if len(v) else np.nan,median_effect=np.median(v) if len(v) else np.nan,positive=int(sum(v>0)),P=exact(v));r['family']=family(r);rows.append(r)
 out=pd.DataFrame(rows);out['FDR']=np.nan
 for fam,z in out[out.P.notna()].groupby('family'):out.loc[z.index,'FDR']=multipletests(z.P,method='fdr_bh')[1]
 return out
sys.path.insert(0,str(P/'phase66_doublet_repair_DSC1_20261005/code'));from current_metadata import load_current_cells
m=load_current_cells().set_index('cell_id');mp=pd.read_csv(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');gm=mp[mp.dataset=='SCPCP000017'].drop_duplicates('symbol').set_index('symbol').feature;ce=[];matches=[];matchroster=[];matchpat=[];cellout=[];roster=[]
for mode,path in [('raw',P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad'),('corrected',P/'phase58_F1F5_round1_20261004/objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==m.index.tolist();ii=np.flatnonzero(m.dataset.eq('SCPCP000017')&m.lineage_stratum.eq('Fibro')&~m.doublet_positive_any);d=m.iloc[ii][['dataset','sample','donor_record','site','subgroup_id','class_id']].copy();lib=np.asarray(a.X[ii].sum(1)).ravel();jj=a.var_names.get_indexer(gm.loc[genes]);assert (jj>=0).all();val=np.log1p(a.X[ii][:,jj].toarray()/lib[:,None]*1e4);d['library']=lib;d['log_library']=np.log1p(lib)
 for name,gs in panels.items():d[name]=val[:,[genes.index(g) for g in gs]].mean(1)
 d['balance']=d.III_V-d.IV;cellout.append(d.assign(matrix=mode).reset_index())
 for minimum in [1,5]:
  for refname in ['F2_C1','other_F2','other_Fibro']:
   for donor,dg in d.groupby('donor_record'):
    vv=[]
    for site,z in dg.groupby('site'):
     l=z[z.subgroup_id=='F2_C3'];r=z[z.subgroup_id!='F2_C3']
     if refname=='F2_C1':r=r[r.subgroup_id=='F2_C1']
     elif refname=='other_F2':r=r[r.class_id=='F2']
     if len(l):roster.append({'matrix':mode,'minimum':minimum,'reference':refname,'donor':donor,'site':site,'target':len(l),'control':len(r),'evaluable':min(len(l),len(r))>=minimum})
     if min(len(l),len(r))>=minimum:vv.append(l[['III_V','IV','balance']].mean().to_numpy()-r[['III_V','IV','balance']].mean().to_numpy())
    if vv:
     for panel,effect in zip(['III_V','IV','balance'],np.mean(vv,axis=0)):ce.append({'matrix':mode,'minimum':minimum,'reference':refname,'donor':donor,'panel':panel,'effect':effect})
 # Standardization of matching variables is source-local and equal-donor.
 cv=d[['COL5_reference','log_library']].to_numpy();w=1/d.groupby('donor_record').donor_record.transform('size').to_numpy();w/=w.sum();mu=w@cv;sd=np.sqrt(w@((cv-mu)**2));assert (sd>0).all();d['match_COL5_z']=(cv[:,0]-mu[0])/sd[0];d['match_library_z']=(cv[:,1]-mu[1])/sd[1]
 for cal in ['caliper0.5','no_caliper']:
  for (donor,site),z in d.groupby(['donor_record','site']):
   l=z[z.subgroup_id=='F2_C3'];r=z[z.subgroup_id!='F2_C3'];nv=0;effects=[]
   for cid,t in l.sort_index().iterrows():
    dist=(r[['match_COL5_z','match_library_z']].to_numpy()-t[['match_COL5_z','match_library_z']].to_numpy(float));eligible=np.all(abs(dist)<=.5,axis=1) if cal=='caliper0.5' else np.ones(len(r),bool);ix=np.flatnonzero(eligible)
    if len(ix)<3:continue
    sel=ix[np.argsort(np.sum(dist[ix]**2,axis=1),kind='stable')[:3]];matched=r.iloc[sel];nv+=1;effects.append(t.processing-matched.processing.mean());row={'matrix':mode,'caliper':cal,'donor':donor,'site':site,'target_cell':cid,'control_cells':';'.join(matched.index),'effect':t.processing-matched.processing.mean(),'COL5_z_difference':t.match_COL5_z-matched.match_COL5_z.mean(),'library_z_difference':t.match_library_z-matched.match_library_z.mean(),'target_processing':t.processing,'control_processing':matched.processing.mean()};matches.append(row)
   if len(l):matchroster.append({'matrix':mode,'caliper':cal,'donor':donor,'site':site,'C3_cells':len(l),'control_cells':len(r),'matched_C3_cells':nv})
   if effects:matchpat.append({'matrix':mode,'caliper':cal,'donor':donor,'site':site,'effect':np.mean(effects)})
 print('CELL_COMPLETE',mode,flush=True)
pd.concat(cellout).to_csv(O/'tables/cell_scores.tsv.gz',sep='\t',index=False);ce=pd.DataFrame(ce);ce.to_csv(O/'tables/composition_patient_effects.tsv',sep='\t',index=False);st=summarize(ce,['matrix','minimum','reference','panel'],'effect',lambda r:f"{r['matrix']}|min{r['minimum']}|composition9");st.to_csv(O/'tables/composition_statistics.tsv',sep='\t',index=False);pd.DataFrame(roster).to_csv(O/'tables/composition_roster.tsv',sep='\t',index=False);pd.DataFrame(matches).to_csv(O/'tables/matched_cell_differences.tsv',sep='\t',index=False);pd.DataFrame(matchroster).to_csv(O/'tables/matching_coverage.tsv',sep='\t',index=False);ms=pd.DataFrame(matchpat);ms.to_csv(O/'tables/matching_site_effects.tsv',sep='\t',index=False);md=ms.groupby(['matrix','caliper','donor']).effect.mean().reset_index();md.to_csv(O/'tables/matching_patient_effects.tsv',sep='\t',index=False);summarize(md,['matrix','caliper'],'effect',lambda r:'processing_matching_all_variants').to_csv(O/'tables/matching_statistics.tsv',sep='\t',index=False)
# Spatial inputs and scores.
sys.path.insert(0,str(P/'phase46_spatial_functional_20261001/code'));from common import feature_map,corr
base=pd.read_csv(P/'phase65_spatial_structure_20261005/tables/spatial_gene_and_program_scores.tsv.gz',sep='\t');sp=[];cov=[]
for sm,t in base.groupby('sample',sort=True):
 path=P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad';a=ad.read_h5ad(path);assert a.obs_names.tolist()==t.spot_id.tolist();fm=feature_map(a);den=np.asarray(a.X.sum(1)).ravel();assert np.array_equal(den,t.library_UMI);t=t[['spot_id','sample','patient_id','site_group','event_month','library_UMI','x','y','valid_library']].copy();present=[g for g in genes if g in fm];v=np.log1p(a.X[:,[fm[g] for g in present]].toarray()/np.where(den>0,den,np.nan)[:,None]*1e4)
 for j,g in enumerate(present):t['gene__'+g]=v[:,j]
 for name,gs in dict(panels,**contexts).items():
  ok=set(gs)<=set(present);t[name]=t[['gene__'+g for g in gs]].mean(1) if ok else np.nan;cov.append({'sample':sm,'panel':name,'complete':ok,'missing':';'.join(set(gs)-set(present))})
 t['balance']=t.III_V-t.IV;t['log_library']=np.log1p(t.library_UMI);sp.append(t);print('SPATIAL_EXTRACT',sm,flush=True)
sp=pd.concat(sp,ignore_index=True);assert len(sp)==53309 and sp.patient_id.nunique()==19;old=base.set_index('spot_id');assert np.allclose(sp.III_V,old.loc[sp.spot_id,'drop_dominant_pair__fibrillar_I_III_V'],equal_nan=True);assert np.allclose(sp.IV,old.loc[sp.spot_id,'full__collagen_IV'],equal_nan=True);sp.to_csv(O/'tables/spatial_scores.tsv.gz',sep='\t',index=False);pd.DataFrame(cov).to_csv(O/'tables/spatial_coverage.tsv',sep='\t',index=False)
basecols=['log_library','fibro_context','mural_context','endothelial_context','osteogenic_context','matrix_CAF'];sampleassoc=[]
for sm,z in sp.groupby('sample'):
 for ctl in ['broad_CAF','broad_CAF_plus_C1']:
  cols=basecols+['x','y']+(['C1_proxy'] if ctl.endswith('plus_C1') else [])
  for endpoint in ['III_V','IV','balance']:
   rho=corr(z.C3_proxy,z[endpoint],z[cols].to_numpy());sampleassoc.append({'sample':sm,'patient_id':z.patient_id.iloc[0],'site_group':z.site_group.iloc[0],'event_month':z.event_month.iloc[0],'control':ctl,'endpoint':endpoint,'rho':rho})
sampleassoc=pd.DataFrame(sampleassoc);sampleassoc.to_csv(O/'tables/spatial_section_associations.tsv',sep='\t',index=False);pa=[]
for scope,sub in [('all_patients',sampleassoc),('primary_only',sampleassoc[sampleassoc.site_group=='primary_tumor'])]:
 keys=['control','endpoint'];events=sub.groupby(['patient_id','site_group','event_month']+keys,dropna=False).rho.mean().reset_index();sites=events.groupby(['patient_id','site_group']+keys,dropna=False).rho.mean().reset_index();patients=sites.groupby(['patient_id']+keys,dropna=False).rho.mean().reset_index();patients['scope']=scope;pa.append(patients)
pa=pd.concat(pa);pa.to_csv(O/'tables/spatial_patient_associations.tsv',sep='\t',index=False);summarize(pa,['scope','control','endpoint'],'rho',lambda r:'spatial_partial_association12').to_csv(O/'tables/spatial_association_statistics.tsv',sep='\t',index=False)
# Leave whole patients out, equal weight per patient/site/event/section; no CV P-values.
def weights(t):
 keys=['patient_id','site_group','event_month','sample'];a=t.groupby(keys,dropna=False).size().rename('n').reset_index();a['sections']=a.groupby(keys[:3],dropna=False)['sample'].transform('nunique');e=a[keys[:3]].drop_duplicates();e['events']=e.groupby(keys[:2],dropna=False).event_month.transform('size');si=e[keys[:2]].drop_duplicates();si['sites']=si.groupby('patient_id').site_group.transform('size');a=a.merge(e,on=keys[:3]).merge(si,on=keys[:2]);a['weight']=1/(a.n*a.sections*a.events*a.sites);return t[keys].merge(a[keys+['weight']],on=keys,how='left',sort=False).weight.to_numpy()
cvrows=[];predictions=[];models=[]
for scope,sub in [('all_patients',sp),('primary_only',sp[sp.site_group=='primary_tumor'])]:
 sub=sub[sub.valid_library].copy();required=basecols+['C3_proxy','C1_proxy','III_V','IV','balance'];sub=sub.dropna(subset=required);w=weights(sub)
 for patient in sorted(sub.patient_id.unique()):
  ho=sub.patient_id.to_numpy()==patient;tr=~ho;wt=w[tr];wh=w[ho];wt=wt/np.sum(wt)*sub.loc[tr,'patient_id'].nunique();wh=wh/wh.sum();outcomes=sub[['III_V','IV','balance']].to_numpy();Y=outcomes[tr];Yh=outcomes[ho]
  for control in ['broad_CAF','broad_CAF_plus_C1']:
   cols=basecols+(['C1_proxy'] if control.endswith('plus_C1') else []);pp=[]
   for added in [False,True]:
    cc=cols+(['C3_proxy'] if added else []);X=sub[cc].to_numpy();mu=np.average(X[tr],axis=0,weights=wt);sd=np.maximum(np.sqrt(np.average((X[tr]-mu)**2,axis=0,weights=wt)),.001);Xt=np.column_stack([np.ones(tr.sum()),(X[tr]-mu)/sd]);Xh=np.column_stack([np.ones(ho.sum()),(X[ho]-mu)/sd]);pen=np.eye(Xt.shape[1]);pen[0,0]=0;coef=np.linalg.solve(Xt.T@(Xt*wt[:,None])+pen,Xt.T@(Y*wt[:,None]));pred=Xh@coef;pp.append(pred)
    for k,endpoint in enumerate(['III_V','IV','balance']):models.append({'scope':scope,'heldout_patient':patient,'control':control,'added_C3':added,'endpoint':endpoint,'terms':'intercept;'+';'.join(cc),'coefficients':json.dumps(coef[:,k].tolist()),'training_means':json.dumps(mu.tolist()),'training_sds':json.dumps(sd.tolist())})
   for k,endpoint in enumerate(['III_V','IV','balance']):
    yy=Yh[:,k];mean=np.sum(wh*yy);var=np.sum(wh*(yy-mean)**2);b=np.sum(wh*(yy-pp[0][:,k])**2);ext=np.sum(wh*(yy-pp[1][:,k])**2);cvrows.append({'scope':scope,'patient_id':patient,'control':control,'endpoint':endpoint,'spots':int(ho.sum()),'baseline_MSE':b,'extended_MSE':ext,'outcome_variance':var,'delta_R2':(b-ext)/var if var>1e-12 else np.nan});predictions.append(pd.DataFrame({'scope':scope,'patient_id':patient,'control':control,'endpoint':endpoint,'spot_id':sub.loc[ho,'spot_id'].to_numpy(),'weight':wh,'observed':yy,'baseline_prediction':pp[0][:,k],'extended_prediction':pp[1][:,k]}))
  print('CV_PATIENT',scope,patient,flush=True)
cv=pd.DataFrame(cvrows);cv.to_csv(O/'tables/spatial_LOPO_prediction_effects.tsv',sep='\t',index=False);pd.concat(predictions).to_csv(O/'tables/spatial_LOPO_predictions.tsv.gz',sep='\t',index=False);pd.DataFrame(models).to_csv(O/'resources/spatial_LOPO_models.tsv',sep='\t',index=False)
rr=[]
for k,z in cv.groupby(['scope','control','endpoint']):
 v=z.delta_R2.dropna();rr.append(dict(zip(['scope','control','endpoint'],k),N=len(v),mean_delta_R2=v.mean(),median_delta_R2=v.median(),positive=int((v>0).sum()),minimum=v.min(),maximum=v.max(),P_not_computed='Overlapping CV training folds; descriptive predictive gain only'))
pd.DataFrame(rr).to_csv(O/'tables/spatial_prediction_summary.tsv',sep='\t',index=False);result={'status':'completed','cell_records':len(ce),'composition_tests':len(st),'matching_cell_records':len(matches),'spatial_spots':len(sp),'spatial_patients':sp.patient_id.nunique(),'spatial_prior_scores_reproduced':True,'spatial_CV_rows':len(cv),'panels_do_not_share_outcome_predictor_genes':not bool(set(III+IV)&set(sum([panels['C3_proxy'],panels['C1_proxy'],panels['matrix_CAF']]+list(contexts.values()),[])))};assert result['panels_do_not_share_outcome_predictor_genes'];(O/'COMPLETE.json').write_text(json.dumps(result,indent=2));print('COMPLETE',json.dumps(result),flush=True)
