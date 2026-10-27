from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,json,ast,time
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase79_C3_subdivision_value_20261010';src=(O/'code/analyze_subdivision.py').read_text(encoding='utf-8-sig');tree=ast.parse(src)
for node in tree.body:
 if isinstance(node,ast.FunctionDef) and node.name in ['exact','summarize','weights']:exec(compile(ast.Module(body=[node],type_ignores=[]),'<frozen_helpers>','exec'))
sys.path.insert(0,str(P/'phase46_spatial_functional_20261001/code'));from common import corr
amendment={'at_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'trigger':'All spatial sections do not assay RUNX2 and MYH11; strict original panels produce0evaluable patients; first run exited1 while trying to concatenate empty predictions','original_complete_panel_result':'Unevaluable, retained in original spatial_association_statistics.tsv and original protocol/coverage/log','additional_analysis':'Assay-limited sensitivity, defined before any evaluable spatial outcome calculation','osteogenic_available5':['SP7','ALPL','BGLAP','IBSP','IFITM5'],'contractile_available2':['ACTA2','TAGLN'],'rule':'Use measured genes explicitly; do not fill missing RUNX2/MYH11 with zeros; osteogenic5 matches the preexisting Phase65 measured-gene context','primary_control':'depth+four measured contexts+matrixCAF6','sensitivity_control':'primary control+contractile_available2','inference':'same outcomes and testing procedures as frozen protocol; one separate12-test family for this sensitivity; full-panel spatial confirmation remains unavailable','prediction':'All and primary separately; entire patient held out; no P on CV fold gains'};(O/'resources/SPATIAL_ASSAY_AMENDMENT.json').write_text(json.dumps(amendment,indent=2))
sp=pd.read_csv(O/'tables/spatial_scores.tsv.gz',sep='\t',dtype={'event_month':str});sp['osteogenic_available5']=sp[['gene__'+g for g in amendment['osteogenic_available5']]].mean(1);sp['contractile_available2']=sp[['gene__ACTA2','gene__TAGLN']].mean(1);old=pd.read_csv(P/'phase65_spatial_structure_20261005/tables/spatial_gene_and_program_scores.tsv.gz',sep='\t',usecols=['spot_id','context__osteogenic_context']).set_index('spot_id');assert np.allclose(sp.osteogenic_available5,old.loc[sp.spot_id,'context__osteogenic_context'],equal_nan=True);sp.to_csv(O/'tables/spatial_scores_assay_limited.tsv.gz',sep='\t',index=False)
basecols=['log_library','fibro_context','mural_context','endothelial_context','osteogenic_available5','matrix_CAF'];sampleassoc=[];cnames=['broad_CAF_available','broad_CAF_plus_contractile2']
for sm,z in sp.groupby('sample'):
 for ctl in cnames:
  cols=basecols+['x','y']+(['contractile_available2'] if ctl.endswith('contractile2') else [])
  for endpoint in ['III_V','IV','balance']:
   rho=corr(z.C3_proxy,z[endpoint],z[cols].to_numpy());sampleassoc.append({'sample':sm,'patient_id':z.patient_id.iloc[0],'site_group':z.site_group.iloc[0],'event_month':z.event_month.iloc[0],'control':ctl,'endpoint':endpoint,'rho':rho})
sampleassoc=pd.DataFrame(sampleassoc);sampleassoc.to_csv(O/'tables/assay_limited_spatial_section_associations.tsv',sep='\t',index=False);pa=[]
for scope,sub in [('all_patients',sampleassoc),('primary_only',sampleassoc[sampleassoc.site_group=='primary_tumor'])]:
 keys=['control','endpoint'];events=sub.groupby(['patient_id','site_group','event_month']+keys,dropna=False).rho.mean().reset_index();sites=events.groupby(['patient_id','site_group']+keys,dropna=False).rho.mean().reset_index();patients=sites.groupby(['patient_id']+keys,dropna=False).rho.mean().reset_index();patients['scope']=scope;pa.append(patients)
pa=pd.concat(pa);pa.to_csv(O/'tables/assay_limited_spatial_patient_associations.tsv',sep='\t',index=False);summarize(pa,['scope','control','endpoint'],'rho',lambda r:'assay_limited_spatial_association12').to_csv(O/'tables/assay_limited_spatial_association_statistics.tsv',sep='\t',index=False)
cvrows=[];predictions=[];models=[]
for scope,sub in [('all_patients',sp),('primary_only',sp[sp.site_group=='primary_tumor'])]:
 sub=sub[sub.valid_library].copy();required=basecols+['C3_proxy','contractile_available2','III_V','IV','balance'];sub=sub.dropna(subset=required);w=weights(sub);assert sub.patient_id.nunique()==(19 if scope=='all_patients' else 9)
 for patient in sorted(sub.patient_id.unique()):
  ho=sub.patient_id.to_numpy()==patient;tr=~ho;wt=w[tr];wh=w[ho];wt=wt/wt.sum()*sub.loc[tr,'patient_id'].nunique();wh=wh/wh.sum();outcomes=sub[['III_V','IV','balance']].to_numpy();Y=outcomes[tr];Yh=outcomes[ho]
  for control in cnames:
   cols=basecols+(['contractile_available2'] if control.endswith('contractile2') else []);pp=[]
   for added in [False,True]:
    cc=cols+(['C3_proxy'] if added else []);X=sub[cc].to_numpy();mu=np.average(X[tr],axis=0,weights=wt);sd=np.maximum(np.sqrt(np.average((X[tr]-mu)**2,axis=0,weights=wt)),.001);Xt=np.column_stack([np.ones(tr.sum()),(X[tr]-mu)/sd]);Xh=np.column_stack([np.ones(ho.sum()),(X[ho]-mu)/sd]);pen=np.eye(Xt.shape[1]);pen[0,0]=0;coef=np.linalg.solve(Xt.T@(Xt*wt[:,None])+pen,Xt.T@(Y*wt[:,None]));pred=Xh@coef;pp.append(pred)
    for k,endpoint in enumerate(['III_V','IV','balance']):models.append({'scope':scope,'heldout_patient':patient,'control':control,'added_C3':added,'endpoint':endpoint,'terms':'intercept;'+';'.join(cc),'coefficients':json.dumps(coef[:,k].tolist()),'training_means':json.dumps(mu.tolist()),'training_sds':json.dumps(sd.tolist())})
   for k,endpoint in enumerate(['III_V','IV','balance']):
    yy=Yh[:,k];mean=wh@yy;var=wh@((yy-mean)**2);b=wh@((yy-pp[0][:,k])**2);ext=wh@((yy-pp[1][:,k])**2);cvrows.append({'scope':scope,'patient_id':patient,'control':control,'endpoint':endpoint,'spots':int(ho.sum()),'baseline_MSE':b,'extended_MSE':ext,'outcome_variance':var,'delta_R2':(b-ext)/var if var>1e-12 else np.nan});predictions.append(pd.DataFrame({'scope':scope,'patient_id':patient,'control':control,'endpoint':endpoint,'spot_id':sub.loc[ho,'spot_id'].to_numpy(),'weight':wh,'observed':yy,'baseline_prediction':pp[0][:,k],'extended_prediction':pp[1][:,k]}))
  print('CV',scope,patient,flush=True)
cv=pd.DataFrame(cvrows);cv.to_csv(O/'tables/assay_limited_spatial_LOPO_prediction_effects.tsv',sep='\t',index=False);pd.concat(predictions).to_csv(O/'tables/assay_limited_spatial_LOPO_predictions.tsv.gz',sep='\t',index=False);pd.DataFrame(models).to_csv(O/'resources/assay_limited_spatial_LOPO_models.tsv',sep='\t',index=False);rr=[]
for k,z in cv.groupby(['scope','control','endpoint']):
 v=z.delta_R2.dropna();rr.append(dict(zip(['scope','control','endpoint'],k),N=len(v),mean_delta_R2=v.mean(),median_delta_R2=v.median(),positive=int((v>0).sum()),minimum=v.min(),maximum=v.max(),P_not_computed='Overlapping CV training folds; descriptive predictive gain only'))
pd.DataFrame(rr).to_csv(O/'tables/assay_limited_spatial_prediction_summary.tsv',sep='\t',index=False)
result={'status':'completed_with_recorded_assay_limitation','cell_composition_and_matching_completed':True,'strict_full_panel_spatial':'not_evaluable_RUNX2_MYH11_unassayed','assay_limited_spatial_completed':True,'spatial_spots':len(sp),'valid_library_spots':int(sp.valid_library.sum()),'spatial_patients':19,'primary_patients':9,'spatial_CV_rows':len(cv),'previous_spatial_context_reproduced':True,'no_master_edits':True};(O/'COMPLETE.json').write_text(json.dumps(result,indent=2));print('COMPLETE',json.dumps(result),flush=True)
