"""SYNTAX_CHECK: python -m py_compile code/integrate_identity.py
RUN_COMMAND: python code/integrate_identity.py
EXPECTED_SUCCESS_MARKER: PHASE6_IDENTITY_INTEGRATION_COMPLETE
"""
from pathlib import Path
import json,datetime
import pandas as pd,numpy as np,anndata as ad
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=r/'phase6_nonpdx_cnv_20260916/pilot_integration_check_20260916';out=p/'identity';(out/'objects').mkdir(exist_ok=True)
manifest=pd.read_csv(p/'review/input_summary.tsv',sep='\t');final=[];summaries=[];toolqc=[];clusters=[];objects={}
for row in manifest.to_dict('records'):
 sample=row['sample'];o=pd.read_csv(p/'review'/(sample+'_cell_evidence.tsv.gz'),sep='\t').set_index('cell_id',drop=False)
 qualities={};states={}
 for method,normal in [('infercnv','low_signal'),('copykat','diploid'),('scevan','normal')]:
  loc=p/method/sample;completed=(loc/'COMPLETE').exists()
  ec=loc/'job.exit_code';states[method]='complete' if completed else ('failed' if ec.exists() and ec.read_text().strip()!='0' else 'not_run')
  if completed:
   d=pd.read_csv(loc/'cell_results.tsv',sep='\t').set_index('cell_id')
   assert d.index.is_unique and set(d.index)<=set(o.index)
   o[method+'_call']=d.call.reindex(o.index).fillna('not_run');o[method+'_signal']=d.signal_strength.reindex(o.index)
   used=sorted(d.loc[d.cnv_role.str.startswith('ref_'),'cnv_role'].unique());metrics={};ok=bool(used)
   for role in used+[x.replace('ref_','holdout_') for x in used]:
    g=d[(d.cnv_role==role)&(d.call!='not_evaluable')];frac=float((g.call==normal).mean()) if len(g) else 0
    metrics[role]={'evaluable':len(g),'normal_fraction':frac};ok=ok and len(g)>=10 and frac>=.9
   held=d[d.cnv_role.str.startswith('holdout_')&(d.call!='not_evaluable')]
   ok=ok and len(held)>=10
   toolqc.append({'sample':sample,'method':method,'state':states[method],'reference_check_pass':bool(ok),'heldout_evaluable':len(held),'heldout_normal_fraction':float((held.call==normal).mean()) if len(held) else 0,'lineage_checks':json.dumps(metrics),'note':'Reference and calibration consistency, not proof of stromal specificity'})
   qualities[method]=bool(ok)
  else:
   o[method+'_call']=np.where(o.cnv_role=='not_selected','not_run',states[method]);o[method+'_signal']=np.nan;qualities[method]=False
   toolqc.append({'sample':sample,'method':method,'state':states[method],'reference_check_pass':False,'note':'Missing/failed method is not CNV negative'})
 secondary='copykat' if qualities['copykat'] else ('scevan' if qualities['scevan'] else 'copykat')
 o['decision_cnv2_source']=secondary;o['decision_cnv2_call']=o[secondary+'_call']
 ic=o.infercnv_call;ck=o[secondary+'_call'].replace({'normal':'diploid','tumor':'aneuploid'})
 low=(ic=='low_signal')&(ck=='diploid');high=(ic=='elevated_signal')&(ck=='aneuploid')
 evaluated=ic.isin(['low_signal','elevated_signal','ambiguous_signal'])&ck.isin(['diploid','aneuploid'])
 o['cnv_concordance']='not_jointly_evaluable';o.loc[low,'cnv_concordance']='both_low';o.loc[high,'cnv_concordance']='both_high';o.loc[evaluated&~low&~high,'cnv_concordance']='discordant_or_ambiguous'
 o['cluster_high_fraction']=np.nan
 for c,z in o.groupby('cluster'):
  ix=z.index;ev=evaluated.loc[ix];frac=float(high.loc[ix][ev].mean()) if ev.any() else np.nan
  o.loc[ix,'cluster_high_fraction']=frac
  clusters.append({'sample':sample,'cluster':c,'lineage':z.lineage_reviewed.iloc[0],'cells':len(z),'jointly_evaluable':int(ev.sum()),'both_low':int(low.loc[ix].sum()),'both_high':int(high.loc[ix].sum()),'concordant_high_fraction':frac})
 good=qualities['infercnv'] and qualities[secondary];o['reference_QC_pass']=good
 o['infercnv_reference_QC_pass']=qualities['infercnv'];o['copykat_reference_QC_pass']=qualities['copykat']
 immune_pct=sum(o['pctUMI_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC'])
 competing=((o.det_Mural_review>=3)&(o.pctUMI_Mural>.2))|((o.det_OC_review>=2)&(o.pctUMI_OC_specific>.1))|((o.det_T_NK_review>=3)&(o.pctUMI_T_NK>.2))|((o.det_Myeloid_review>=4)&(o.pctUMI_Myeloid>.2))|((o.det_Endo_review>=3)&(o.pctUMI_Endo>.2))|((o.det_Epi_review>=3)&(o.pctUMI_Epi>.2))|((o.det_Myogenic_review>=2)&(o.pctUMI_Myogenic>.2))
 clean=~o.mixed_lineage_fraction_flag&~o.CTA_context&(o.det_CTA_review<2)
 nonmalcluster=o.cluster_high_fraction.fillna(1)<.2
 fib=(o.lineage_reviewed=='FibLike')&(o.det_Fib_review>=2)
 mesfib=(o.lineage_reviewed=='MesAmb')&(o.det_Fib_review>=3)&(o.det_Osteo_review<2)&(o.det_Mural_review<2)
 o['caf_strict']=fib&low&clean&~competing&nonmalcluster&(o.n_genes>=500)&good
 one_good_low=((ic=='low_signal')&qualities['infercnv'])|((ck=='diploid')&qualities[secondary])
 no_opposing=(ic!='elevated_signal')&(ck!='aneuploid')
 # A single-method sensitivity set still requires assessable cluster risk from a valid method.
 risk=o.cluster_high_fraction.copy()
 if qualities['infercnv']:
  for c,z in o.groupby('cluster'):
   if risk.loc[z.index].isna().all():
    calls=z.infercnv_call;ev=calls.isin(['low_signal','elevated_signal','ambiguous_signal'])
    if ev.any():risk.loc[z.index]=float((calls[ev]=='elevated_signal').mean())
 elif qualities[secondary]:
  for c,z in o.groupby('cluster'):
   if risk.loc[z.index].isna().all():
    calls=z[secondary+'_call'].replace({'normal':'diploid','tumor':'aneuploid'});ev=calls.isin(['diploid','aneuploid'])
    if ev.any():risk.loc[z.index]=float((calls[ev]=='aneuploid').mean())
 o['caf_extended']=o.caf_strict|((fib|mesfib)&one_good_low&no_opposing&clean&~competing&(risk.fillna(1)<.2))
 o['final_identity']='Identity_pending'
 other=o.lineage_reviewed.isin(['T_NK','B_Plasma','Myeloid','DC','Mast','Erythroid','Endo','Epi','OC_like','Immune_ambiguous'])&~high&~o.mixed_lineage_fraction_flag
 o.loc[other,'final_identity']='Other'
 o.loc[(o.lineage_reviewed=='Mural')&low&clean&good,'final_identity']='Mural_SMC'
 malignant=high&o.lineage_reviewed.isin(['FibLike','MesAmb','OsteoLike','Cycling_pending','Chondrogenic_pending','Myogenic_pending'])&~o.mixed_lineage_fraction_flag&good
 o.loc[malignant,'final_identity']='Malignant_OS_candidate'
 o.loc[o.caf_strict|(mesfib&low&o.caf_extended&good),'final_identity']='Fib_candidate'
 o['unresolved_issues']=''
 for mask,label in [(o.cnv_concordance=='discordant_or_ambiguous','CNV_discordance;'),(o.cnv_concordance=='not_jointly_evaluable','CNV_not_jointly_evaluable;'),(o.mixed_lineage_fraction_flag,'competing_lineage;'),(o.CTA_context,'CTA_cluster;'),(o.n_genes<500,'low_complexity;')]:o.loc[mask,'unresolved_issues']+=label
 if not good:o['unresolved_issues']+='reference_QC_not_passed;'
 o['CAF_note']='tumor-associated fibroblast candidates; no metastatic-specific state claim'
 o.to_csv(out/(sample+'_identity.tsv.gz'),sep='\t',index=False);final.append(o)
 summaries.append({'sample':sample,'short':row['short'],'patient_id':row['patient_id'],'site':row['site'],'clinical_group':row['clinical_group'],'retained':len(o),'cnv_input':int((o.cnv_role!='not_selected').sum()),'infercnv_state':states['infercnv'],'copykat_state':states['copykat'],'scevan_state':states['scevan'],'selected_secondary_method':secondary,'reference_QC_pass':good,'caf_strict':int(o.caf_strict.sum()),'caf_extended':int(o.caf_extended.sum()),'CNV_discordant_or_ambiguous':int((o.cnv_concordance=='discordant_or_ambiguous').sum()),**{k:int(v) for k,v in o.final_identity.value_counts().items()}})
 if o.caf_extended.any():
  a=ad.read_h5ad(r/'phase5_nonpdx_qc_20260916/results'/sample/'annotated.h5ad');assert a.obs_names.tolist()==o.index.tolist()
  for flag in ['caf_strict','caf_extended']:
   if not o[flag].any():continue
   ix=o[flag].to_numpy();obj=a[ix].copy();obs=o.loc[ix].drop(columns='cell_id').copy()
   for col in obs.select_dtypes(include='object'):obs[col]=obs[col].fillna('').astype(str)
   obj.obs=obs;objects.setdefault((str(o.dataset.iloc[0]),flag),[]).append(obj)
  del a
 print('IDENTITY_INTEGRATED',sample,summaries[-1]['caf_strict'],summaries[-1]['caf_extended'],flush=True)
allcells=pd.concat(final,ignore_index=True);assert len(allcells)==int(manifest.retained.sum()) and allcells.cell_id.is_unique
allcells.to_csv(out/'all_cell_identity_evidence.tsv.gz',sep='\t',index=False)
for flag in ['caf_strict','caf_extended']:
 allcells.loc[allcells[flag],['dataset','sample','short','patient_id','site','group','cell_id','barcode','final_identity','unresolved_issues']].to_csv(out/(flag+'_cells.tsv'),sep='\t',index=False)
for (dataset,flag),objs in objects.items():
 obj=ad.concat(objs,join='inner',merge='first');obj.uns['gene_join']='within-dataset shared Ensembl IDs';obj.uns['identity_status']='candidate sets; CAF state analysis pending';obj.write_h5ad(out/'objects'/(dataset+'_'+flag+'.h5ad'),compression='gzip')
pd.DataFrame(summaries).fillna(0).to_csv(out/'sample_identity_summary.tsv',sep='\t',index=False)
pd.DataFrame(toolqc).to_csv(out/'tool_reference_QC.tsv',sep='\t',index=False)
pd.DataFrame(clusters).to_csv(out/'cluster_CNV_concordance.tsv',sep='\t',index=False)
state={'state':'INITIAL_DUAL_CNV_INTEGRATION_COMPLETE','completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(summaries),'strict_CAF_candidates':int(allcells.caf_strict.sum()),'extended_CAF_candidates':int(allcells.caf_extended.sum()),'samples_passing_both_reference_checks':sum(x['reference_QC_pass'] for x in summaries),'failed_or_missing_tool_jobs':sum(x['state']!='complete' for x in toolqc if x['method'] in ['infercnv','copykat']),'next_action':'review discordant/failed methods and CAF sets before state analysis; no automatic biological conclusion'}
(out/'integration_parameters.json').write_text(json.dumps({'reference_rules_file':'resources/reference_protocol.json','min_reference_and_holdout_per_used_lineage':10,'min_normal_fraction':.9,'max_cluster_concordant_high_fraction':.2,'strict_min_genes':500,'competing_lineage_UMI_percent':.2,'OC_specific_UMI_percent':.1,'no_patient_group_inference':True},indent=2))
(p/'STATUS.json').write_text(json.dumps(state,indent=2));(out/'COMPLETE').write_text('PHASE6_IDENTITY_INTEGRATION_COMPLETE\n');print(json.dumps(state),flush=True)
