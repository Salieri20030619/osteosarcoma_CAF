from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,itertools
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=P/'phase17_all_cohort_expanded_20260921';O=P/'phase18_layer_evidence_20260921'
for f in ['code','tables','delivery','logs']:(O/f).mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
a=ad.read_h5ad(I/'objects/expanded_integrated.h5ad')
surfaces=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates'];genes=surfaces+['PTPRC','PECAM1','EPCAM']
sy=a.var.symbol.astype(str).to_numpy();x=a.layers['counts'].tocsc();v=np.column_stack([np.asarray(x[:,sy==g].sum(1)).ravel() for g in genes]);assert all((sy==g).any() for g in genes)
dump=(v[:,-3:]<1).all(1);pairs=[(g,) for g in surfaces]+list(itertools.combinations(surfaces,2));gate={'+'.join(pair):dump&(v[:,[genes.index(g) for g in pair]]>=1).all(1) for pair in pairs}
base=read(I/'tables/surface_screen_all_QC.tsv.gz')
sms=sorted(a.obs['sample'].astype(str).unique());smcodes=pd.Categorical(a.obs['sample'],categories=sms).codes
base=base.set_index(['sample','markers']);fulltarget=a.obs.myCAF_customer_positive.to_numpy()
for name,gg in gate.items():
 counts=np.bincount(smcodes[gg&fulltarget],minlength=len(sms))
 for sm,n in zip(sms,counts):assert int(base.loc[(sm,name),'target_in_gate'])==int(n),(sm,name,n)
print('FULL_SOURCE_SURFACE_NUMERATOR_CONCORDANCE_PASS',flush=True)
markers=[];effects=[];screens=[];panelranks=[];D={}
for layer in ['strict','supported','expanded']:
 d=read(I/'tables'/f'{layer}_cell_labels.tsv.gz');d['cluster']=d.cluster.astype(str)
 d['known_patient']=~d.patient_id.fillna('').astype(str).str.lower().isin(['','nan','na','unknown','unconfirmed','not_established'])
 d['analysis_donor']=np.where(d.known_patient,d.dataset+'|'+d.patient_id,d.dataset+'|'+d['sample']+'|UNKNOWN_PATIENT')
 ix=a.obs_names.get_indexer(d.cell_id);assert (ix>=0).all();z=a.X[ix].tocsr();D[layer]=d
 for cl in sorted(d.cluster.unique(),key=int):
  mask=d.cluster.eq(cl).to_numpy();ins=np.asarray(z[mask].mean(0)).ravel();outs=np.asarray(z[~mask].mean(0)).ravel();delta=ins-outs
  top=np.argsort(delta)[-20:][::-1]
  top=[j for j in top if not str(sy[j]).startswith(('MT-','RPL','RPS'))][:15]
  for j in top:markers.append({'layer':layer,'cluster':cl,'gene_id':a.var_names[j],'symbol':sy[j],'mean_inside':ins[j],'mean_outside':outs[j],'delta_mean_logCP10k':delta[j],'interpretation':'cell-weighted discovery rank, descriptive, not a DE test'})
  for donor,di in d[d.known_patient].groupby('analysis_donor').groups.items():
   di=np.array(list(di));mi=di[mask[di]];mo=di[~mask[di]]
   if len(mi)<10 or len(mo)<10:continue
   dd=np.asarray(z[mi][:,top].mean(0)-z[mo][:,top].mean(0)).ravel()
   for j,value in zip(top,dd):effects.append({'layer':layer,'cluster':cl,'donor_record':donor,'dataset':d.loc[di[0],'dataset'],'gene_id':a.var_names[j],'symbol':sy[j],'delta':value,'inside_cells':len(mi),'outside_cells':len(mo)})
  target=np.zeros(a.n_obs,bool);target[ix[mask&d.myCAF_customer_positive.to_numpy()]]=True
  targets=np.bincount(smcodes[target],minlength=len(sms));rawrows=[]
  for name,gg in gate.items():
   counts=np.bincount(smcodes[target&gg],minlength=len(sms))
   for sm,n,nt in zip(sms,counts,targets):
    b=base.loc[(sm,name)];rawrows.append({'sample':sm,'markers':name,'target_cells':int(nt),'target_in_gate':int(n),'gate_cells':int(b.gate_cells),'QC_cells':int(b.QC_cells),'CNV_high_in_gate':int(b.CNV_high_in_gate),'mural_in_gate':int(b.mural_in_gate)})
  rr=pd.DataFrame(rawrows);md=a.obs[['sample','dataset','patient_id','site']].drop_duplicates('sample').copy();md['known_patient']=~md.patient_id.astype(str).str.lower().isin(['unconfirmed','unknown','nan','','na']);md['donor_record']=md.dataset.astype(str)+'|'+md.patient_id.astype(str)
  md['site_group']=np.where(md.site.astype(str).str.contains('lung',case=False),'lung',np.where(md.site.astype(str).str.contains('bone|femur|tibia|humerus|fibula|pelvis|ilium|sacrum|radius|ulna',case=False),'bone','other_or_unknown'))
  rr=rr.merge(md,on='sample',validate='many_to_one');rr=rr[rr.known_patient]
  sd=rr.groupby(['donor_record','site_group','markers'])[['target_cells','target_in_gate','gate_cells','QC_cells','CNV_high_in_gate','mural_in_gate']].sum().reset_index()
  sd['coverage']=sd.target_in_gate/sd.target_cells.replace(0,np.nan);sd['known_target_fraction']=sd.target_in_gate/sd.gate_cells.replace(0,np.nan);sd['CNV_high_fraction']=sd.CNV_high_in_gate/sd.gate_cells.replace(0,np.nan);sd['mural_fraction']=sd.mural_in_gate/sd.gate_cells.replace(0,np.nan)
  sd['layer']=layer;sd['cluster']=cl;screens.append(sd)
  for name,g in sd[sd.target_cells>=10].groupby('markers'):
   # donor/site records are observations, not automatically independent replicates.
   panelranks.append({'layer':layer,'cluster':cl,'markers':name,'evaluable_donor_site_records':len(g),'known_donor_records':g.donor_record.nunique(),'median_coverage':g.coverage.median(),'median_known_target_fraction':g.known_target_fraction.median(),'median_CNV_high_fraction':g.CNV_high_fraction.median(),'median_mural_fraction':g.mural_fraction.median(),'ranking_score':(g.coverage*g.known_target_fraction.fillna(0)).median(),'evidence_level':'S1_RNA_only','protein_validation':False})
  print('CLUSTER_DONE',layer,cl,flush=True)
save(pd.DataFrame(markers),'descriptive_gene_candidates.tsv');E=pd.DataFrame(effects);save(E,'gene_candidate_donor_effects.tsv')
gr=[]
for key,g in E.groupby(['layer','cluster','gene_id','symbol']):gr.append(dict(zip(['layer','cluster','gene_id','symbol'],key),evaluable_donor_records=len(g),datasets=g.dataset.nunique(),positive_records=int((g.delta>0).sum()),median_effect=g.delta.median(),positive_fraction=float((g.delta>0).mean()),evidence_level='G2_exploratory_recurrent' if len(g)>=3 and g.dataset.nunique()>=2 and (g.delta>0).mean()>=.75 and g.delta.median()>=.1 else 'G1_discovery_only'))
save(pd.DataFrame(gr),'gene_candidate_evidence.tsv');save(pd.concat(screens),'cluster_surface_donor_screen.tsv.gz');save(pd.DataFrame(panelranks).sort_values(['layer','cluster','ranking_score'],ascending=[True,True,False]),'cluster_surface_ranking.tsv')
qa={'source_surface_numerators_match_full_source':True,'layers':3,'clusters':sum(d.cluster.nunique() for d in D.values()),'surface_combinations':136,'surface_denominator':'all QC cells of 51 samples contributing expanded candidates; zero-target samples within this scope retained','unknown_patients_excluded_from_donor_summaries':True,'gene_candidates':len(markers),'marker_claims':'RNA only, no protein or metastasis specificity validation','source_objects_unchanged':True}
(O/'delivery/QA_remote.json').write_text(json.dumps(qa,indent=2));(O/'REMOTE_COMPLETE').write_text('MARKER_AND_SURFACE_EXPLORATION_COMPLETE\n');print(json.dumps(qa),flush=True)
