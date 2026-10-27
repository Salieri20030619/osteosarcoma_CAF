from pathlib import Path
import os,json,datetime,traceback,itertools
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,stats
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
B=P/'phase28b_CAF_refined_20260926';R=P/'phase27_identity_review_20260926';O=P/'phase29_myCAF_followup_20260926'
# Reuse only function definitions of the newly reviewed CAF implementation.
src=(B/'code/cluster.py').read_text().split("lock=(O/")[0]
ns={};exec(compile(src,str(B/'code/cluster.py'),'exec'),ns)
O=P/'phase29_myCAF_followup_20260926'
for d in ['code','tables','objects','figures','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
ns['O']=O
model=ns['model'];profiles=ns['profiles'];plots=ns['plots'];clean=ns['clean']
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(p,v):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str));t.replace(p)
def status(s,**kw):writej(O/'STATUS.json',{'status':s,'at':now(),**kw});print(s,kw,flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def signflip(d):
 d=np.asarray(d,float);d=d[np.isfinite(d)];n=len(d)
 if n<2:return np.nan,'not_tested_no_biological_replication'
 observed=abs(d.mean())
 if n<=15:
  signs=np.array(list(itertools.product([-1,1],repeat=n)),dtype=np.int8)
  return float(np.mean(np.abs(signs@d/n)>=observed-1e-12)),'exact_paired_sign_flip'
 rng=np.random.default_rng(20260926);null=np.abs(rng.choice([-1,1],size=(19999,n))@d/n)
 return float((1+(null>=observed-1e-12).sum())/20000),'paired_sign_flip_MC_20000'
def unpaired(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float)
 if len(a)<2 or len(b)<2:return np.nan,'not_tested_no_within_group_biological_replication'
 r=stats.permutation_test((a,b),lambda x,y:np.mean(x)-np.mean(y),permutation_type='independent',n_resamples=19999,alternative='two-sided',random_state=20260926)
 return float(r.pvalue),'independent_permutation_up_to_19999'
def apply_fdr(t):
 t=t.copy();t['FDR']=np.nan
 if len(t):
  for _,g in t.groupby('hypothesis_family',dropna=False):
   ix=g.index[g.P.notna()]
   if len(ix):t.loc[ix,'FDR']=multipletests(t.loc[ix,'P'],method='fdr_bh')[1]
 return t
rc=1
try:
 annotation=read(B/'tables/CAF_cluster_annotation_reviewed.tsv').fillna('')
 a=ad.read_h5ad(B/'objects/CAF_integrated.h5ad')
 amap=annotation.set_index('cluster');amap.index=amap.index.astype(str)
 assert set(a.obs.cluster.astype(str))==set(amap.index)
 a.obs['CAF_type']=a.obs.cluster.astype(str).map(amap.CAF_type)
 a.obs['identity_review']=a.obs.cluster.astype(str).map(amap.identity_review)
 main=a.obs.identity_review.eq('retain_CAF_candidate')
 my=main&a.obs.CAF_type.eq('myCAF_like')
 writej(O/'PROTOCOL.json',{'source':'phase28b reviewed CAF cluster annotation','myCAF_selection':'reviewed myCAF-like CAF clusters; no per-cell marker-number gate','CNV':'phase27 frozen native/provider evidence, no conflict/unassessed cells in main','myCAF_reclustering':'refit HVG, PCA, Harmony, neighbors and Leiden within myCAF; original global clusters not reused as subclusters','statistics':'patient-level within-study comparisons; exact paired sign flip up to N15, Monte Carlo beyond; independent label permutation where both groups have replication; BH per same-class hypothesis family','small_N':'exploratory; biological samples only; directions, effects and consistency emphasized','unknown_clinical':'no inferred A from bone site; B1_history_supported is not contemporaneous pairing; no cross-study independent N pooling','surface':'RNA detectability and target fractions, not measured FACS purity or protein validation','active_S525_recovery':'does not block this frozen analysis'})
 a.obs.to_csv(O/'tables/CAF_reviewed_labels.tsv.gz',sep='\t',index_label='cell_id')
 assert my.sum()>2,'No review-supported myCAF cluster available'
 b=a[my].copy();b.obs['parent_CAF_cluster']=b.obs.cluster.astype(str);b.X=b.layers['counts'].copy()
 b=model(b);profile=profiles(b,'myCAF');plots(b,profile)
 b.obs.to_csv(O/'tables/myCAF_cell_labels.tsv.gz',sep='\t',index_label='cell_id');clean(b);b.write_h5ad(O/'objects/myCAF_integrated.h5ad',compression='gzip')
 save(b.obs.groupby(['dataset','sample','donor_record','clinical_group','cluster'],observed=True).size().reset_index(name='cells'),'myCAF_sample_cluster_counts.tsv')
 panels=[c for c in b.obs if c.startswith('program_') and any(z in c for z in ['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer'])]
 programs=[];markerrows=[]
 # Within each cohort, compare subcluster with other myCAF cells in the same participant.
 for ds,sub in b.obs.groupby('dataset',observed=True):
  for cl in sorted(b.obs.cluster.unique(),key=lambda x:int(x)):
   donor_delta=[];effect_records=[]
   for dn,g in sub.groupby('donor_record',observed=True):
    target=g.cluster.eq(cl);rest=~target
    if not target.any() or not rest.any():continue
    row={'dataset':ds,'donor_record':dn,'cluster':cl,'target_cells':int(target.sum()),'other_myCAF_cells':int(rest.sum())}
    for panel in panels:row[panel]=float(g.loc[target,panel].mean()-g.loc[rest,panel].mean())
    donor_delta.append(row)
    it=b.obs.index.get_indexer(g.index[target]);ir=b.obs.index.get_indexer(g.index[rest])
    tc=np.asarray(b.layers['counts'][it].sum(0)).ravel();rcount=np.asarray(b.layers['counts'][ir].sum(0)).ravel()
    td=b.obs.iloc[it].full_library_counts.sum();rd=b.obs.iloc[ir].full_library_counts.sum()
    effect_records.append(np.log2(tc/td*1e4+.1)-np.log2(rcount/rd*1e4+.1))
   if not donor_delta:continue
   dd=pd.DataFrame(donor_delta);save(dd,f'paired_program_effects_{ds}_C{cl}.tsv')
   for panel in panels:
    vals=dd[panel].to_numpy();pv,method=signflip(vals)
    programs.append(dict(dataset=ds,cluster=cl,panel=panel,independent_participants=len(vals),median_within_patient_delta=float(np.median(vals)),mean_delta=float(np.mean(vals)),positive=int((vals>0).sum()),negative=int((vals<0).sum()),P=pv,test=method,hypothesis_family=ds+'|myCAF_subcluster_programs',evidence='exploratory_discovery'))
   delta=np.vstack(effect_records);med=np.median(delta,axis=0);pos=(delta>0).sum(0)
   # All genes retain descriptive patient-level effects; no cell-level P values are assigned.
   frame=pd.DataFrame({'dataset':ds,'cluster':cl,'ensembl':b.var_names,'symbol':b.var.symbol.to_numpy(),'paired_participants':len(delta),'median_delta_log2CP10k':med,'positive_participants':pos,'negative_participants':(delta<0).sum(0)})
   save(frame,f'all_gene_patient_effects_{ds}_C{cl}.tsv.gz')
   markerrows.extend(frame.sort_values('median_delta_log2CP10k',ascending=False).head(50).to_dict('records'))
 save(apply_fdr(pd.DataFrame(programs)),'myCAF_patient_program_statistics.tsv');save(pd.DataFrame(markerrows),'myCAF_patient_top_gene_effects.tsv')
 # Clinical contrasts: aggregate repeats to the same participant before statistical calculations.
 clinic=read(R/'tables/clinical_mapping_v3.tsv').set_index('sample');clinical_rows=[];paired_rows=[];comparisons=[]
 for scope,ob in [('CAF',a.obs.loc[main]),('myCAF',b.obs)]:
  for sm,g in ob.groupby('sample',observed=True):
   cr=clinic.loc[sm];row={'scope':scope,'sample':sm,'dataset':cr.dataset,'donor_record':cr.donor_record,'site':str(cr.site),'clinical_group':str(cr.clinical_group),'treatment':str(cr.treatment),'cells':len(g)}
   row.update({panel:float(g[panel].mean()) for panel in panels});clinical_rows.append(row)
 clinical_frame=pd.DataFrame(clinical_rows);save(clinical_frame,'clinical_sample_programs.tsv')
 pairs=read(R/'tables/bone_lung_pairs_v3.tsv')
 for scope in ['CAF','myCAF']:
  cf=clinical_frame[clinical_frame.scope==scope]
  for ds,dsf in cf.groupby('dataset'):
   pair_effect=[]
   for r in pairs[pairs.dataset==ds].to_dict('records'):
    q=dsf[dsf.donor_record==r['donor_record']]
    bone=q[q['sample'].isin(r['bone_samples'].split(';'))];lung=q[q['sample'].isin(r['lung_samples'].split(';'))]
    if bone.empty or lung.empty:continue
    rec=dict(scope=scope,dataset=ds,donor_record=r['donor_record'],bone_samples=';'.join(bone['sample']),lung_samples=';'.join(lung['sample']),pair_type=r['pair_type'])
    # equal sample weighting within participant/site, then equal participant inference
    for panel in panels:rec[panel]=float(lung[panel].mean()-bone[panel].mean())
    paired_rows.append(rec);pair_effect.append(rec)
   if pair_effect:
    q=pd.DataFrame(pair_effect)
    for panel in panels:
     vals=q[panel].to_numpy();pv,method=signflip(vals)
     comparisons.append(dict(scope=scope,dataset=ds,contrast='same_participant_lung_minus_bone_timing_unresolved',panel=panel,n_pairs=len(vals),n_group1=len(vals),n_group0=len(vals),median_difference=float(np.median(vals)),positive=int((vals>0).sum()),negative=int((vals<0).sum()),P=pv,test=method,hypothesis_family=ds+'|'+scope+'|paired_site_programs',interpretation='exploratory_site_difference_not_metastasis_causation'))
   group=dsf.clinical_group.str.lower();aa=group.str.startswith('a_');bb=group.str.contains('metastatic.*primary|metastatic_bone_primary',regex=True)&~group.str.contains('b1')
   for label,m1,m0 in [('reported_any_metastasis_bone_minus_reported_nonmetastasis_bone',bb,aa)]:
    q1=dsf[m1].groupby('donor_record')[panels].mean();q0=dsf[m0].groupby('donor_record')[panels].mean()
    # prevent same donor appearing as both unpaired groups
    overlap=set(q1.index)&set(q0.index);q1=q1.drop(index=list(overlap));q0=q0.drop(index=list(overlap))
    if q1.empty or q0.empty:continue
    for panel in panels:
     v1=q1[panel].to_numpy();v0=q0[panel].to_numpy();pv,method=unpaired(v1,v0)
     comparisons.append(dict(scope=scope,dataset=ds,contrast=label,panel=panel,n_pairs=0,n_group1=len(v1),n_group0=len(v0),median_difference=float(np.median(v1)-np.median(v0)),positive=np.nan,negative=np.nan,P=pv,test=method,hypothesis_family=ds+'|'+scope+'|any_metastasis_programs',interpretation='reported_any_metastasis_not_verified_lung_B1'))
 save(pd.DataFrame(paired_rows),'bone_lung_patient_program_effects.tsv');save(apply_fdr(pd.DataFrame(comparisons)),'clinical_program_statistics.tsv')
 # Treatment labels describe exposure; only within-cohort/site comparisons avoid pooling unrelated cohorts.
 treatment_rows=[]
 for scope in ['CAF','myCAF']:
  cf=clinical_frame[clinical_frame.scope==scope].copy()
  def treatment_label(x):
   x=str(x).lower()
   if x=='no' or 'no chemotherapy' in x or 'untreated' in x:return 'reported_untreated'
   if x=='yes' or 'chemotherapy' in x or ' times ' in x:return 'reported_treated'
   return 'unknown_or_surgery_only'
  cf['treatment_group']=cf.treatment.map(treatment_label)
  cf['site_group']=np.where(cf.site.str.contains('lung',case=False),'lung',np.where(cf.site.str.contains('bone|femur|tibia|humerus|ulna|calcaneus',case=False,regex=True),'bone','other'))
  for (ds,site),g in cf.groupby(['dataset','site_group']):
   t=g[g.treatment_group.eq('reported_treated')].groupby('donor_record')[panels].mean()
   u=g[g.treatment_group.eq('reported_untreated')].groupby('donor_record')[panels].mean()
   shared=set(t.index)&set(u.index);t=t.drop(index=list(shared));u=u.drop(index=list(shared))
   if t.empty or u.empty:continue
   for panel in panels:
    pv,method=unpaired(t[panel],u[panel])
    treatment_rows.append(dict(scope=scope,dataset=ds,site_group=site,panel=panel,n_treated=len(t),n_untreated=len(u),median_difference=float(t[panel].median()-u[panel].median()),P=pv,test=method,hypothesis_family=ds+'|'+site+'|'+scope+'|treatment_programs',interpretation='exploratory_association_not_treatment_effect'))
 save(apply_fdr(pd.DataFrame(treatment_rows)),'treatment_program_statistics.tsv')
 # Report per-participant/subcluster denominators and evidence-grade sensitivity without numerical cell-count entry gates.
 denom=b.obs.groupby('donor_record',observed=True).size()
 prop=b.obs.groupby(['dataset','donor_record','cluster'],observed=True).size().reset_index(name='cells')
 prop['myCAF_denominator']=prop.donor_record.map(denom);prop['fraction_within_myCAF']=prop.cells/prop.myCAF_denominator
 save(prop,'myCAF_participant_subcluster_proportions.tsv')
 sensitivity=b.obs.groupby(['dataset','donor_record','cluster','CNV_evidence','identity_overlap_flag'],observed=True)[panels].agg(['mean','count'])
 sensitivity.columns=['__'.join(c) for c in sensitivity.columns];save(sensitivity.reset_index(),'myCAF_CNV_and_osteogenic_overlap_sensitivity.tsv')

 # Surface detection against full called-cell background for every contributing sample.
 surfaces=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
 gates=[(x,) for x in surfaces]+list(itertools.combinations(surfaces,2));surface_rows=[]
 myids=set(b.obs_names);targetmap=b.obs.cluster.astype(str).to_dict();allcls=sorted(b.obs.cluster.unique(),key=lambda x:int(x))
 for sm in sorted(b.obs['sample'].unique()):
  bg=ad.read_h5ad(R/'objects'/f'{sm}_marker_counts.h5ad');measured=set(bg.uns['measured_genes']);arr=bg.X.toarray();genes=list(bg.var_names)
  dump=(arr[:,[genes.index(g) for g in ['PTPRC','PECAM1','EPCAM']]]==0).all(1)
  targets=np.array([targetmap.get(c,'not_myCAF') for c in bg.obs_names])
  for pair in gates:
   if not set(pair)<=measured:continue
   mask=dump&(arr[:,[genes.index(g) for g in pair]]>0).all(1)
   for cl in allcls:
    target=targets==cl;nt=int(target.sum())
    if not nt:continue
    surface_rows.append(dict(sample=sm,dataset=bg.obs.dataset.iloc[0],donor_record=bg.obs.donor_record.iloc[0],cluster=cl,markers='+'.join(pair),full_called_cells=len(bg),target_cells=nt,RNA_gate_cells=int(mask.sum()),target_in_RNA_gate=int((mask&target).sum()),CNV_alteration_or_conflict_in_gate=int((mask&bg.obs.CNV_evidence.isin(['method_conflict','single_method_alteration','concordant_alteration']).to_numpy()).sum()),mural_in_gate=int((mask&bg.obs.provisional_lineage.eq('Mural').to_numpy()).sum()),known_doublet_in_gate=int((mask&bg.obs.known_doublet.to_numpy()).sum())))
 st=pd.DataFrame(surface_rows);save(st,'surface_RNA_sample_screen.tsv.gz')
 dt=st.groupby(['dataset','donor_record','cluster','markers'],as_index=False)[['full_called_cells','target_cells','RNA_gate_cells','target_in_RNA_gate','CNV_alteration_or_conflict_in_gate','mural_in_gate','known_doublet_in_gate']].sum()
 dt['target_coverage']=dt.target_in_RNA_gate/dt.target_cells
 dt['target_fraction_of_RNA_gate']=dt.target_in_RNA_gate/dt.RNA_gate_cells.replace(0,np.nan)
 dt['CNV_risk_fraction_of_RNA_gate']=dt.CNV_alteration_or_conflict_in_gate/dt.RNA_gate_cells.replace(0,np.nan)
 save(dt,'surface_RNA_participant_screen.tsv')
 rank=dt.groupby(['cluster','markers'],as_index=False).agg(participant_records=('donor_record','size'),median_target_coverage=('target_coverage','median'),median_target_fraction_of_RNA_gate=('target_fraction_of_RNA_gate','median'),median_CNV_risk_fraction=('CNV_risk_fraction_of_RNA_gate','median'))
 rank['ranking_score']=rank.median_target_coverage*rank.median_target_fraction_of_RNA_gate
 rank=rank.sort_values(['cluster','ranking_score'],ascending=[True,False]);save(rank,'surface_candidate_descriptive_ranking.tsv')
 qa={'at':now(),'reviewed_CAF_candidates':int(main.sum()),'myCAF_cells':b.n_obs,'myCAF_subclusters':b.obs.cluster.nunique(),'myCAF_sample_records':b.obs['sample'].nunique(),'myCAF_donor_records_not_crossstudy_independent_N':b.obs.donor_record.nunique(),'patient_program_tests':len(programs),'clinical_program_rows':len(comparisons),'paired_patient_scope_rows':len(paired_rows),'surface_candidates_screened':len(surfaces),'final_metastasis_specific_protein_panel_validated':False,'all_prior_results_preserved':True}
 writej(O/'delivery/FOLLOWUP_QA.json',qa);(O/'ANALYSIS_COMPLETE').write_text(now()+'\n');status('MYCAF_RECLUSTERING_PATIENT_COMPARISONS_SURFACE_SCREEN_COMPLETE',**qa);rc=0
except Exception as e:status('FAILED',error=str(e));traceback.print_exc()
finally:(O/'exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)
