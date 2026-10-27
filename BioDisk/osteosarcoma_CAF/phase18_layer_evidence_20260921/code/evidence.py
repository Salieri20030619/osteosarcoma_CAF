from pathlib import Path
import json,hashlib,itertools
import numpy as np,pandas as pd
B=Path('C:/codex_test');I=B/'CAF_expanded_20260921';O=B/'CAF_evidence_20260921'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
layers=['strict','supported','expanded'];programs=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer']
protocol={'layers':layers,'scope':'all existing independently fitted clusters: strict10, supported18, expanded24; no new clustering parameter search',
 'evidence_levels':{'L0':'not evaluable','L1':'expression/descriptive evidence only','L2':'exploratory recurrence: >=3 known donor records in >=2 datasets, >=75% positive within-donor effects and median effect >=0.1 logCP10k','L3':'L2 plus agreement in >=3 shared known donor records with another layer: >=75% both positive, median effects>=0.1 in both layers','L4':'independent external or experimental validation; unavailable this round'},
 'working_rubric':'Exploratory bookkeeping, not a calibrated confidence score or statistical significance. Levels record positive program enrichment only; negative/contradictory evidence has separate status.',
 'replicates':'unknown patient IDs displayed per sample and excluded from donor recurrence; known repeated samples pooled within source patient; cross-source independent N unresolved',
 'program_effect':'within known donor cluster minus other cells in same analysis layer; >=10 cells on each side; >=80% gene coverage; source full-gene-space logCP10k means',
 'cross_layer':'fixed-label comparison counts only if shared population is <95% of anchor layer; identical-cell copies never upgrade evidence. Alternatively separately refit best cluster must have Jaccard>=0.3 on shared cells and >=3 shared known donors with >=75% both positive and median effect>=0.1 in each layer. These are internal same-source sensitivity checks, not independent validation.',
 'clinical':'ScPCA known same-patient bone/lung pairs descriptive; GSE162454 any-metastasis Yes2/No4 descriptive; unknown clinical states not imputed; zero cluster fraction only within an evaluable layer population; no inference or FDR',
 'identity':'CNV conflict/no-low, reference failure, ambient absent, mural/osteo flags recorded as proportions, not collapsed into program evidence',
 'limitations':'program scores used to characterize groups can overlap clustering genes; same cells/sources recur between layers, so no external validation; source-specific donor records are not verified independent pooled patients'}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
cover=read(I/'tables/panel_coverage.tsv').set_index(['sample','panel']).coverage
audit=read(B/'CAF_client_steps_20260920/tables/01_clinical_sample_audit_updated.tsv').set_index('sample')
qc=read(I/'tables/sample_layer_counts.tsv').set_index('sample')
D={};fix=[]
for layer in layers:
 d=read(I/'tables'/f'{layer}_cell_labels.tsv.gz');d['cluster']=d.cluster.astype(str);d['cell_id']=d.cell_id.astype(str)
 d['known_patient']=~d.patient_id.fillna('').astype(str).str.lower().isin(['','nan','na','unknown','unconfirmed','not_established'])
 d['analysis_donor']=np.where(d.known_patient,d.dataset+'|'+d.patient_id,d.dataset+'|'+d['sample']+'|UNKNOWN_PATIENT')
 for sm,g in d[~d.known_patient].groupby('sample'):fix.append({'layer':layer,'sample':sm,'previous_donor_record':g.donor_record.iloc[0],'analysis_record':g.analysis_donor.iloc[0],'cells':len(g),'counts_as_biological_replicate':False})
 for pr in programs:
  cv=d['sample'].map({sm:float(v) for (sm,p),v in cover.items() if p==pr});d.loc[cv<.8,pr+'_score']=np.nan
 d['site_class']=np.where(d.site.str.contains('lung',case=False,na=False),'lung',np.where(d.site.str.contains('bone|femur|tibia|humerus|fibula|pelvis|ilium|sacrum|radius|ulna',case=False,na=False),'bone','other_or_unknown'))
 d.to_csv(O/'tables'/f'{layer}_analysis_labels.tsv.gz',sep='\t',index=False);D[layer]=d
save(pd.DataFrame(fix),'unknown_patient_correction.tsv')
rows=[];summary=[];cross=[];overlap=[];counts=[];profiles=[];coupling=[]
for layer,d in D.items():
 for cl,g in d.groupby('cluster',sort=False):
  known=g[g.known_patient];n=known.groupby('analysis_donor').size()
  smry={'layer':layer,'cluster':cl,'cells':len(g),'sample_records':g['sample'].nunique(),'known_donor_records':len(n),'known_donor_records_ge10':int((n>=10).sum()),'datasets':g.dataset.nunique(),'unknown_patient_cells':int((~g.known_patient).sum()),'max_known_donor_fraction':float(n.max()/len(g)) if len(n) else np.nan,'strict_cells':int(g.strict.sum()),'supported_cells':int(g.supported.sum()),'CNV_conflict_fraction':g.CNV_state.eq('low_high_conflict').mean(),'CNV_no_low_fraction':g.CNV_state.eq('no_low_evidence').mean(),'reference_not_passed_fraction':(~g.reference_QC_pass).mean(),'ambient_missing_fraction':(~g.ambient_checked).mean(),'mural_fraction':g.mural_flag.mean(),'osteo_fraction':g.osteo_flag.mean()}
  for pr in programs:smry[pr+'_positive_cells']=int(g[pr+'_positive'].sum())
  summary.append(smry)
  for donor,gg in g.groupby('analysis_donor'):
   for pr in programs:profiles.append({'layer':layer,'cluster':cl,'donor_record':donor,'known_patient':bool(gg.known_patient.iloc[0]),'cells':len(gg),'program':pr,'mean_score':gg[pr+'_score'].mean()})
 for donor,g in d[d.known_patient].groupby('analysis_donor'):
  for cl in d.cluster.unique():
   inside=g[g.cluster==cl];outside=g[g.cluster!=cl]
   for pr in programs:
    x=inside[pr+'_score'].dropna();y=outside[pr+'_score'].dropna();ok=len(x)>=10 and len(y)>=10
    rows.append({'layer':layer,'cluster':cl,'donor_record':donor,'dataset':g.dataset.iloc[0],'program':pr,'inside_cells':len(x),'outside_cells':len(y),'eligible':ok,'delta':x.mean()-y.mean() if ok else np.nan})
 for other,od in D.items():
  if other==layer:continue
  shared=d[d.cell_id.isin(od.cell_id)]
  matched=shared[['cell_id','cluster']].merge(od[['cell_id','cluster']],on='cell_id',validate='one_to_one',suffixes=('_anchor','_other'))
  for cl in d.cluster.unique():
   q=matched[matched.cluster_anchor==cl];candidates=[]
   for oc,gg in q.groupby('cluster_other'):
    union=int((matched.cluster_anchor.eq(cl)|matched.cluster_other.eq(oc)).sum())
    candidates.append((len(gg)/union,oc,len(gg),len(gg)/len(q)))
   best=max(candidates) if candidates else (np.nan,'none',0,np.nan)
   overlap.append({'anchor_layer':layer,'anchor_cluster':cl,'other_layer':other,'anchor_total_cells':int((d.cluster==cl).sum()),'shared_anchor_cells':len(q),'best_other_cluster':best[1],'Jaccard_on_shared_population':best[0],'intersection_cells':best[2],'anchor_shared_purity':best[3]})
  for donor,g in shared[shared.known_patient].groupby('analysis_donor'):
   for cl in d.cluster.unique():
    ins=g[g.cluster==cl];outs=g[g.cluster!=cl]
    for pr in programs:
     x=ins[pr+'_score'].dropna();y=outs[pr+'_score'].dropna()
     if len(x)>=10 and len(y)>=10:cross.append({'layer':layer,'cluster':cl,'other_layer':other,'donor_record':donor,'program':pr,'delta_other':x.mean()-y.mean(),'inside_cells_other':len(x),'outside_cells_other':len(y)})
 # Sample/cluster grid includes zero clusters, even where the entire layer contains zero cells.
 for sm,qcr in qc.iterrows():
  g=d[d['sample']==sm];total=len(g)
  for cl in d.cluster.unique():
   gg=g[g.cluster==cl]
   counts.append({'layer':layer,'cluster':cl,'sample':sm,'QC_cells':int(qcr.QC_cells),'layer_cells':total,'cluster_cells':len(gg),'cluster_fraction_of_QC':len(gg)/qcr.QC_cells,'cluster_fraction_of_layer':len(gg)/total if total else np.nan,'myCAF_cells':int(gg.myCAF_customer_positive.sum()),'layer_zero_is_identity_acceptance_not_absence':total==0})
 # Descriptive association between panel scores across donor records, never cells as replicate.
 eligible=pd.DataFrame([r for r in profiles if r['layer']==layer and r['known_patient'] and r['cells']>=10])
 for cl,g in eligible.groupby('cluster'):
  pp=g.pivot(index='donor_record',columns='program',values='mean_score')
  for pr in ['Osteogenic_customer','Mineral_customer']:
   pair=pp[['myCAF_customer',pr]].dropna();r=pair.corr(method='spearman').iloc[0,1] if len(pair)>=3 and (pair.nunique()>1).all() else np.nan
   coupling.append({'layer':layer,'cluster':cl,'program_pair':'myCAF_customer vs '+pr,'donor_records':len(pair),'Spearman_descriptive':r,'P_value':'not_computed','interpretation':'across donor records; shared gene content and selection can induce association'})
E=pd.DataFrame(rows);C=pd.DataFrame(cross);S=pd.DataFrame(summary)
save(E,'within_donor_program_effects.tsv');save(C,'cross_layer_fixed_cluster_effects.tsv');save(S,'cluster_identity_and_support.tsv');save(pd.DataFrame(overlap),'independent_refit_correspondence.tsv');save(pd.DataFrame(counts),'sample_cluster_abundance.tsv');save(pd.DataFrame(profiles),'donor_cluster_programs.tsv');save(pd.DataFrame(coupling),'donor_program_coupling.tsv')
ledger=[];cons=[];refit_cons=[]
cellsets={k:set(v.cell_id) for k,v in D.items()}
changed={(k,j):len(cellsets[k]&cellsets[j])<.95*len(cellsets[k]) for k in layers for j in layers}
for (layer,cl,pr),g in E.groupby(['layer','cluster','program']):
 ev=g[g.eligible];n=len(ev);pos=int((ev.delta>0).sum());neg=int((ev.delta<0).sum());median=ev.delta.median();frac=pos/n if n else np.nan
 grade='L0' if n==0 else 'L1';state='not_evaluable' if n==0 else ('positive_direction_mixed' if median>0 else 'current_data_do_not_support_enrichment')
 recurrent=n>=3 and ev.dataset.nunique()>=2 and frac>=.75 and median>=.1
 if recurrent:grade='L2';state='exploratory_recurrent_enrichment'
 best_other='';best_n=0;best_fraction=np.nan
 for other in layers:
  if other==layer:continue
  ce=C[(C.layer==layer)&(C.cluster==cl)&(C.program==pr)&(C.other_layer==other)]
  v=ev.merge(ce,on=['layer','cluster','program','donor_record'],how='inner',validate='one_to_one')
  nn=len(v);both=int(((v.delta>0)&(v.delta_other>0)).sum());same=int((np.sign(v.delta)==np.sign(v.delta_other)).sum())
  changed_population=changed[(layer,other)]
  passed=changed_population and nn>=3 and both/nn>=.75 and v.delta.median()>=.1 and v.delta_other.median()>=.1
  cons.append({'layer':layer,'cluster':cl,'program':pr,'other_layer':other,'paired_donor_records':nn,'both_positive':both,'same_direction':same,'anchor_median':v.delta.median(),'other_median':v.delta_other.median(),'population_changed':changed_population,'cross_layer_positive_pass':passed})
  if recurrent and passed:grade='L3';state='cross_layer_exploratory_support';best_other=other;best_n=nn;best_fraction=both/nn
  match=next(z for z in overlap if z['anchor_layer']==layer and z['anchor_cluster']==cl and z['other_layer']==other)
  oe=E[(E.layer==other)&(E.cluster==match['best_other_cluster'])&(E.program==pr)&E.eligible]
  rv=ev[['donor_record','delta']].merge(oe[['donor_record','delta']],on='donor_record',suffixes=('_anchor','_other'),validate='one_to_one')
  rn=len(rv);rp=int(((rv.delta_anchor>0)&(rv.delta_other>0)).sum())
  refitpass=match['Jaccard_on_shared_population']>=.3 and rn>=3 and rp/rn>=.75 and rv.delta_anchor.median()>=.1 and rv.delta_other.median()>=.1
  refit_cons.append({'layer':layer,'cluster':cl,'program':pr,'other_layer':other,'matched_other_cluster':match['best_other_cluster'],'Jaccard_on_shared_population':match['Jaccard_on_shared_population'],'paired_donor_records':rn,'both_positive':rp,'anchor_median':rv.delta_anchor.median(),'other_median':rv.delta_other.median(),'refit_positive_pass':refitpass})
  if recurrent and refitpass:grade='L3';state='cross_layer_exploratory_support';best_other=other+'_separately_refit';best_n=rn;best_fraction=rp/rn
 r=S[(S.layer==layer)&(S.cluster==cl)].iloc[0]
 ledger.append({'layer':layer,'cluster':cl,'program':pr,'evidence_level':grade,'evidence_status':state,'evaluable_known_donor_records':n,'datasets':ev.dataset.nunique(),'positive_records':pos,'negative_records':neg,'positive_fraction':frac,'median_within_donor_delta':median,'supporting_other_layer':best_other,'cross_layer_paired_records':best_n,'cross_layer_both_positive_fraction':best_fraction,'mural_fraction':r.mural_fraction,'CNV_conflict_fraction':r.CNV_conflict_fraction,'CNV_no_low_fraction':r.CNV_no_low_fraction,'identity_interpretation':'operational candidate; evaluate risk columns separately','clinical_evidence':'see clinical tables; not implied by program level','external_validation':'none','evidence_source':'within_donor_program_effects.tsv; cross_layer_consistency.tsv'})
L=pd.DataFrame(ledger);save(L,'program_evidence_ledger.tsv');save(pd.DataFrame(cons),'cross_layer_consistency.tsv');save(pd.DataFrame(refit_cons),'separately_refit_program_consistency.tsv')
# Clinical contrasts are evaluated separately from program recurrence.
pairs=[];clinical_group=[];clinical_ledger=[]
known_pairs={'SCPCP000017|CCG1106_006','SCPCP000017|CCG1106_009','SCPCP000023|SJOS013768','SCPCP000023|SJOS063833'}
for layer,d in D.items():
 for donor in sorted(known_pairs):
  g=d[d.analysis_donor==donor]
  for cl in d.cluster.unique():
   for subset in ['all','myCAF']:
    gg=g if subset=='all' else g[g.myCAF_customer_positive]
    b=gg[gg.site_class=='bone'];l=gg[gg.site_class=='lung'];bc=b[b.cluster==cl];lc=l[l.cluster==cl]
    pairs.append({'layer':layer,'cluster':cl,'donor_record':donor,'set':subset,'bone_layer_cells':len(b),'lung_layer_cells':len(l),'bone_cluster_cells':len(bc),'lung_cluster_cells':len(lc),'abundance_eligible':len(b)>=10 and len(l)>=10,'lung_minus_bone_fraction':len(lc)/len(l)-len(bc)/len(b) if len(l)>=10 and len(b)>=10 else np.nan,'interpretation':'same-patient site contrast; treatment and time may differ; not metastasis causation'})
    for pr in programs:
     x=lc[pr+'_score'].dropna();y=bc[pr+'_score'].dropna()
     clinical_group.append({'contrast':'ScPCA_paired_site','layer':layer,'cluster':cl,'donor_or_group':donor,'set':subset,'program':pr,'n_lung_cells':len(x),'n_bone_cells':len(y),'evaluable':len(x)>=10 and len(y)>=10,'effect':x.mean()-y.mean() if len(x)>=10 and len(y)>=10 else np.nan,'direction':'lung_minus_bone'})
 # GSE162454: any metastasis status, not verified lung-metastasis group.
 met=read(B/'CAF_client_steps_20260920/tables/01_GSE162454_cross_publication_clinical.tsv').set_index('sample').metastasis_any
 gd=d[d.dataset=='GSE162454'].copy();gd['metastasis_any']=gd['sample'].map(met)
 for cl in d.cluster.unique():
  for pr in programs:
   records=[]
   for sm,g in gd[gd.cluster==cl].groupby('sample'):
    v=g[pr+'_score'].dropna()
    if len(v)>=10:records.append({'sample':sm,'group':met.loc[sm],'value':v.mean()})
   rr=pd.DataFrame(records,columns=['sample','group','value']);yes=rr[rr.group=='Yes'].value;no=rr[rr.group=='No'].value
   clinical_group.append({'contrast':'GSE162454_any_metastasis','layer':layer,'cluster':cl,'donor_or_group':'Yes_vs_No','set':'all','program':pr,'n_yes_donors':len(yes),'n_no_donors':len(no),'evaluable':len(yes)>=1 and len(no)>=1,'effect':yes.median()-no.median() if len(yes) and len(no) else np.nan,'direction':'median_Yes_minus_No','interpretation':'any-metastasis label; organ/timing unresolved; maximum2 vs4, descriptive only'})
 for sm,g in gd.groupby('sample'):
  for cl in d.cluster.unique():
   clinical_ledger.append({'layer':layer,'sample':sm,'donor_record':g.analysis_donor.iloc[0],'metastasis_any':met.loc[sm],'cluster':cl,'layer_cells':len(g),'cluster_cells':int((g.cluster==cl).sum()),'cluster_fraction':float((g.cluster==cl).mean())})
save(pd.DataFrame(pairs),'paired_site_cluster_abundance.tsv');save(pd.DataFrame(clinical_group),'clinical_program_effects.tsv');save(pd.DataFrame(clinical_ledger),'GSE162454_any_metastasis_sample_abundance.tsv')
pairdf=pd.DataFrame(pairs);ce=pd.DataFrame(clinical_group);cards=[]
for r in S.to_dict('records'):
 layer=r['layer'];cl=r['cluster'];le=L[(L.layer==layer)&(L.cluster==cl)];p=pairdf[(pairdf.layer==layer)&(pairdf.cluster==cl)&(pairdf['set']=='myCAF')&pairdf.abundance_eligible]
 r['highest_program_level']=max(le.evidence_level) if len(le) else 'L0';r['L2_L3_programs']=';'.join(le[le.evidence_level.isin(['L2','L3'])].program)
 r['paired_site_evaluable_records']=len(p);r['paired_site_lung_enriched_records']=int((p.lung_minus_bone_fraction>0).sum());r['paired_site_lung_depleted_records']=int((p.lung_minus_bone_fraction<0).sum());r['clinical_support_level']='C1_descriptive_site_or_any_metastasis_only' if len(p) or ce[(ce.layer==layer)&(ce.cluster==cl)&ce.evaluable].shape[0] else 'C0_not_evaluable'
 r['identity_risk_tags']=';'.join(name for name,val in [('mural_enriched',r['mural_fraction']>=.2),('CNV_conflict_enriched',r['CNV_conflict_fraction']>=.2),('CNV_unassessed_enriched',r['CNV_no_low_fraction']>=.2),('reference_QC_gap',r['reference_not_passed_fraction']>=.2),('ambient_gap',r['ambient_missing_fraction']>=.2)] if val)
 cards.append(r)
save(pd.DataFrame(cards),'cluster_evidence_cards.tsv')
q={'analysis_layers':3,'clusters':len(S),'program_claims':len(L),'unknown_patient_policy_corrected':True,'known_cross_cohort_independent_N':None,'inference_tests':False,'external_validation':False,'source_labels_unchanged':True,'level_counts':L.evidence_level.value_counts().to_dict()}
assert len(S)==52 and len(L)==260
(O/'delivery/QA_local.json').write_text(json.dumps(q,indent=2));print(json.dumps(q),flush=True)
