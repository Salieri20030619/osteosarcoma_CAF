from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,itertools,math,gc
import pandas as pd,numpy as np,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase31_customer_completion_20260926';M=P/'phase30_two_layer_classification_20260926';D=O/'clinical_internal';D.mkdir(exist_ok=True)
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
def bh(p):
 p=np.asarray(p,float);valid=np.isfinite(p);v=np.where(valid,p,1);ix=np.argsort(v);q=np.empty(len(v));q[ix]=np.minimum(1,np.minimum.accumulate((v[ix]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);q[~valid]=np.nan;return q
def test(a,b,paired):
 a=np.atleast_2d(a);b=np.atleast_2d(b);ng=a.shape[1];n1=len(a);n2=len(b)
 effect=b.mean(0)-a.mean(0);pv=np.full(ng,np.nan);positive=np.full(ng,np.nan);negative=positive.copy();consistency=positive.copy()
 if paired:
  delta=b-a;positive=(delta>0).sum(0);negative=(delta<0).sum(0);consistency=positive/n1;median=np.median(delta,0)
  if n1<2:return effect,median,pv,positive,negative,consistency,'no_valid_independent_replication'
  exact=n1<=15;rng=np.random.default_rng(20260926)
  weights=np.array(list(itertools.product([-1.,1.],repeat=n1)),np.float32)/n1 if exact else rng.choice([-1.,1.],size=(65536,n1)).astype(np.float32)/n1
  vals=delta
 else:
  median=np.median(b,0)-np.median(a,0)
  if n1<2 or n2<2:return effect,median,pv,positive,negative,consistency,'no_valid_independent_replication'
  vals=np.vstack([a,b]);nt=n1+n2;total=math.comb(nt,n2);exact=total<=65536
  choices=list(itertools.combinations(range(nt),n2)) if exact else [np.random.default_rng(20260926+i).choice(nt,n2,replace=False) for i in range(65536)]
  weights=np.full((len(choices),nt),-1/n1,np.float32)
  for i,ix in enumerate(choices):weights[i,list(ix)]=1/n2
 counts=np.zeros(ng,int);observed=np.abs(effect)
 for j in range(0,ng,512):
  simulated=np.abs(weights@vals[:,j:j+512]);counts[j:j+512]=(simulated>=observed[j:j+512]-1e-7).sum(0)
 pv=counts/len(weights) if exact else (counts+1)/(len(weights)+1)
 return effect,median,pv,positive,negative,consistency,('exact_' if exact else 'MC65536_')+('patient_signflip' if paired else 'patient_permutation')
clinical=read(O/'resources/clinical_mapping_frozen.tsv').set_index('sample')
contrasts=[
 ('clinical_group','A_reported_nonmetastasis','B1_history_supported_timing_unresolved','A_reported__B1_history'),
 ('clinical_group','A_reported_nonmetastasis','B2','A_reported__B2'),
 ('clinical_group','A_candidate_reported_no_metastasis_timing_unknown','B1_history_supported_timing_unresolved','A_timing_unknown__B1_history'),
 ('clinical_group','A_candidate_reported_no_metastasis_timing_unknown','B2','A_timing_unknown__B2'),
 ('clinical_group','B1_history_supported_timing_unresolved','B2','B1_history__B2'),
 ('anatomic_site','Bone','Lung','bone__lung_site'),
 ('treatment','No','Yes','reported_treatment_No__Yes'),
 ('treatment','Surgery','Chemotherapy, surgery','Surgery__Chemotherapy_surgery')]
writej(D/'CONTRASTS_FROZEN.json',{'at':now(),'contrasts':contrasts,'direction':'second minus first','no_unknown_reassignment':True,'paired_identity':'same source study and donor; contemporaneous sampling not assumed','replication':'P unavailable below two independent donors per unpaired group or below two paired donors','no_cell_count_entry_threshold':True,'permutation_precision':'exact paired up to 15 donors; exact group allocation if <=65536 assignments; otherwise fixed 65536 Monte Carlo, seeded before results','gene_space':'common measured 15391 genes, raw count pseudobulk normalized by original full-library counts; no expression significance prefilter','FDR_family':'layer, analysis, source study, contrast, paired/unpaired, feature family; includes all population targets within each family'})
allres=[];avail=[];roster=[];pairedcounts=[]
try:
 for layer,expected in [('supported',1224),('expanded',8174)]:
  a=ad.read_h5ad(M/f'objects/{layer}_integrated.h5ad');obs=read(M/f'tables/{layer}_classified_cells.tsv.gz').set_index('cell_id').loc[a.obs_names].copy()
  assert len(obs)==expected and obs.index.is_unique and 'counts' in a.layers
  for col in ['donor_record','dataset','site','clinical_group','treatment','primary_or_metastasis']:
   obs[col]=obs['sample'].map(clinical[col])
  assert obs[['donor_record','dataset','site']].notna().all().all()
  obs['anatomic_site']=obs.site.str.lower().map(lambda s:'Lung' if 'lung' in s else ('Bone' if any(k in s for k in ['bone','femur','tibia','humerus','calcaneus','ulna']) else 'Other_or_unspecified'))
  my=read(M/f'tables/{layer}_myCAF_cell_labels.tsv.gz').set_index('cell_id')
  obs['my_subgroup']=my.cluster.reindex(obs.index).map(lambda x:str(int(x)) if pd.notna(x) else 'not_myCAF')
  prog=[c for c in obs if c.startswith('program_') and not c.endswith('descriptive')]
  x=a.layers['counts'].tocsr();den=obs.full_library_counts.to_numpy(float);assert np.all(den>0)
  geneids=a.var_names.astype(str).tolist();symbols=a.var['symbol'].astype(str).tolist() if 'symbol' in a.var else geneids
  names=geneids+prog;symbols=symbols+prog;families=['genes']*len(geneids)+['programs']*len(prog)
  def aggregate(idx):
   idx=np.asarray(idx,int);return np.r_[np.log2(np.asarray(x[idx].sum(0)).ravel()/den[idx].sum()*1e4+.1),obs.iloc[idx][prog].mean(0).to_numpy()]
  def record(v1,v2,meta,paired,ids1,ids2,ncells1,ncells2):
   avail.append(dict(**meta,paired=paired,n_first=len(v1),n_second=len(v2),cells_first=ncells1,cells_second=ncells2,status='evaluated' if len(v1) and len(v2) else 'no_cells_in_one_or_both_conditions'))
   if not len(v1) or not len(v2):return
   effect,median,pv,pos,neg,cons,method=test(np.stack(v1),np.stack(v2),paired)
   frame=pd.DataFrame(dict(feature=names,symbol=symbols,feature_family=families,mean_effect=effect,median_effect=median,P=pv,positive_donors=pos,negative_donors=neg,positive_fraction=cons))
   for k,v in dict(**meta,paired=paired,n_first=len(v1),n_second=len(v2),test_method=method).items():frame[k]=v
   allres.append(frame)
   for condition,ids in [('first',ids1),('second',ids2)]:
    for donor in ids:roster.append(dict(**meta,paired=paired,condition=condition,donor_record=donor))
  # myCAF subgroups versus other myCAF, matched within exact raw anatomical site then averaged per donor.
  for study in sorted(obs.dataset.unique()):
   for target in sorted(my.cluster.astype(str).unique(),key=int):
    ds={};cells1=cells2=0
    use=obs.dataset.eq(study)&obs.my_subgroup.ne('not_myCAF')
    for (donor,site),sub in obs.loc[use].groupby(['donor_record','site']):
     i1=obs.index.get_indexer(sub.index[sub.my_subgroup.ne(target)]);i2=obs.index.get_indexer(sub.index[sub.my_subgroup.eq(target)])
     if len(i1) and len(i2):
      ds.setdefault(donor,[]).append((aggregate(i1),aggregate(i2)));cells1+=len(i1);cells2+=len(i2)
    ids=sorted(ds);v1=[np.mean([r[0] for r in ds[z]],0) for z in ids];v2=[np.mean([r[1] for r in ds[z]],0) for z in ids]
    record(v1,v2,dict(layer=layer,analysis='myCAF_internal',study=study,target='subgroup_'+target,contrast='subgroup_vs_other_myCAF'),True,ids,ids,cells1,cells2)
  # Whole candidate layer, main classes and every actual myCAF subgroup, with no minimum cell gate.
  targets={'all_CAF_candidates':np.ones(len(obs),bool)}
  targets.update({c:obs.class_id.eq(c).to_numpy() for c in sorted(obs.class_id.unique())})
  targets.update({'my_subgroup_'+str(c):obs.my_subgroup.eq(str(c)).to_numpy() for c in sorted(my.cluster.unique())})
  for study in sorted(obs.dataset.unique()):
   study_mask=obs.dataset.eq(study).to_numpy()
   for field,level1,level2,cname in contrasts:
    for target,tmask in targets.items():
     selections=[]
     for lev in [level1,level2]:
      sub=obs.loc[study_mask&tmask&obs[field].eq(lev).to_numpy()]
      selections.append({d:obs.index.get_indexer(z.index) for d,z in sub.groupby('donor_record')})
     s1,s2=selections;shared=sorted(set(s1)&set(s2))
     meta=dict(layer=layer,analysis='clinical_expression',study=study,target=target,contrast=cname)
     record([aggregate(s1[z]) for z in shared],[aggregate(s2[z]) for z in shared],meta,True,shared,shared,sum(len(s1[z]) for z in shared),sum(len(s2[z]) for z in shared))
     ids1=sorted(set(s1)-set(shared));ids2=sorted(set(s2)-set(shared))
     record([aggregate(s1[z]) for z in ids1],[aggregate(s2[z]) for z in ids2],meta,False,ids1,ids2,sum(len(s1[z]) for z in ids1),sum(len(s2[z]) for z in ids2))
  # Descriptive site proportions include zero target cells, while zero denominator remains undefined.
  pairdefs=read(O/'resources/bone_lung_pairs_frozen.tsv')
  for r in pairdefs.to_dict('records'):
   for site,field in [('Bone','bone_samples'),('Lung','lung_samples')]:
    sub=obs.loc[obs['sample'].isin(str(r[field]).split(';'))]
    for target,tmask in targets.items():
     count=int((tmask & obs.index.isin(sub.index)).sum())
     pairedcounts.append(dict(layer=layer,donor_record=r['donor_record'],site=site,target=target,target_cells=count,layer_cells=len(sub),fraction=count/len(sub) if len(sub) else np.nan,pair_type=r['pair_type']))
  writej(O/'CLINICAL_STATUS.json',{'status':'RUNNING','finished_layer':layer,'at':now(),'result_blocks':len(allres)})
  del a,x;gc.collect()
 res=pd.concat(allres,ignore_index=True) if allres else pd.DataFrame()
 if len(res):
  family=['layer','analysis','study','contrast','paired','feature_family']
  res['FDR']=res.groupby(family,dropna=False,sort=False).P.transform(bh)
  res['nominal_P_lt_0_05']=res.P.lt(.05);res['exploratory_small_N']=np.where(res.paired,res.n_first,(res.n_first+res.n_second))<10
  for (layer,analysis),frame in res.groupby(['layer','analysis'],sort=False):
   frame.to_csv(D/f'{layer}_{analysis}_all_feature_statistics.tsv.gz',sep='\t',index=False)
  res.loc[res.feature_family.eq('programs')].to_csv(D/'all_program_statistics.tsv',sep='\t',index=False)
 pd.DataFrame(avail).to_csv(D/'comparison_availability.tsv',sep='\t',index=False)
 pd.DataFrame(roster).to_csv(D/'independent_donor_roster.tsv.gz',sep='\t',index=False)
 pd.DataFrame(pairedcounts).to_csv(D/'paired_bone_lung_counts_and_proportions.tsv',sep='\t',index=False)
 qa={'at':now(),'result_rows':len(res),'availability_records':len(avail),'supported_cells':1224,'expanded_cells':8174,'no_primary_cells_deleted':True,'clinical_labels_from_frozen_source_mapping':True,'no_minimum_cell_entry_gate':True,'all_P_in_0_1':bool(res.P.dropna().between(0,1).all()),'all_FDR_in_0_1':bool(res.FDR.dropna().between(0,1).all()),'paired_myCAF_available_donors_must_be_read_from_roster':True,'clinical_proportion_outputs_descriptive_no_invented_metastasis_conclusion':True}
 assert qa['all_P_in_0_1'] and qa['all_FDR_in_0_1']
 writej(D/'QA.json',qa);(D/'COMPLETE').write_text(now());(D/'exit_code').write_text('0\n');writej(O/'CLINICAL_STATUS.json',dict(status='COMPUTED_REVIEW_REQUIRED',**qa))
except Exception as ex:
 traceback.print_exc();writej(O/'CLINICAL_STATUS.json',{'status':'FAILED','at':now(),'error':str(ex)});(D/'exit_code').write_text('1\n');raise
