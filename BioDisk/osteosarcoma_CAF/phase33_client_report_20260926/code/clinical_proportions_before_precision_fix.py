from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,itertools,math,traceback
import numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase33_client_report_20260926';M=P/'phase30_two_layer_classification_20260926';Q=P/'phase31_customer_completion_20260926'
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
clinical=read(Q/'resources/clinical_mapping_frozen.tsv').set_index('sample')
clinical['anatomic_site']=clinical.site.str.lower().map(lambda s:'Lung' if 'lung' in s else ('Bone' if any(k in s for k in ['bone','femur','tibia','humerus','calcaneus','ulna']) else 'Other_or_unspecified'))
full={sm:json.loads((Q/'decontx'/sm/'parameters.json').read_text())['input_cells'] for sm in clinical.index}
contrasts=json.loads((Q/'clinical_internal/CONTRASTS_FROZEN.json').read_text())['contrasts']
contrasts.append(['clinical_group','A_candidate_reported_no_metastasis_timing_unknown','metastatic_bone_primary_lung_and_timing_unconfirmed','A_timing_unknown__any_metastasis_bone_unconfirmed_lung'])
writej(O/'resources/PROPORTION_PROTOCOL.json',{'at':now(),'contrasts':contrasts,'all107_records_in_capture_denominator':True,'denominators':['all full called cells','within frozen candidate layer','within frozen F2 candidate class'],'zero_numerator_retained':True,'zero_denominator_not_interpreted_as_zero_proportion':True,'independent_unit':'source-study donor; repeated records aggregated per condition before testing; shared donors only paired and excluded from unpaired subset','no_minimum_cell_gate':True,'P_screen':'original P<0.05, both directions; report FDR in separate same-family tests','FDR_families':'layer, study, contrast, paired, denominator family; all targets within family','precision':'exact sign-flip <=15 pairs; exact unpaired allocations <=65536 otherwise seeded65536 Monte Carlo','capture_proportion_not_tissue_absolute_abundance':True})
cache={}
def weights(n1,n2,paired):
 key=(n1,n2,paired)
 if key in cache:return cache[key]
 if paired:
  exact=n1<=15
  if exact:w=np.array(list(itertools.product([-1.,1.],repeat=n1)),np.float32)/n1
  else:w=np.random.default_rng(20260926).choice([-1.,1.],size=(65536,n1)).astype(np.float32)/n1
 else:
  n=n1+n2;exact=math.comb(n,n2)<=65536
  if exact:ix=np.array(list(itertools.combinations(range(n),n2)))
  else:ix=np.argpartition(np.random.default_rng(20260926).random((65536,n)),n2-1,axis=1)[:,:n2]
  w=np.full((len(ix),n),-1/n1,np.float32);w[np.arange(len(ix))[:,None],ix]=1/n2
 cache[key]=(w,exact);return w,exact
def stat(a,b,paired):
 n1,n2=len(a),len(b)
 if not n1 or not n2:return dict(mean_effect=np.nan,median_effect=np.nan,P=np.nan,test='no_observable_condition')
 effect=float(b.mean()-a.mean());median=float(np.median(b-a) if paired else np.median(b)-np.median(a))
 out=dict(mean_effect=effect,median_effect=median,mean_first=float(a.mean()),mean_second=float(b.mean()),positive_pairs=int((b-a>0).sum()) if paired else None)
 if (paired and n1<2) or (not paired and min(n1,n2)<2):return dict(**out,P=np.nan,test='no_valid_independent_replication')
 val=b-a if paired else np.r_[a,b]
 if np.all(val==0) or (not paired and np.all(val==val[0])):return dict(**out,P=np.nan,test='all_zero_or_invariant_not_tested')
 w,exact=weights(n1,n2,paired);null=np.abs(w@val);k=int((null>=abs(effect)-1e-7).sum())
 return dict(**out,P=k/len(w) if exact else (k+1)/(len(w)+1),test=('exact' if exact else 'MC65536')+('_paired_signflip' if paired else '_patient_permutation'))
def bh(x):
 p=x.to_numpy(float);valid=np.isfinite(p);v=np.where(valid,p,1);ix=np.argsort(v);q=np.empty(len(v));q[ix]=np.minimum(1,np.minimum.accumulate((v[ix]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);q[~valid]=np.nan;return q
rows=[];samples=[];roster=[]
try:
 for layer,expected in [('supported',1224),('expanded',8174)]:
  labels=read(M/f'tables/{layer}_classified_cells.tsv.gz');assert len(labels)==expected
  my=read(M/f'tables/{layer}_myCAF_cell_labels.tsv.gz')
  total=labels.groupby('sample').size().reindex(clinical.index,fill_value=0)
  f2=labels.loc[labels.class_id.eq('F2')].groupby('sample').size().reindex(clinical.index,fill_value=0)
  features=[('all_candidates','full_called_capture',total,pd.Series(full))]
  for cls in ['F1','F2','F3','F4','F5']:
   num=labels.loc[labels.class_id.eq(cls)].groupby('sample').size().reindex(clinical.index,fill_value=0)
   features += [(cls,'full_called_capture',num,pd.Series(full)),(cls,'within_candidate_layer',num,total)]
  for cl in sorted(my.cluster.unique()):
   num=my.loc[my.cluster.eq(cl)].groupby('sample').size().reindex(clinical.index,fill_value=0)
   features.append(('my_subgroup_'+str(cl),'within_F2',num,f2))
  for target,denom,num,den in features:
   data=clinical.copy();data['numerator']=num;data['denominator']=den
   assert data[['numerator','denominator']].notna().all().all() and (data.numerator<=data.denominator).all()
   for sm,r in data.iterrows():samples.append(dict(layer=layer,sample=sm,dataset=r.dataset,donor_record=r.donor_record,clinical_group=r.clinical_group,site=r.site,treatment=r.treatment,target=target,denominator_family=denom,numerator=int(r.numerator),denominator=int(r.denominator),fraction=r.numerator/r.denominator if r.denominator else np.nan))
   for study,ds in data.groupby('dataset'):
    for field,lev1,lev2,cname in contrasts:
     group=[]
     for lev in [lev1,lev2]:
      z=ds.loc[ds[field].eq(lev)].groupby('donor_record')[['numerator','denominator']].sum()
      z=z.loc[z.denominator>0];z['fraction']=z.numerator/z.denominator;group.append(z)
     a,b=group;shared=a.index.intersection(b.index)
     for paired in [True,False]:
      i1=shared if paired else a.index.difference(shared);i2=shared if paired else b.index.difference(shared)
      av=a.loc[i1,'fraction'].to_numpy(float);bv=b.loc[i2,'fraction'].to_numpy(float)
      meta=dict(layer=layer,study=study,target=target,denominator_family=denom,contrast=cname,condition_first=lev1,condition_second=lev2,paired=paired,n_first=len(av),n_second=len(bv))
      result=stat(av,bv,paired);rows.append(dict(**meta,**result))
      for cond,ids,df in [('first',i1,a),('second',i2,b)]:
       for dn in ids:roster.append(dict(**meta,condition=cond,donor_record=dn,numerator=int(df.loc[dn,'numerator']),denominator=int(df.loc[dn,'denominator']),fraction=float(df.loc[dn,'fraction'])))
  writej(O/'PROPORTION_STATUS.json',{'status':'RUNNING','finished_layer':layer,'at':now()})
 res=pd.DataFrame(rows);fam=['layer','study','contrast','paired','denominator_family'];res['FDR']=res.groupby(fam,sort=False,dropna=False).P.transform(bh)
 res['nominal_P_lt_0_05']=res.P.lt(.05);res['family']=res[fam].astype(str).agg('|'.join,axis=1)
 res.to_csv(O/'tables/clinical_proportion_statistics.tsv',sep='\t',index=False)
 pd.DataFrame(samples).to_csv(O/'tables/clinical_proportion_sample_counts.tsv.gz',sep='\t',index=False)
 rr=pd.DataFrame(roster);rr.to_csv(O/'tables/clinical_proportion_donor_roster.tsv.gz',sep='\t',index=False)
 assert res.P.dropna().between(0,1).all() and res.FDR.dropna().between(0,1).all()
 assert not rr.duplicated(['layer','study','target','denominator_family','contrast','paired','condition','donor_record']).any()
 qa={'at':now(),'comparisons':len(res),'valid_tests':int(res.P.notna().sum()),'nominal_tests':int(res.P.lt(.05).sum()),'full_called_total':int(sum(full.values())),'source_records':len(clinical),'counts_and_unique_patient_roster_validated':True,'no_new_cell_filters':True}
 writej(O/'delivery/PROPORTION_QA.json',qa);(O/'proportions.exit_code').write_text('0');writej(O/'PROPORTION_STATUS.json',{'status':'COMPUTED_REVIEW_REQUIRED',**qa})
 print(res.loc[res.P.lt(.05),['layer','study','target','denominator_family','contrast','paired','n_first','n_second','mean_effect','P','FDR']].to_string(index=False))
except Exception as e:
 traceback.print_exc();(O/'proportions.exit_code').write_text('1');writej(O/'PROPORTION_STATUS.json',{'status':'FAILED','error':str(e),'at':now()});raise
