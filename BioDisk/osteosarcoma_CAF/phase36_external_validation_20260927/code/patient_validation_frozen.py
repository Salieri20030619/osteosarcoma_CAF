from pathlib import Path
import json,datetime,itertools,hashlib,sys,traceback
import pandas as pd,numpy as np,anndata as ad
from scipy.stats import spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase36_external_validation_20260927';E=P/'phase35_external_validation_prep_20260927'
p=json.loads((O/'resources/VALIDATION_PROTOCOL_FROZEN.json').read_text());meta=pd.read_csv(O/'resources/sample_patient_manifest_corrected.tsv',sep='\t',keep_default_na=False).set_index('sample');genes=p['gene_features'];panels=p['programs'];markers=p['markers'];rng=np.random.default_rng(20260927)
def now():return datetime.datetime.now().astimezone().isoformat()
def bh(vals):
 a=np.asarray(vals,float);out=np.full(len(a),np.nan);ok=np.where(np.isfinite(a))[0]
 if len(ok):
  ix=ok[np.argsort(a[ok])];z=a[ix]*len(ix)/np.arange(1,len(ix)+1);out[ix]=np.minimum(1,np.minimum.accumulate(z[::-1])[::-1])
 return out
def signflip(d):
 d=np.asarray(d,float);d=d[np.isfinite(d)];n=len(d)
 if n<2 or np.all(d==0):return np.nan,'not_tested_insufficient_or_zero_contrast'
 obs=abs(d.mean());tol=1e-12
 if n<=15:
  bits=((np.arange(2**n,dtype=np.uint64)[:,None]>>np.arange(n,dtype=np.uint64))&1).astype(float)*2-1
  return float(np.mean(np.abs(bits@d/n)>=obs-tol)),'exact_sign_flip'
 hit=0
 for _ in range(16):hit+=int(np.count_nonzero(np.abs(rng.choice([-1.,1.],size=(4096,n))@d/n)>=obs-tol))
 return (hit+1)/65537,'Monte_Carlo_sign_flip_65536'
def permute(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float);a=a[np.isfinite(a)];b=b[np.isfinite(b)]
 if len(a)<2 or len(b)<2 or np.all(np.r_[a,b]==np.r_[a,b][0]):return np.nan,'not_tested_insufficient_or_constant'
 z=np.r_[a,b];n=len(a);obs=abs(b.mean()-a.mean());hits=0;total=0
 for c in itertools.combinations(range(len(z)),n):
  aa=z[list(c)];bb=(z.sum()-aa.sum())/len(b);hits+=abs(bb-aa.mean())>=obs-1e-12;total+=1
 return float(hits/total),'exact_label_permutation_'+str(total)
sample=[];program=[];associations=[];coverage=[]
for idx,(sm,row) in enumerate(meta.iterrows()):
 a=ad.read_h5ad(E/'objects/GSE293065'/f'{sm}_counts_spatial.h5ad');lib=np.asarray(a.X.sum(axis=1)).ravel().astype(float);total=float(lib.sum());assert total>0
 zero=lib==0;symbols=a.var.gene_symbol.astype(str).to_numpy();raw=np.zeros((a.n_obs,len(genes)),dtype=float);present=[]
 for k,g in enumerate(genes):
  ids=np.flatnonzero(symbols==g);present.append(bool(len(ids)))
  raw[:,k]=np.asarray(a.X[:,ids].sum(axis=1)).ravel() if len(ids) else np.nan
 cp=np.divide(raw,lib[:,None],out=np.full_like(raw,np.nan),where=lib[:,None]>0)*10000;log=np.log1p(cp)
 for k,g in enumerate(genes):
  sample.append(dict(sample=sm,patient_id=row.patient_id,site_group=row.site_group,event_month=row.months_from_diagnosis,feature=g,feature_type='gene',counts=float(raw[:,k].sum()) if present[k] else np.nan,library_counts=total,log2_CP10k_plus_0_1=float(np.log2(raw[:,k].sum()/total*10000+.1)) if present[k] else np.nan,spots=a.n_obs))
 scores={}
 for name,gs in panels.items():
  ids=[genes.index(g) for g in gs if present[genes.index(g)]];scores[name]=np.nanmean(log[:,ids],axis=1) if ids else np.full(a.n_obs,np.nan)
  program.append(dict(sample=sm,patient_id=row.patient_id,site_group=row.site_group,event_month=row.months_from_diagnosis,feature=name,feature_type='program',score=float(np.nanmean(scores[name])),spots=a.n_obs))
  coverage.append(dict(sample=sm,program=name,requested=len(gs),measured=len(ids),missing=';'.join(g for g in gs if not present[genes.index(g)]),zero_library_spots=int(zero.sum())))
 for marker in markers:
  m=log[:,genes.index(marker)]
  for name,gs in panels.items():
   use=[genes.index(g) for g in gs if g!=marker and present[genes.index(g)]];s=np.nanmean(log[:,use],axis=1) if use else np.full(a.n_obs,np.nan);valid=np.isfinite(m)&np.isfinite(s)
   rho=float(spearmanr(m[valid],s[valid]).statistic) if valid.sum()>1 and np.ptp(m[valid])>0 and np.ptp(s[valid])>0 else np.nan
   associations.append(dict(sample=sm,patient_id=row.patient_id,site_group=row.site_group,event_month=row.months_from_diagnosis,marker=marker,program=name,rho=rho,evaluable_spots=int(valid.sum()),spot_P_not_used=True,tested_marker_excluded_from_program=True))
 (O/'PATIENT_VALIDATION_STATUS.json').write_text(json.dumps({'at':now(),'stage':'sample_aggregation','samples_finished':idx+1,'samples_total':26}))
sample=pd.DataFrame(sample);program=pd.DataFrame(program);assoc=pd.DataFrame(associations)
sample.to_csv(O/'tables/sample_gene_pseudobulk.tsv',sep='\t',index=False);program.to_csv(O/'tables/sample_program_scores.tsv',sep='\t',index=False);assoc.to_csv(O/'tables/sample_marker_program_correlations.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/program_feature_coverage.tsv',sep='\t',index=False)
events=sample.groupby(['patient_id','site_group','event_month','feature'],dropna=False).agg(counts=('counts',lambda x:x.sum(min_count=1)),library_counts=('library_counts','sum')).reset_index();events['value']=np.log2(events.counts/events.library_counts*10000+.1);events['feature_type']='gene'
progevents=program.groupby(['patient_id','site_group','event_month','feature'],dropna=False).score.mean().reset_index().rename(columns={'score':'value'});progevents['feature_type']='program'
all_events=pd.concat([events[['patient_id','site_group','event_month','feature','feature_type','value']],progevents],ignore_index=True)
values=all_events.groupby(['patient_id','site_group','feature','feature_type']).value.mean().reset_index();values.to_csv(O/'tables/patient_site_feature_values.tsv',sep='\t',index=False);all_events.to_csv(O/'tables/patient_event_feature_values.tsv',sep='\t',index=False)
primary=set(meta.loc[meta.site_group.eq('primary_tumor'),'patient_id']);lung=set(meta.loc[meta.site_group.eq('lung_metastasis'),'patient_id']);shared=sorted(primary&lung);assert len(shared)==4 and len(primary-lung)==5 and len(lung-primary)==6
stats=[];paired_detail=[]
for typ,features in [('gene',genes),('program',list(panels))]:
 for feature in features:
  t=values[values.feature.eq(feature)].pivot(index='patient_id',columns='site_group',values='value');a=t.reindex(shared).primary_tumor.to_numpy();b=t.reindex(shared).lung_metastasis.to_numpy();good=np.isfinite(a)&np.isfinite(b);d=b[good]-a[good];Pval,test=signflip(d)
  stats.append(dict(contrast='paired_lung_minus_primary',feature_type=typ,feature=feature,n_primary=int(good.sum()),n_lung=int(good.sum()),mean_effect=float(np.mean(d)),median_effect=float(np.median(d)),positive_donors=int(np.sum(d>0)),P=Pval,test=test,family='paired_site_'+typ))
  for donor,x,y in zip(np.array(shared)[good],a[good],b[good]):paired_detail.append(dict(patient_id=donor,feature=feature,feature_type=typ,primary=x,lung=y,delta=y-x))
  aa=t.reindex(sorted(primary-lung)).primary_tumor.to_numpy();bb=t.reindex(sorted(lung-primary)).lung_metastasis.to_numpy();aa=aa[np.isfinite(aa)];bb=bb[np.isfinite(bb)];Pval,test=permute(aa,bb)
  stats.append(dict(contrast='independent_lung_minus_primary',feature_type=typ,feature=feature,n_primary=len(aa),n_lung=len(bb),mean_effect=float(np.mean(bb)-np.mean(aa)),median_effect=float(np.median(bb)-np.median(aa)),positive_donors=np.nan,P=Pval,test=test,family='independent_site_'+typ))
stats=pd.DataFrame(stats);stats['FDR']=stats.groupby('family').P.transform(lambda x:bh(x.to_numpy()));stats['nominal_P_lt_0_05']=stats.P.lt(.05);stats.to_csv(O/'tables/site_comparison_all_prespecified_features.tsv',sep='\t',index=False);pd.DataFrame(paired_detail).to_csv(O/'tables/paired_patient_feature_deltas.tsv',sep='\t',index=False)
ae=assoc.groupby(['patient_id','site_group','event_month','marker','program']).rho.mean().reset_index();av=ae.groupby(['patient_id','marker','program']).rho.mean().reset_index();av.to_csv(O/'tables/patient_marker_program_correlations.tsv',sep='\t',index=False)
ar=[]
for (marker,name),t in av.groupby(['marker','program']):
 v=t.rho.dropna().to_numpy();Pval,test=signflip(v);ar.append(dict(marker=marker,program=name,n_patients=len(v),median_rho=float(np.median(v)),mean_rho=float(np.mean(v)),positive_patients=int(np.sum(v>0)),P=Pval,test=test,family='marker_program_'+name))
ar=pd.DataFrame(ar);ar['FDR']=ar.groupby('family').P.transform(lambda x:bh(x.to_numpy()));ar['nominal_P_lt_0_05']=ar.P.lt(.05);ar.to_csv(O/'tables/marker_program_patient_level_tests.tsv',sep='\t',index=False)
roster=[]
for donor in sorted(primary|lung):roster.append(dict(patient_id=donor,role='paired' if donor in shared else 'independent_primary' if donor in primary else 'independent_lung'))
pd.DataFrame(roster).to_csv(O/'tables/site_comparison_patient_roster.tsv',sep='\t',index=False)
qa={'at':now(),'samples':26,'patients':19,'paired_patients':len(shared),'independent_primary':len(primary-lung),'independent_lung':len(lung-primary),'genes':len(genes),'programs':len(panels),'site_tests':len(stats),'marker_program_tests':len(ar),'no_spot_level_P_used':True,'no_new_QC_gate':True,'metadata_correction_used':True,'protocol_sha256':hashlib.sha256((O/'resources/VALIDATION_PROTOCOL_FROZEN.json').read_bytes()).hexdigest()}
(O/'delivery/PATIENT_VALIDATION_QA.json').write_text(json.dumps(qa,indent=2));(O/'PATIENT_VALIDATION_FINISHED').write_text(now());print(json.dumps(qa));print(stats[stats.feature.isin(markers)].to_string(index=False))
