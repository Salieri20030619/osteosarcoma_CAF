
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,itertools,json
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'audit_phase86_nonsignificance_20261012_v1';B=P/'phase86_functional_bone_lung_20261012'
def save(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
f=pd.read_csv(B/'tables/all_functional_bone_lung_statistics.tsv',sep='\t');c=pd.read_csv(B/'tables/cell_scores_and_original_labels.tsv.gz',sep='\t');th=pd.read_csv(B/'tables/state_thresholds.tsv',sep='\t');clinical=pd.read_csv(B/'tables/clinical_manifest.tsv',sep='\t')
keys=['scope','dataset','matrix','program'];threshold=th.set_index(keys).threshold.to_dict();patients=[];samples=[]
for kk,t in th.groupby(keys):
 scope,ds,matrix,pr=kk;g=c[c.dataset==ds]
 if scope=='Fibro_reviewed':g=g[g.review_lineage_code=='Fibro']
 col=matrix+'__'+pr;cut=threshold[kk]
 for (donor,site),v in g.groupby(['donor_record','strict_site']):
  ss=[]
  for sm,q in v.groupby('sample'):
   score=q[col].mean();frac=q[col].gt(cut).mean();ss.append((score,frac))
   samples.append(dict(scope=scope,dataset=ds,matrix=matrix,program=pr,donor_record=donor,site=site,sample=sm,cells=len(q),mean_score=score,high_fraction=frac))
  patients.append(dict(scope=scope,dataset=ds,matrix=matrix,program=pr,donor_record=donor,site=site,mean_score=np.mean([z[0] for z in ss]),high_fraction=np.mean([z[1] for z in ss]),samples=len(ss)))
a=pd.DataFrame(patients);save(a,'17_equal_sample_patient_endpoints.tsv');ss=pd.DataFrame(samples).merge(clinical[['sample','provider_disease_timing','treatment']],on='sample',validate='many_to_one');save(ss,'18_sample_scores_and_clinical_context.tsv')
cache={};rows=[]
for r in f.itertuples():
 g=a[(a.scope==r.scope)&(a.dataset==r.dataset)&(a.matrix==r.matrix)&(a.program==r.program)]
 def ids(s):return [] if pd.isna(s) else s.split(';')
 b=g[g.site=='Bone'].set_index('donor_record');l=g[g.site=='Lung'].set_index('donor_record');ib=ids(r.bone_patients);il=ids(r.lung_patients);x=b.loc[ib,r.endpoint].to_numpy(float);y=l.loc[il,r.endpoint].to_numpy(float);n,m=len(x),len(y);paired=r.design=='paired'
 effect=float(y.mean()-x.mean()) if n and m else np.nan;pp=np.nan
 if n>=2 and m>=2:
  z=y-x if paired else np.r_[x,y]
  if np.any(z!=0) if paired else np.any(z!=z[0]):
   key=(paired,n,m)
   if key not in cache:
    if paired:w=np.array(list(itertools.product([-1.,1.],repeat=n)))/n
    else:
     ii=list(itertools.combinations(range(n+m),m));w=np.full((len(ii),n+m),-1/n)
     for j,ix in enumerate(ii):w[j,list(ix)]=1/m
    cache[key]=w
   null=np.abs(cache[key]@z);pp=float(np.mean(null>=abs(effect)-1e-12))
 row={k:getattr(r,k) for k in ['scope','dataset','matrix','program','minimum_cells','design','endpoint','bone_N','lung_N','family']};row.update(original_effect=r.mean_lung_minus_bone,original_P=r.P,equal_sample_effect=effect,equal_sample_P=pp,weighting_diagnostic='all_samples_equal_within_patient_site;original_rosters_and_high_thresholds_fixed');rows.append(row)
r=pd.DataFrame(rows);r['equal_sample_FDR']=np.nan
for _,g in r.groupby('family'):
 ix=g.index[g.equal_sample_P.notna()]
 if len(ix):r.loc[ix,'equal_sample_FDR']=multipletests(r.loc[ix,'equal_sample_P'],method='fdr_bh')[1]
r['nominal_sensitivity']=r.equal_sample_P.lt(.05);save(r,'19_all800_equal_sample_weight_sensitivity.tsv')
main=r[(r.matrix=='raw')&(r.minimum_cells==1)&(r.design=='independent')&(r.endpoint=='mean_score')];save(main,'20_primary_equal_sample_weight_sensitivity.tsv')
# Clinical context counts restricted to actual primary patient rosters, rather than entire source manifests.
rosters=[]
for q in main[main.program=='collagen_assembly'].itertuples():
 old=f[(f.scope==q.scope)&(f.dataset==q.dataset)&(f.program==q.program)&(f.matrix=='raw')&(f.minimum_cells==1)&(f.design=='independent')&(f.endpoint=='mean_score')].iloc[0]
 for site,col in [('Bone','bone_patients'),('Lung','lung_patients')]:
  g=ss[(ss.scope==q.scope)&(ss.dataset==q.dataset)&(ss.matrix=='raw')&(ss.program==q.program)&(ss.site==site)&ss.donor_record.isin(ids(old[col]))]
  rosters.append(g)
rr=pd.concat(rosters);save(rr,'21_primary_sample_roster_clinical.tsv')
counts=rr.groupby(['scope','dataset','site','provider_disease_timing','treatment'],dropna=False).agg(samples=('sample','nunique'),donors=('donor_record','nunique'),cells=('cells','sum')).reset_index();save(counts,'22_primary_clinical_context_counts.tsv')
# Same original marker evidence remains intact, with its own target and matrix scale.
mk=pd.read_csv(P/'phase82_allCAF_crossdataset_markers_20261012/tables/primary_same_layer_crossdataset_markers.tsv',sep='\t');save(mk,'23_original14_crossdataset_markers.tsv')
summary={'all_valid':int(r.equal_sample_P.notna().sum()),'nominal':int(r.nominal_sensitivity.sum()),'FDR_lt005':int(r.equal_sample_FDR.lt(.05).sum()),'min_P':r.equal_sample_P.min(),'primary_raw_Fibro_direction_changes':int(((main.scope=='Fibro_reviewed')&(main.dataset=='SCPCP000017')&(main.original_effect*main.equal_sample_effect<0)).sum()),'purpose':'Post hoc diagnostic of within-patient sample weighting, not replacement primary analysis; tiny biopsies receive equal weight in sensitivity.'}
(O/'SAMPLE_WEIGHT_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
print('NOMINAL');print(r[r.nominal_sensitivity][['scope','dataset','matrix','program','minimum_cells','design','endpoint','original_effect','original_P','equal_sample_effect','equal_sample_P','equal_sample_FDR']].to_string(index=False))
print('MULTISAMPLE_LUNG');print(ss[(ss.scope=='Fibro_reviewed')&(ss.matrix=='raw')&ss.donor_record.str.endswith('CCG1106_007')&ss.program.isin(['collagen_assembly','contraction','IL6_JAK_STAT3'])].to_string(index=False))
print('PRIMARY_CLINICAL');print(counts[counts.scope=='Fibro_reviewed'].to_string(index=False))
