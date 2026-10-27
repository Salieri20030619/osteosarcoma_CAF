
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json,math,itertools,hashlib
from scipy.stats import spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase86_functional_bone_lung_20261012';O=P/'audit_phase86_nonsignificance_20261012_v1'
def read(n):return pd.read_csv(B/'tables'/n,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
def dump(d,n): (O/n).write_text(json.dumps(d,ensure_ascii=False,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),encoding='utf-8')
def ids(s):return [] if pd.isna(s) else s.split(';')
f=read('all_functional_bone_lung_statistics.tsv');pat=read('patient_site_programs.tsv');c=read('cell_scores_and_original_labels.tsv.gz');cov=read('program_coverage.tsv');programs=json.loads((B/'resources/PROGRAMS.json').read_text())
inputs=[B/'tables'/n for n in ['all_functional_bone_lung_statistics.tsv','patient_site_programs.tsv','cell_scores_and_original_labels.tsv.gz','program_coverage.tsv']]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs};dump(hashes,'INPUT_SHA256.json')
# Enumerate the exact null independently for every valid original functional comparison.
rows=[];patient_rows=[];cache={};allerr=[]
for r in f.itertuples():
 g=pat[(pat.scope==r.scope)&(pat.dataset==r.dataset)&(pat.matrix==r.matrix)&(pat.program==r.program)&(pat.cells>=r.minimum_cells)]
 b=g[g.site=='Bone'].set_index('donor_record');l=g[g.site=='Lung'].set_index('donor_record')
 ib=ids(r.bone_patients);il=ids(r.lung_patients);x=b.loc[ib,r.endpoint].to_numpy(float);y=l.loc[il,r.endpoint].to_numpy(float)
 paired=r.design=='paired'
 assert set(ib)==set(il) if paired else not(set(ib)&set(il))
 assert len(x)==r.bone_N and len(y)==r.lung_N
 if paired:
  valid=np.isfinite(x)&np.isfinite(y);x=x[valid];y=y[valid]
 else:x=x[np.isfinite(x)];y=y[np.isfinite(y)]
 nx,ny=len(x),len(y);effect=y.mean()-x.mean() if nx and ny else np.nan
 z=np.r_[x,y];null=np.array([]);theoretical=np.nan
 if nx>=2 and ny>=2:
  key=(paired,nx,ny)
  if key not in cache:
   if paired:w=np.array(list(itertools.product([-1.,1.],repeat=nx)))/nx
   else:
    perm=list(itertools.combinations(range(nx+ny),ny));w=np.full((len(perm),nx+ny),-1/nx)
    for i,ix in enumerate(perm):w[i,list(ix)]=1/ny
   cache[key]=w
  null=np.abs(cache[key]@(y-x if paired else z));obs=abs(effect)
  theoretical=(2/2**nx if paired else (2 if nx==ny else 1)/math.comb(nx+ny,ny))
  actual_floor=float(np.mean(null>=null.max()-1e-12))
  variable=np.any(y-x!=0) if paired else np.any(z!=z[0])
  pp=float(np.mean(null>=obs-1e-12)) if variable else np.nan
 else:actual_floor=np.nan;pp=np.nan
 if pd.notna(r.P):allerr.append(abs(pp-r.P));assert abs(pp-r.P)<1e-10
 else:assert np.isnan(pp)
 if paired:
  d=y-x;sd=float(np.std(d,ddof=1)) if nx>1 else np.nan
  standardized=effect/sd if sd>0 else np.nan
  auc=np.nan
 else:
  sd=float(np.sqrt(((nx-1)*np.var(x,ddof=1)+(ny-1)*np.var(y,ddof=1))/(nx+ny-2))) if nx>1 and ny>1 else np.nan
  standardized=effect/sd if sd>0 else np.nan
  auc=float(np.mean((y[:,None]>x).astype(float)+.5*(y[:,None]==x))) if nx and ny else np.nan
 row=r._asdict();row.update(P_recomputed=pp,permutations=len(null),design_P_lower_bound=theoretical,pool_min_attainable_P=actual_floor,design_cannot_reach_005=bool(theoretical>=.05) if np.isfinite(theoretical) else None,pool_cannot_reach_005=bool(actual_floor>=.05) if np.isfinite(actual_floor) else None,bone_SD=np.std(x,ddof=1) if nx>1 else np.nan,lung_SD=np.std(y,ddof=1) if ny>1 else np.nan,standardized_effect=standardized,lung_greater_bone_pair_fraction=auc)
 rows.append(row)
 if r.matrix=='raw' and r.minimum_cells==1 and r.endpoint=='mean_score':
  for site,frame,ii in [('Bone',b,ib),('Lung',l,il)]:
   for donor in ii:
    v=frame.loc[donor];patient_rows.append(dict(scope=r.scope,dataset=r.dataset,design=r.design,program=r.program,site=site,donor=donor,cells=v.cells,score=v.mean_score,detected_features=v.detected_features_mean,mito=v.mitochondrial_fraction_mean,samples=v.samples))
a=pd.DataFrame(rows);save(a,'01_all800_test_diagnostics.tsv');save(pd.DataFrame(patient_rows),'02_primary_patient_roster_scores.tsv')
main=a[(a.matrix=='raw')&(a.minimum_cells==1)&(a.design=='independent')&(a.endpoint=='mean_score')].copy();save(main,'03_primary_effects_and_resolution.tsv')
design=a.groupby(['scope','dataset','matrix','minimum_cells','design','endpoint','bone_N','lung_N'],dropna=False).agg(rows=('program','size'),valid_P=('P','count'),min_P=('P','min'),floor=('design_P_lower_bound','min'),pool_floor_min=('pool_min_attainable_P','min')).reset_index();save(design,'04_sample_design_resolution.tsv')
# Actual cell and donor attrition at fixed identity stages; do not reclassify.
meta=pd.read_csv(P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz',sep='\t',low_memory=False)
keep=meta.dataset.isin(['SCPCP000017','SCPCP000023'])&meta.strict_site.isin(['Bone','Lung'])
dbl=meta.doublet_positive_any.astype(str).str.lower().isin(['true','1','1.0'])
attr=[]
for stage,mask in [('fixed_candidate_pool',keep),('doublet_negative',keep&~dbl),('Fibro_reviewed',keep&~dbl&meta.review_lineage_code.eq('Fibro'))]:
 t=meta[mask].groupby(['dataset','strict_site']).agg(cells=('cell_id','size'),donors=('donor_record','nunique'),samples=('sample','nunique')).reset_index();t['stage']=stage;attr.append(t)
save(pd.concat(attr),'05_cell_patient_attrition.tsv')
save(meta[keep&~dbl].groupby(['dataset','strict_site','review_lineage_code','analysis_class']).agg(cells=('cell_id','size'),patients=('donor_record','nunique')).reset_index(),'06_lineage_class_composition.tsv')
cov['coverage_fraction']=cov.measured/cov.defined;save(cov.groupby(['dataset','program']).agg(samples=('sample','size'),minimum_coverage=('coverage_fraction','min'),maximum_coverage=('coverage_fraction','max'),defined_genes=('defined','first')).reset_index(),'07_program_coverage_summary.tsv')
# Check score aggregation against exported cell scores, including unequal sample weighting.
aggregation=[];sample_effects=[]
for r in pat.itertuples():
 t=c[(c.dataset==r.dataset)&(c.donor_record==r.donor_record)&(c.strict_site==r.site)]
 if r.scope=='Fibro_reviewed':t=t[t.review_lineage_code=='Fibro']
 v=t[r.matrix+'__'+r.program];vm=v.mean();assert len(t)==r.cells
 assert abs(vm-r.mean_score)<1e-10 or pd.isna(vm) and pd.isna(r.mean_score)
 samplemean=t.groupby('sample')[r.matrix+'__'+r.program].mean().mean()
 aggregation.append(dict(scope=r.scope,dataset=r.dataset,matrix=r.matrix,program=r.program,donor_record=r.donor_record,site=r.site,cells=len(t),samples=t['sample'].nunique(),score_error=abs(vm-r.mean_score),mean_score=vm,equal_sample_score=samplemean,weighting_difference=samplemean-vm))
ag=pd.DataFrame(aggregation);save(ag,'08_patient_aggregation_audit.tsv')
for r in main.itertuples():
 t=ag[(ag.scope==r.scope)&(ag.dataset==r.dataset)&(ag.matrix=='raw')&(ag.program==r.program)]
 x=t[(t.site=='Bone')&t.donor_record.isin(ids(r.bone_patients))].equal_sample_score.to_numpy();y=t[(t.site=='Lung')&t.donor_record.isin(ids(r.lung_patients))].equal_sample_score.to_numpy()
 sample_effects.append(dict(scope=r.scope,dataset=r.dataset,program=r.program,original_effect=r.mean_lung_minus_bone,equal_sample_effect=y.mean()-x.mean(),note='descriptive_weighting_diagnostic_not_new_primary_test'))
save(pd.DataFrame(sample_effects),'09_equal_sample_weighting_effects.tsv')
# Technical covariance, correlation among gene sets, and raw/corrected effect stability.
tech=read('program_technical_covariation.tsv');save(tech,'10_technical_covariation_original.tsv')
co=[]
for (scope,ds),g in pat[pat.matrix=='raw'].groupby(['scope','dataset']):
 wide=g.pivot(index=['donor_record','site'],columns='program',values='mean_score');corr=wide.corr(method='spearman')
 for i,p1 in enumerate(corr.columns):
  for p2 in corr.columns[i+1:]:co.append(dict(scope=scope,dataset=ds,program1=p1,program2=p2,spearman=corr.loc[p1,p2],shared_genes=len(set(programs[p1])&set(programs[p2])),gene_jaccard=len(set(programs[p1])&set(programs[p2]))/len(set(programs[p1])|set(programs[p2]))))
save(pd.DataFrame(co),'11_program_correlation_overlap.tsv')
keys=['scope','dataset','program','minimum_cells','design','endpoint']
rc=a[a.matrix=='raw'].merge(a[a.matrix=='corrected'],on=keys,suffixes=('_raw','_corrected'));rc['effect_direction_same']=rc.mean_lung_minus_bone_raw*rc.mean_lung_minus_bone_corrected>0;save(rc,'12_raw_corrected_diagnostics.tsv')
# Gene-wise contribution to the identical original program mean: effect cancellation diagnostic.
mz=np.load(P/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz');mi=pd.Index(mz['cells'].astype(str));mg={v:i for i,v in enumerate(mz['genes'].astype(str))}
vals=np.log1p(mz['raw']/mz['raw_library'][:,None]*1e4);vals[~mz['assayed'].astype(bool)]=np.nan
contrib=[];cancell=[]
for r in main.itertuples():
 if not all(g in mg for g in programs[r.program]):continue
 t=c[c.dataset==r.dataset]
 if r.scope=='Fibro_reviewed':t=t[t.review_lineage_code=='Fibro']
 js=[mg[g] for g in programs[r.program]];means={}
 for site,ii in [('Bone',ids(r.bone_patients)),('Lung',ids(r.lung_patients))]:
  vv=[]
  for donor in ii:
   cells=t[(t.strict_site==site)&(t.donor_record==donor)].cell_id;ix=mi.get_indexer(cells);assert (ix>=0).all();v=vals[ix][:,js]
   # Score contribution per gene uses each cell's measured-gene denominator.
   n=np.isfinite(v).sum(1);vv.append(np.nansum(v/n[:,None],axis=0)/len(v))
  means[site]=np.array(vv).mean(0)
 delta=means['Lung']-means['Bone'];assert abs(delta.sum()-r.mean_lung_minus_bone)<1e-10
 for gene,d in zip(programs[r.program],delta):contrib.append(dict(scope=r.scope,dataset=r.dataset,program=r.program,gene=gene,mean_score_effect_contribution=d))
 totalabs=np.abs(delta).sum();cancell.append(dict(scope=r.scope,dataset=r.dataset,program=r.program,genes=len(delta),positive_genes=int((delta>0).sum()),negative_genes=int((delta<0).sum()),net_effect=delta.sum(),sum_absolute_contributions=totalabs,cancellation_fraction=1-abs(delta.sum())/totalabs if totalabs else np.nan))
save(pd.DataFrame(contrib),'13_gene_effect_contributions.tsv');save(pd.DataFrame(cancell),'14_program_cancellation.tsv')
# Main phenotype comparisons across data sources and original marker selection scope.
cross=main[main.dataset=='SCPCP000017'].merge(main[main.dataset=='SCPCP000023'],on=['scope','program'],suffixes=('_17','_23'));cross['same_direction']=cross.mean_lung_minus_bone_17*cross.mean_lung_minus_bone_23>0;save(cross,'15_cross_source_directions.tsv')
all14=json.loads((P/'CAF_ALL14_RESEARCH_CURRENT.json').read_text())['retained_genes'];save(pd.DataFrame([dict(gene=g,program_membership=';'.join(pr for pr,genes in programs.items() if g in genes),retained=True) for g in all14]),'16_all14_program_overlap.tsv')
ass=read('all_marker_functional_associations.tsv.gz')
summary={'all_registered':len(f),'valid_tests':int(f.P.notna().sum()),'nominal_tests':int(f.P.lt(.05).sum()),'global_min_P':float(f.P.min()),'independently_reproduced_P':len(allerr),'max_P_error':max(allerr),'sample_size_floor_blocked_valid_tests':int((a.P.notna()&a.design_cannot_reach_005.eq(True)).sum()),'pool_floor_blocked_valid_tests':int((a.P.notna()&a.pool_cannot_reach_005.eq(True)).sum()),'attainable_valid_tests':int((a.P.notna()&a.pool_cannot_reach_005.eq(False)).sum()),'nonevaluable_reasons':f[f.P.isna()].groupby('test').size().to_dict(),'min_P_by_scope_dataset':f.groupby(['scope','dataset']).P.min().to_dict(),'main_Fibro_LOPO_stable':int(main[main.scope=='Fibro_reviewed'].LOPO_direction_retained.eq(True).sum()),'main_Fibro_cross_source_same_direction':int(cross[cross.scope=='Fibro_reviewed'].same_direction.sum()),'all14_in_any_program':sum(any(g in genes for genes in programs.values()) for g in all14),'marker_association_valid_tests':int(ass.P.notna().sum()),'marker_association_nominal_tests':int(ass.P.lt(.05).sum()),'aggregation_max_error':float(ag.score_error.max()),'multiple_sample_patient_site_records':int(ag.samples.gt(1).sum()),'original_inputs_unchanged':all(hashlib.sha256(Path(k).read_bytes()).hexdigest()==v for k,v in hashes.items())}
summary['min_P_by_scope_dataset']={'|'.join(k):v for k,v in summary['min_P_by_scope_dataset'].items()}
dump(summary,'DIAGNOSTIC_SUMMARY.json')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print('PRIMARY_FIBRO');print(main[main.scope=='Fibro_reviewed'][['dataset','program','bone_N','lung_N','mean_lung_minus_bone','standardized_effect','P','LOPO_direction_retained','pool_min_attainable_P']].to_string(index=False))
print('DESIGNS');print(design[design.matrix=='raw'][['scope','dataset','minimum_cells','design','endpoint','bone_N','lung_N','valid_P','min_P','floor']].to_string(index=False))
