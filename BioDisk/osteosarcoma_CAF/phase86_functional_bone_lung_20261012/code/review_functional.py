from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,itertools,math
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase86_functional_bone_lung_20261012'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def exact(x,y,paired):
 x=np.array(x,float);y=np.array(y,float)
 if paired:
  ix=np.isfinite(x)&np.isfinite(y);d=y[ix]-x[ix]
  if len(d)<2 or not np.any(d!=0):return np.nan
  null=np.array([abs(np.mean(d*np.array(s))) for s in itertools.product([-1,1],repeat=len(d))]);obs=abs(d.mean())
 else:
  x=x[np.isfinite(x)];y=y[np.isfinite(y)];z=np.r_[x,y]
  if min(len(x),len(y))<2 or not np.any(z!=z[0]):return np.nan
  obs=abs(y.mean()-x.mean());ss=z.sum();null=np.array([abs(z[list(ix)].sum()/len(y)-(ss-z[list(ix)].sum())/len(x)) for ix in itertools.combinations(range(len(z)),len(y))])
 return float(np.mean(null>=obs-1e-12))
assert (O/'logs/analysis.exit_code').read_text().strip()=='0'
done=json.loads((O/'COMPUTATION_COMPLETE.json').read_text());assert done['original_inputs_unchanged']
f=read('all_functional_bone_lung_statistics.tsv');p=read('patient_site_programs.tsv');a=read('all_marker_functional_associations.tsv.gz');mp=read('marker_patient_effects.tsv.gz');cells=read('cell_scores_and_original_labels.tsv.gz')
assert not cells.doublet_positive_any.astype(str).str.lower().eq('true').any()
assert set(f.dataset)=={'SCPCP000017','SCPCP000023'}
errors=[]
for r in f.itertuples():
 g=p[(p.scope==r.scope)&(p.dataset==r.dataset)&(p.matrix==r.matrix)&(p.program==r.program)&(p.cells>=r.minimum_cells)]
 def ids(s):return [] if pd.isna(s) else str(s).split(';')
 b=g[g.site=='Bone'].set_index('donor_record');l=g[g.site=='Lung'].set_index('donor_record');ib=ids(r.bone_patients);il=ids(r.lung_patients)
 x=b.loc[ib,r.endpoint].to_numpy();y=l.loc[il,r.endpoint].to_numpy()
 assert len(x)==r.bone_N and len(y)==r.lung_N
 assert (set(ib)==set(il)) if r.design=='paired' else not (set(ib)&set(il))
 if r.test.startswith('exact_'):
  value=exact(x,y,r.design=='paired');assert abs(value-r.P)<1e-10,(r.program,value,r.P);errors.append(abs(value-r.P))
for df in [f,a]:
 for fam,g in df.groupby('family'):
  t=g[g.P.notna()]
  if len(t):assert np.max(np.abs(multipletests(t.P,method='fdr_bh')[1]-t.FDR))<1e-10
qa_markers=[]
# Diverse raw-program and cell-scope marker checks, including retained candidates.
selected=pd.concat([a[a.P.notna()].groupby(['scope','dataset','program'],group_keys=False).head(1),a[a.retained14&a.P.notna()].groupby(['gene','scope','dataset'],group_keys=False).head(1)]).drop_duplicates(['scope','dataset','program','gene'])
for r in selected.itertuples():
 v=mp[(mp.scope==r.scope)&(mp.dataset==r.dataset)&(mp.program==r.program)&(mp.gene==r.gene)].high_minus_low.to_numpy()
 assert len(v)==r.N and abs(v.mean()-r.mean_high_minus_low)<1e-10
 if r.test.startswith('exact_'):
  ev=exact(np.zeros(len(v)),v,True);assert abs(ev-r.P)<1e-10;qa_markers.append(abs(ev-r.P))
reg=read('ALL14_retention_registry.tsv');assert len(reg)==14 and reg.retain.all()
assert a.loc[a.component_of_score,'P'].isna().all()
main=f[(f.matrix=='raw')&(f.minimum_cells==1)&(f.design=='independent')].copy()
keys=['scope','program','endpoint']
s17=main[main.dataset=='SCPCP000017'];s23=main[main.dataset=='SCPCP000023']
cross=s17.merge(s23,on=keys,suffixes=('_17','_23'),validate='one_to_one')
cross['same_direction']=cross.mean_lung_minus_bone_17*cross.mean_lung_minus_bone_23>0
cross['both_nominal_same']=cross.same_direction&cross.P_17.lt(.05)&cross.P_23.lt(.05)
cross.to_csv(O/'tables/cross_source_functional_evidence.tsv',sep='\t',index=False)
# Preserve both directions; a positive high-state marker is a possible sorting proxy.
chain=a.merge(main[['scope','dataset','program','endpoint','mean_lung_minus_bone','P','FDR','bone_N','lung_N']],on=['scope','dataset','program'],suffixes=('_marker','_clinical'),validate='many_to_many')
chain['exploratory_positive_proxy']=chain.P_marker.lt(.05)&chain.mean_high_minus_low.gt(0)&chain.P_clinical.lt(.05)&chain.surface_annotation
chain.to_csv(O/'tables/clinical_function_surface_evidence_chain.tsv.gz',sep='\t',index=False)
short=chain[chain.exploratory_positive_proxy].copy();short.to_csv(O/'tables/exploratory_surface_proxies.tsv',sep='\t',index=False)
receipt=dict(status='PASS',functional_exact_P_checks=len(errors),marker_exact_P_checks=len(qa_markers),max_P_error=max(errors+qa_markers+[0]),all_FDR_families_reproduced=True,all14_retained=True,no_incomplete_source_used=True,no_score_component_used_as_independent_marker=True,clinical_patient_rosters_checked=True)
(O/'QA/INDEPENDENT_QA.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print('QA',json.dumps(receipt))
cols=['scope','dataset','program','endpoint','bone_N','lung_N','mean_lung_minus_bone','P','FDR','LOPO_direction_retained']
print('NOMINAL_CLINICAL\n'+main[main.nominal][cols].to_string(index=False))
print('CROSS_NOMINAL\n'+cross[cross.both_nominal_same][keys+['mean_lung_minus_bone_17','P_17','mean_lung_minus_bone_23','P_23']].to_string(index=False))
cols2=['scope','dataset','program','endpoint','gene','N','mean_high_minus_low','P_marker','FDR_marker','mean_lung_minus_bone','P_clinical','FDR_clinical']
print('PROXIES_COUNTS',short.groupby(['scope','dataset','program','endpoint']).size().to_dict())
print('FIBRO_PROXIES\n'+short[short.scope=='Fibro_reviewed'].sort_values(['program','P_marker']).groupby(['program','endpoint'],group_keys=False).head(12)[cols2].to_string(index=False))
print('RETAINED14\n'+reg.to_string(index=False))
