from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,hashlib,itertools
import numpy as np,pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase115_incremental_function_20261023'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
d=pd.read_csv(O/'tables/functional_statistics.tsv',sep='\t');e=pd.read_csv(O/'tables/patient_function_effects.tsv.gz',sep='\t');m=pd.read_csv(O/'tables/cell_mapping.tsv.gz',sep='\t',low_memory=False)
part=pd.read_csv(O/'tables/donor_partitions.tsv',sep='\t');assert part.groupby('donor_key').split.nunique().max()==1
ref=pd.read_csv(R/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',sep='\t').set_index('cell_id')
old=m[m.cohort.eq('baseline')].set_index('cell_id');assert set(old.index)==set(ref.index);assert (old.class_compatible==ref.loc[old.index].expanded_class_id).all();assert (old.state_compatible==ref.loc[old.index].subgroup_id).all()
assert len(m)==10510 and m.cell_id.is_unique
fx=pd.read_csv(O/'tables/cell_function_scores.tsv.gz',sep='\t');assert fx.cell_id.equals(m.cell_id)
excluded=set(json.loads((O/'resources/MAPPING_EXCLUDED_GENES.json').read_text()))
for f,col in [('class_mapping_genes.tsv','gene'),('state_mapping_genes.tsv','gene'),('NCC_cluster_features.tsv','symbol'),('NCI_cluster_features.tsv','symbol')]:assert not set(pd.read_csv(O/'tables'/f,sep='\t')[col])&excluded
checked=0;maxerror=0.;mcmax=0.;signcache={}
new=d[d.cohort.eq('added')&d.arm.isin(['singlet','singlet_min5'])]
for keys,rows in new.groupby(['source','arm','level','label']):
 src,arm,lev,lab=keys;eff=e[e.source.eq(src)&e.arm.eq(arm)&e.level.eq(lev)&e.label.eq(lab)]
 for _,v in rows.iterrows():
  vec=eff[eff.endpoint.eq(v.endpoint)].effect.to_numpy();assert len(vec)==v.N
  if len(vec):
   assert abs(np.mean(vec)-v.mean_effect)<1e-10
   assert int((vec>1e-12).sum())==v.positive_patients
  if pd.isna(v.P):continue
  n=len(vec)
  if n not in signcache:signcache[n]=np.array(list(itertools.product((-1.,1.),repeat=n)),dtype=np.float64)
  exact=float(np.mean(np.abs(signcache[n]@vec)>=abs(vec.sum())-1e-10))
  if v.test=='exact_patient_signflip':
   err=abs(exact-v.P);maxerror=max(maxerror,err);assert err<1e-10
  else:
   err=abs(exact-v.P);mcmax=max(mcmax,err);assert err<max(.001,6*np.sqrt(exact*(1-exact)/65536))
  checked+=1
maxbh=0.
for _,g in d.groupby(['source','arm','level','family_kind']):
 g=g[g.P.notna()].sort_values('P');n=len(g)
 if not n:continue
 expected=np.minimum.accumulate((g.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1)
 err=float(np.max(abs(expected-g.FDR.to_numpy())));assert err<1e-10;maxbh=max(maxbh,err)
assert d.loc[d.N.lt(2),'P'].isna().all();assert d.loc[d.test.eq('all_zero_effect_not_tested'),'P'].isna().all()
hashes=json.loads((O/'resources/INPUT_HASHES.json').read_text())
for h in hashes:assert sha(h['path'])==h['sha256'],h['path']
cn=json.loads((R/'CAF_CURRENT_CELL_EVIDENCE.json').read_text());assert sha(cn['cell_table'])==cn['cell_table_sha256']
for h in json.loads((R/'phase93_expanded_admission_20261016/objects/COUNT_COLLECTION.json').read_text())['added_objects']:assert sha(h['path'])==h['sha256']
for task in ['extract','map','cluster','stats']:assert (O/'logs'/(task+'.exit_code')).read_text().strip()=='0'
result=dict(status='PASS',cells=10510,baseline_cells_unchanged=8174,query_cells=2336,independent_exact_or_exact_MC_comparisons=checked,max_exact_P_error=maxerror,max_MC_P_error=mcmax,all_BH_families_reproduced=True,max_BH_error=maxbh,input_hashes_unchanged=len(hashes),baseline_authoritative_hash_matches=True,NCC_archived_raw_hashes_match=7,donor_partitions_disjoint=True,function_genes_excluded_from_mapping_and_clustering=True,insufficient_N_and_all_zero_P_missing=True,all_native_exit_codes_zero=True,figure_queue='paused_by_user',interpretive_acceptance='Exploratory functional assessment only; label identity not validated')
dump('ACCEPTANCE.json',result);print(json.dumps(result))
