from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd,json,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=R/'phase115_incremental_function_20261023';O=R/'phase116_original33_analogs_20261023'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
d=read('all_33_by_19_comparisons.tsv');old=read('original_33_subgroup_dictionary.tsv');cells=read('cell_original33_retrieval.tsv.gz');assert len(d)==627 and d.old_state.nunique()==33 and d.query_group.nunique()==19
assert len(cells)==2336 and cells.cell_id.is_unique
orig=pd.read_csv(R/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv',sep='\t');assert old.subgroup_id.tolist()==orig.subgroup_id.tolist();assert old.working_name.equals(orig.working_name)
e=read('signature_patient_effects.tsv.gz');s=read('all33_signature_tests.tsv');checked=0;errmax=0;signs={}
for keys,g in e.groupby(['source','query_group','old_state','arm','minimum']):
 source,group,state,arm,minimum=keys;t=s[(s.source==source)&(s.query_group==group)&(s.old_state==state)&(s.arm==arm)&(s.minimum==minimum)].iloc[0];v=g.effect.to_numpy();n=len(v)
 assert n==t.N and abs(v.mean()-t.mean_effect)<1e-10 and int((v>1e-12).sum())==t.positive_patients
 if pd.isna(t.P):continue
 if n not in signs:signs[n]=1-2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8)
 p=float(np.mean(np.abs(signs[n]@v)>=abs(v.sum())-1e-10));err=abs(p-t.P);assert err<1e-10;errmax=max(err,err);checked+=1
bh=0.
for frame,groupcols in [(s,['source','arm','minimum']),(read('inherited_pathway_newgroup_statistics.tsv'),['source','minimum'])]:
 assert frame.loc[frame.N.lt(2),'P'].isna().all()
 for _,g in frame.groupby(groupcols):
  g=g[g.P.notna()].sort_values('P');n=len(g)
  if not n:continue
  expected=np.minimum.accumulate((g.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1);er=np.max(abs(expected-g.FDR.to_numpy()));assert er<1e-10;bh=max(bh,float(er))
for (_,state),g in d.groupby(['source','old_state']):assert (g.sort_values('marker_similarity',ascending=False).marker_rank_for_old_state.to_numpy()==np.arange(1,len(g)+1)).all()
hashes=json.loads((O/'resources/INPUT_HASHES.json').read_text())
for h in hashes:assert sha(h['path'])==h['sha256']
assert (O/'logs/compare.exit_code').read_text().strip()=='0'
qa=dict(status='PASS_EXPLORATORY_SUBGROUP_RETRIEVAL',old_subgroups=33,query_groups=19,all_pairwise_comparisons=627,query_cells=2336,no_parent_class_gate=True,original_names_and_labels_preserved=True,signature_tests_independently_recomputed=checked,max_exact_P_error=errmax,max_BH_error=bh,input_hashes_unchanged=len(hashes),no_new_scientific_figures=True,interpretation='Analog retrieval and conditional functional evidence; not one-to-one biological identity validation')
(O/'ACCEPTANCE.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa))
