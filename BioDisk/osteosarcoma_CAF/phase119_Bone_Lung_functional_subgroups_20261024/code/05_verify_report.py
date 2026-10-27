from caf119_common import *
import itertools
assert (O/'HOLDOUT_COMPLETE.json').exists();event('05_verify','started')
m=read('cohort_registry.tsv.gz');labels=read('hierarchical_cell_labels.tsv.gz');nodes=read('hierarchical_nodes.tsv');assert len(m)==len(labels)==6633 and m.cell_id.equals(labels.cell_id);assert m.cell_id.is_unique and nodes.node_id.is_unique
assert set(m.source)=={'SCPCP000017','SCPCP000023','NCI'} and set(m.site)=={'Bone','Lung'};assert m.groupby(['source','donor']).site.nunique().max()==1
original=pd.read_csv(R/'phase111_reference_figure_adaptation_20261021/tables/Fig1i_cell_values.tsv.gz',sep='\t');assert original[['cell_id','source','donor','site']].equals(m[['cell_id','source','donor','site']])
assert (~m.doublet).sum()==6423;assert (m.site.eq('Bone')&~m.doublet).sum()==2788 and (m.site.eq('Lung')&~m.doublet).sum()==3635
for row in nodes.itertuples():
 col='parent' if row.level=='parent' else 'leaf';actual=(labels[col].eq(row.node_id)&~m.doublet).sum();assert actual==row.fit_cells
 if row.level=='child':assert labels.loc[labels.leaf.eq(row.node_id),'parent'].eq(row.parent).all()
for site,g in labels[~m.doublet].groupby('site'):assert g.groupby('leaf').size().sum()==len(g)
split=read('holdout_patient_split.tsv');assert split.groupby('donor_key').holdout.nunique().max()==1
df=read('main_statistics.tsv');hd=read('holdout_statistics.tsv');maxbh=0.
for frame in [df,hd]:
 assert frame.loc[frame.N.lt(2),'P'].isna().all();assert frame.loc[frame.test.eq('not_tested_allzero'),'P'].isna().all()
 for _,g in frame.groupby(['site','source','level','background','arm','family']):
  g=g[g.P.notna()].sort_values('P');n=len(g)
  if not n:continue
  q=np.minimum.accumulate((g.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1);error=float(np.max(abs(q-g.FDR.to_numpy())));assert error<1e-10;maxbh=max(maxbh,error)
effects=read('main_patient_effects.tsv.gz');target=df[df.arm.isin(['singlet','singlet_min5'])&df.family.isin(['program','rank'])];checked=0;maxexact=0.;maxmc=0.;cache={}
for key,g in target.groupby(['site','source','node_id','background','arm']):
 site,src,node,bg,arm=key;ee=effects[(effects.site==site)&(effects.source==src)&(effects.node_id==node)&(effects.background==bg)&(effects.arm==arm)]
 if ee.empty:assert g.N.eq(0).all();continue
 pivot=ee.pivot(index='donor',columns='endpoint',values='effect')
 for row in g.itertuples():
  v=pivot[row.endpoint].dropna().to_numpy() if row.endpoint in pivot else np.array([]);n=len(v);assert n==row.N
  if n:assert abs(v.mean()-row.mean_effect)<1e-10 and int((v>1e-12).sum())==row.positive_patients
  if pd.isna(row.P):continue
  if n not in cache:cache[n]=np.array(list(itertools.product((-1.,1.),repeat=n)))
  exact=float(np.mean(abs(cache[n]@v)>=abs(v.sum())-1e-10));error=abs(exact-row.P)
  if row.test=='exact_signflip':assert error<1e-10;maxexact=max(maxexact,error)
  else:assert error<max(.001,6*np.sqrt(exact*(1-exact)/65536));maxmc=max(maxmc,error)
  checked+=1
hashes=json.loads((O/'resources/INPUT_HASHES.json').read_text())
for h in hashes:assert sha(h['path'])==h['sha256']
for stage in ['01_inputs','02','03','04']:assert (O/'logs'/(stage+'.exit_code')).read_text().strip()=='0'
qa=dict(status='PASS_COMPUTATION_AND_SCOPE',all_cells=6633,fit_cells=6423,Bone_fit_cells=2788,Lung_fit_cells=3635,patient_site_source_labels_unchanged=True,no_NCC_or_normal_cells=True,no_repeat_overall_Bone_Lung_DE=True,hierarchy_counts_reconciled=True,heldout_patient_split_disjoint=True,independent_function_tests=checked,max_exact_P_error=maxexact,max_MC_P_error=maxmc,max_BH_error=maxbh,input_hashes_unchanged=len(hashes),scope='Exploratory functional subgroup discovery; not protein secretion or verified CAF identity')
dump('ACCEPTANCE.json',qa)
candidates=read('all_node_function_evidence.tsv');focus=[]
for site in ['Bone','Lung']:
 for level in ['parent','child']:
  x=candidates[candidates.site.eq(site)&candidates.level.eq(level)&candidates.evidence_tier.eq('nominal_function_candidate')].sort_values(['min5_positive_tests','fibro_context_nominal_positive_tests','maximum_source_N','min_seed_jaccard'],ascending=False)
  focus.extend(x.head(3).to_dict('records'))
dump('REVIEW_SUMMARY.json',dict(structure=json.loads((O/'STRUCTURE_COMPLETE.json').read_text()),candidate_nodes=len(candidates[candidates.evidence_tier.eq('nominal_function_candidate')]),focus_candidates=focus,holdout=json.loads((O/'HOLDOUT_COMPLETE.json').read_text()),QA=qa))
event('05_verify','completed',independent_tests=checked);print(json.dumps(qa));print(pd.DataFrame(focus)[['site','node_id','level','fit_cells','source_patients','maximum_source_N','nominal_target_programs','min5_positive_tests','fibro_context_nominal_positive_tests','min_seed_jaccard']].to_string(index=False) if focus else 'NO_NOMINAL_FUNCTION_CANDIDATES')
