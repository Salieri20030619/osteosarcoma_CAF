from caf119_hierarchy import *
from caf119_statistics import evaluate
assert (O/'FUNCTION_COMPLETE.json').exists();event('04_holdout','started');X,m,genes=load();main=read('hierarchical_cell_labels.tsv.gz');nodes=read('hierarchical_nodes.tsv');score=read('cell_function_scores.tsv.gz');splits={}
for (site,src),g in m[~m.doublet].groupby(['site','source']):
 ds=sorted(g.donor_key.unique(),key=lambda d:hashlib.sha256(('Phase119_holdout|'+d).encode()).hexdigest());n=max(1,round(len(ds)*.2)) if len(ds)>=3 else 0
 for i,d in enumerate(ds):splits[d]='holdout' if i<n else 'discovery'
m['holdout']=m.donor_key.map(splits).eq('holdout');save(m[['site','source','donor','donor_key','holdout']].drop_duplicates(),'holdout_patient_split.tsv');assert m.groupby('donor_key').holdout.nunique().max()==1
pred=main.copy();bridges=[];diag=[]
for site in ['Bone','Lung']:
 allrows=np.flatnonzero(m.site.eq(site));fitrows=np.flatnonzero(m.site.eq(site)&~m.doublet&~m.holdout);held=np.flatnonzero(m.site.eq(site)&~m.doublet&m.holdout);disc,dnodes=discover(X,m,genes,site,fitrows,allrows,'holdout_'+site)
 for col in ['parent','leaf']:
  translation={}
  for label in disc.loc[fitrows,col].unique():
   rr=fitrows[disc.loc[fitrows,col].eq(label).to_numpy()];vc=main.loc[rr,col].value_counts();f=vc.iloc[0]/len(rr);mapped=vc.index[0] if f>=.5 else 'unaligned_'+label;translation[label]=mapped;bridges.append(dict(site=site,level=col,discovery_label=label,mapped_main_label=mapped,training_cells=len(rr),overlap=float(f)))
  pred.loc[allrows,col]=disc[col].map(translation).fillna('unaligned').to_numpy()
 for donor,g in m.iloc[held].groupby('donor_key'):
  ix=g.index;diag.append(dict(site=site,source=g.source.iloc[0],donor=g.donor.iloc[0],cells=len(ix),parent_agreement=float(pred.loc[ix,'parent'].eq(main.loc[ix,'parent']).mean()),leaf_agreement=float(pred.loc[ix,'leaf'].eq(main.loc[ix,'leaf']).mean())))
save(pred,'heldout_projected_labels.tsv.gz');save(pd.DataFrame(bridges),'discovery_to_main_bridges.tsv');save(pd.DataFrame(diag),'heldout_structure_diagnostic.tsv');df=evaluate(m,pred,nodes,score,'holdout',holdout=True)
dump('resources/HOLDOUT_METHOD.json',dict(discovery_only_HVG_PCA_Harmony_clustering=True,query_projection='Same-source30-neighbor learned displacement then15-neighbor distance vote',native_Harmony_query_API=False,bridge='Training-cell majority overlap >=0.5 otherwise unaligned',discovery_child_labels_with_unaligned_parent='Recorded; within-parent comparisons with unmatched parents can be unassessed',diagnostic_vs_full_labels='Transductive agreement only, not ground-truth accuracy',no_independent_external_validation=True))
dump('HOLDOUT_COMPLETE.json',dict(status='PASS_COMPUTATION_WITH_REPLICATION_LIMITS',heldout_cells=int((m.holdout&~m.doublet).sum()),heldout_patient_records=m[m.holdout&~m.doublet].groupby(['site','source']).donor.nunique().to_dict().__str__(),valid_tests=int(df.P.notna().sum()),no_test_patient_used_in_discovery=True));event('04_holdout','completed',heldout_cells=int((m.holdout&~m.doublet).sum()))
