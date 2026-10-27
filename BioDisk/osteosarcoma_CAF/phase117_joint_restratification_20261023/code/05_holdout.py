from caf117_function_utils import *
from sklearn.metrics import adjusted_rand_score
assert (O/'FUNCTION_COMPLETE.json').exists();event('05_holdout','started');X,m,genes=load();full=read('joint_cell_labels.tsv.gz').set_index('cell_id');m['full_joint_cluster']=m.cell_id.map(full.joint_cluster);fitmask=~m.doublet
split={}
for src,g in m[fitmask].groupby('source'):
 ds=sorted(g.donor_key.unique(),key=lambda x:hashlib.sha256(('holdout117|'+x).encode()).hexdigest());n=max(1,round(len(ds)*.2)) if len(ds)>=3 else 0
 for i,d in enumerate(ds):split[d]='holdout' if i<n else 'discovery'
m['validation_split']=m.donor_key.map(split);train=np.flatnonzero(fitmask&m.validation_split.eq('discovery'));test=np.flatnonzero(fitmask&m.validation_split.eq('holdout'));save(m[['source','donor','donor_key','validation_split']].drop_duplicates(),'holdout_donor_split.tsv')
assert not set(m.iloc[train].donor_key)&set(m.iloc[test].donor_key)
pc,hvg=representation(X,m,genes,train,'holdout_discovery');selection=json.loads((O/'PILOT_ACCEPTANCE.json').read_text())['chosen'];tr=pc[train]
if selection!='unintegrated':tr=harmony(tr,m.iloc[train],int(selection[-1]))
dl,_=clusters(tr,.5,20261023);mapping=[];bridge={}
for c in sorted(set(dl)):
 ix=train[dl==c];vc=m.iloc[ix].full_joint_cluster.value_counts();fraction=vc.iloc[0]/len(ix);label=vc.index[0] if fraction>=.5 else 'D'+c+'_unaligned';bridge[c]=label;mapping.append(dict(discovery_cluster=c,mapped_joint_cluster=label,training_cells=len(ix),training_overlap_fraction=float(fraction)))
save(pd.DataFrame(mapping),'holdout_discovery_to_joint_bridge.tsv')
qt=pc[test]
if selection!='unintegrated':qt=project(pc[train],tr,pc[test],m.iloc[train].source.to_numpy(),m.iloc[test].source.to_numpy())
pred=assign_labels(tr,dl,qt);predj=np.array([bridge[c] for c in pred]);hm=m.iloc[test].copy();hm['joint_cluster']=predj;hm['discovery_cluster']=pred;hm['full_data_label_agreement']=hm.joint_cluster.eq(hm.full_joint_cluster)
save(hm[['cell_id','source','donor','joint_cluster','discovery_cluster','full_joint_cluster','full_data_label_agreement']],'heldout_cell_predictions.tsv.gz');save(hm.groupby(['source','donor']).agg(cells=('cell_id','size'),agreement_with_full_data_label=('full_data_label_agreement','mean')).reset_index(),'heldout_structural_agreement.tsv')
dump('resources/HOLDOUT_PROJECTION_METHOD.json',dict(method='Discovery-only HVG/PCA and Harmony; query projected with discovery PCA; for each query cell, interpolate learned corrected-minus-raw displacement from30 nearest discovery cells of the same source in PCA space; then15-NN distance-weighted discovery labels',native_Harmony_query_transform=False,no_query_data_used_to_fit_HVG_PCA_or_Harmony=True,joint_name_bridge='Training-cell majority overlap >=0.5; otherwise D_unaligned',test_label_comparison='Agreement with full-data labels is a transductive structural diagnostic, not ground-truth accuracy',function_test='Heldout patients assigned from discovery; expression values never used to select outcome-specific threshold',biological_validation='Exploratory within-study donor holdout, not independent external experiment'))
v,names,family=endpoints(m,X,genes);df=evaluate(hm,v,names,family,'joint_cluster',['singlet'],'holdout_function')
main=read('joint_function_statistics.tsv');main=main[main.arm.eq('singlet')];compare=main.merge(df,on=['source','joint_cluster','endpoint','family'],suffixes=('_full','_holdout'),how='left');compare['same_effect_direction']=np.sign(compare.mean_effect_full)==np.sign(compare.mean_effect_holdout);save(compare,'full_vs_holdout_function_comparison.tsv')
dump('HOLDOUT_COMPLETE.json',dict(status='PASS_COMPUTATION_WITH_REPLICATION_LIMITS',discovery_cells=len(train),heldout_cells=len(test),heldout_source_patient_records=m.iloc[test].groupby('source').donor.nunique().to_dict(),discovery_groups=len(set(dl)),unaligned_groups=sum('unaligned' in v for v in bridge.values()),valid_tests=int(df.P.notna().sum()),patient_split_disjoint=True,feature_selection_discovery_only=True,independent_external_validation=False))
event('05_holdout','completed',discovery_cells=len(train),heldout_cells=len(test))
