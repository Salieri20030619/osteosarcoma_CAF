from caf117_common import *
from sklearn.metrics import adjusted_rand_score
import scanpy as sc,anndata as ad
assert json.loads((O/'FREEZE_COMPLETE.json').read_text())['status']=='PASS'
event('02_structure','started');X,m,genes=load();fit=np.flatnonzero(~m.doublet);q=np.flatnonzero(m.doublet)
pc,hvg=representation(X,m,genes,fit,'joint');raw=pc[fit];fm=m.iloc[fit].reset_index(drop=True);representations={'unintegrated':raw};metrics=[]
def quality(emb,name):
 good=np.flatnonzero(fm.truth_lineage.ne('Unassessed'));z=fm.iloc[good].reset_index(drop=True);nn=neighbors(emb[good],15)
 lab=z.truth_lineage.to_numpy();source=z.source.to_numpy();purity=(lab[nn]==lab[:,None]).mean(1);cross=((source[nn]!=source[:,None])&(lab[nn]==lab[:,None])).mean(1)
 df=z[['source','donor','truth_lineage']].copy();df['purity']=purity;df['same_lineage_cross_source_fraction']=cross
 patient=df.groupby(['source','donor','truth_lineage'])[['purity','same_lineage_cross_source_fraction']].mean().reset_index();r=patient.groupby('truth_lineage')[['purity','same_lineage_cross_source_fraction']].mean().reset_index();r['representation']=name
 save(patient,name+'_pilot_patient_metrics.tsv');return r
metrics.append(quality(raw,'unintegrated'))
for theta in [1,2]:
 name='harmony_theta'+str(theta);event('02_structure','pilot_running',representation=name)
 z=harmony(raw,fm,theta);representations[name]=z;np.savez_compressed(O/'objects'/(name+'.npz'),cell_id=fm.cell_id.to_numpy(str),embedding=z);metrics.append(quality(z,name))
metrics=pd.concat(metrics,ignore_index=True);base=metrics[metrics.representation.eq('unintegrated')].set_index('truth_lineage');decisions=[];chosen='unintegrated'
for name in ['harmony_theta1','harmony_theta2']:
 cur=metrics[metrics.representation.eq(name)].set_index('truth_lineage');loss=base.purity-cur.purity;passed=bool(loss.mean()<=.05 and loss.max()<=.10)
 decisions.append(dict(representation=name,macro_purity_loss=float(loss.mean()),worst_lineage_loss=float(loss.max()),pass_operational_gate=passed))
 if chosen=='unintegrated' and passed:chosen=name
save(metrics,'pilot_representation_metrics.tsv');dump('PILOT_ACCEPTANCE.json',dict(status='PASS_COMPARISON',chosen=chosen,decisions=decisions,rule='First conservative source Harmony preserving known-lineage neighborhood purity; unintegrated fallback',not_biological_identity_validation=True))
event('02_structure','pilot_accepted',chosen=chosen)
z=representations[chosen];full=pc.copy();full[fit]=z
if len(q) and chosen!='unintegrated':full[q]=project(raw,z,pc[q],fm.source.to_numpy(),m.iloc[q].source.to_numpy())
labels,graph=clusters(z,.5,20261023);primary=relabel(labels);sparse.save_npz(O/'objects/joint_neighbor_graph.npz',graph)
out=m[['cell_id','source','donor','cohort','original_class','original_state','doublet','truth_lineage']].copy();out['joint_cluster']='';out.loc[fit,'joint_cluster']=primary
if len(q):out.loc[q,'joint_cluster']=assign_labels(z,primary,full[q])
sens={}
for name,res,seed in [('res03',.3,20261023),('res08',.8,20261023),('seed_repeat',.5,20261024)]:
 a=ad.AnnData(np.zeros((len(fit),0)));sc.tl.leiden(a,adjacency=graph,resolution=res,random_state=seed,key_added='labels',flavor='leidenalg',n_iterations=-1);v=a.obs.labels.astype(str).to_numpy();sens[name]=v;out[name]='';out.loc[fit,name]=v
labs_all,gall=clusters(full,.5,20261023)
out['allcell_sensitivity_cluster']=labs_all
# Compare no integration and alternative integration without selecting for a desired cluster count.
other='harmony_theta2' if chosen=='harmony_theta1' else 'harmony_theta1'
for name,emb in [('unintegrated',raw),(other,representations[other])]:
 v,_=clusters(emb,.5,20261023);sens[name]=v;out[name]='';out.loc[fit,name]=v
stability=[]
for label in sorted(set(primary)):
 mask=primary==label;ids=fit[mask];d=m.iloc[ids];dc=d.groupby(['source','donor']).size();scount=d.source.value_counts();rec=dict(joint_cluster=label,fit_cells=len(ids),source_patient_records=len(dc),sources=d.source.nunique(),dominant_patient_fraction=float(dc.max()/len(d)),dominant_source=scount.index[0],dominant_source_fraction=float(scount.iloc[0]/len(d)),original_cells=int(d.cohort.eq('baseline').sum()),new_cells=int(d.cohort.eq('added').sum()))
 for k,v in sens.items():rec[k+'_best_jaccard']=float(max(np.sum(mask&(v==t))/np.sum(mask|(v==t)) for t in set(v)))
 stability.append(rec)
save(pd.DataFrame(stability),'joint_cluster_stability.tsv');save(out,'joint_cell_labels.tsv.gz');np.savez_compressed(O/'objects/joint_representation.npz',cell_id=m.cell_id.to_numpy(str),raw_pca=pc,integrated=full,fit_rows=fit,features=hvg)
event('02_structure','baseline_control_running')
oldfit=np.flatnonzero(m.cohort.eq('baseline')&~m.doublet);bpc,_=representation(X,m,genes,oldfit,'baseline_control',fixed_features=hvg);bz=bpc[oldfit]
if chosen!='unintegrated':bz=harmony(bz,m.iloc[oldfit],int(chosen[-1]))
bl,_=clusters(bz,.5,20261023);save(pd.DataFrame(dict(cell_id=m.iloc[oldfit].cell_id,baseline_control_cluster=bl,joint_cluster=out.iloc[oldfit].joint_cluster)),'baseline_control_labels.tsv')
# Excluding original markers/core matrix-contractile genes is a sensitivity, not the primary discovery representation.
dictionary=read('original33_dictionary.tsv').fillna('');exclude=set(';'.join(dictionary.stable_markers).split(';'))
progs=json.loads((B/'resources/PROGRAMS.json').read_text())
for k in ['fibrillar_I_III_V','myCAF_customer','F2C3_collagen_support','collagen_crosslinking']:exclude.update(progs[k])
event('02_structure','marker_exclusion_sensitivity_running',excluded=len(exclude));spc,_=representation(X,m,genes,fit,'marker_exclusion',exclude=exclude);sz=spc[fit]
if chosen!='unintegrated':sz=harmony(sz,fm,int(chosen[-1]))
sl,_=clusters(sz,.5,20261023);save(pd.DataFrame(dict(cell_id=fm.cell_id,marker_exclusion_cluster=sl)),'marker_exclusion_labels.tsv')
stable=pd.DataFrame(stability)
for label in stable.joint_cluster:
 mask=primary==label;stable.loc[stable.joint_cluster.eq(label),'marker_exclusion_best_jaccard']=max(np.sum(mask&(sl==t))/np.sum(mask|(sl==t)) for t in set(sl))
save(stable,'joint_cluster_stability.tsv')
dump('STRUCTURE_COMPLETE.json',dict(status='PASS_COMPUTATION',selected_representation=chosen,all_cells=len(m),fit_cells=len(fit),clusters=len(set(primary)),baseline_control_clusters=len(set(bl)),baseline_to_joint_ARI=adjusted_rand_score(bl,out.iloc[oldfit].joint_cluster),sensitivity_ARI={k:adjusted_rand_score(primary,v) for k,v in sens.items()},allcell_ARI=adjusted_rand_score(primary,labs_all[fit]),marker_exclusion_ARI=adjusted_rand_score(primary,sl),functional_genes_allowed_in_primary=True,no_UMAP_or_figures=True,baseline_labels_preserved=True))
event('02_structure','completed',clusters=len(set(primary)),chosen=chosen)
