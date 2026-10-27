from caf119_common import *
from sklearn.metrics import adjusted_rand_score
def discover(X,m,genes,site,fitrows,allrows,prefix):
 event(prefix,'root_fit',site=site,cells=len(fitrows))
 z,raw,pos=fit_rep(X,m,genes,fitrows,allrows,prefix+'_root');labs,graph=graph_and_labels(z[pos],.6);primary=np.array([encode(site,v) for v in labs]);out=m.iloc[allrows][['cell_id','site','source','donor','donor_key','doublet','cohort','original_state','old_lineage','fibro_supported']].copy()
 predicted=predict(z[pos],primary,z);predicted[pos]=primary;out['parent']=predicted;out['leaf']=predicted;out['is_fit']=np.isin(allrows,fitrows);nodes=[];eligibility=[];graphs={}
 alt={}
 for key,res,seed in [('res_low',.3,20261024),('res_high',1.,20261024),('seed',.6,20261025)]:alt[key]=graph_and_labels(z[pos],res,seed,graph)[0]
 alt['raw']=graph_and_labels(raw[pos],.6)[0]
 sparse.save_npz(O/'objects'/(prefix+'_root_graph.npz'),graph)
 for parent in sorted(set(primary)):
  mask=primary==parent;fr=fitrows[mask];ar=allrows[predicted==parent];pm=m.iloc[fr];dc=pm.groupby('donor_key').size();sc=pm.source.value_counts();base=dict(site=site,node_id=parent,parent=parent,level='parent',fit_cells=len(fr),all_assigned_cells=len(ar),source_patients=len(dc),sources=pm.source.nunique(),dominant_patient_fraction=float(dc.max()/len(fr)),dominant_source=sc.index[0],dominant_source_fraction=float(sc.iloc[0]/len(fr)),reference_cells=int(pm.cohort.eq('baseline').sum()),NCI_cells=int(pm.source.eq('NCI').sum()),fibro_supported_cells=int(pm.fibro_supported.sum()),stability_scope='root')
  for key,v in alt.items():base[key+'_jaccard']=best_jaccard(mask,v)
  base['parent_seed_jaccard']=base['seed_jaccard'];nodes.append(base)
  eligible=len(fr)>=80 and len(dc)>=3 and int((dc>=5).sum())>=2;rec=dict(parent=parent,fit_cells=len(fr),source_patients=len(dc),patients_with5=int((dc>=5).sum()),eligible=eligible)
  if not eligible:rec['outcome']='coverage_rule_no_further_split';eligibility.append(rec);continue
  event(prefix,'child_fit',parent=parent,cells=len(fr));cz,craw,cp=fit_rep(X,m,genes,fr,ar,prefix+'_'+parent,n_hvg=2000,npc=30);ll,cgraph=graph_and_labels(cz[cp],.5);rec['child_clusters']=len(set(ll))
  if len(set(ll))==1:rec['outcome']='one_cluster_no_forced_split';eligibility.append(rec);continue
  rec['outcome']='split_once';eligibility.append(rec);actual=np.array([parent+'_s'+v for v in ll]);pr=predict(cz[cp],actual,cz);pr[cp]=actual;out.loc[ar,'leaf']=pr
  aa={}
  for key,res,seed in [('res_low',.3,20261024),('res_high',.8,20261024),('seed',.5,20261025)]:aa[key]=graph_and_labels(cz[cp],res,seed,cgraph)[0]
  aa['raw']=graph_and_labels(craw[cp],.5)[0]
  for child in sorted(set(actual)):
   mask=actual==child;rr=fr[mask];mm=m.iloc[rr];dd=mm.groupby('donor_key').size();ss=mm.source.value_counts();b=dict(site=site,node_id=child,parent=parent,level='child',fit_cells=len(rr),all_assigned_cells=int((pr==child).sum()),source_patients=len(dd),sources=mm.source.nunique(),dominant_patient_fraction=float(dd.max()/len(rr)),dominant_source=ss.index[0],dominant_source_fraction=float(ss.iloc[0]/len(rr)),reference_cells=int(mm.cohort.eq('baseline').sum()),NCI_cells=int(mm.source.eq('NCI').sum()),fibro_supported_cells=int(mm.fibro_supported.sum()),stability_scope='within fixed parent; parent stability listed separately',parent_seed_jaccard=base['seed_jaccard'])
   for key,v in aa.items():b[key+'_jaccard']=best_jaccard(mask,v)
   nodes.append(b)
  sparse.save_npz(O/'objects'/(prefix+'_'+parent+'_graph.npz'),cgraph)
 save(out.reset_index(drop=True),prefix+'_labels.tsv.gz');save(pd.DataFrame(nodes),prefix+'_nodes.tsv');save(pd.DataFrame(eligibility),prefix+'_refinement_eligibility.tsv');dump('resources/'+prefix+'_root_stability.json',{k:float(adjusted_rand_score(primary,v)) for k,v in alt.items()});return out,pd.DataFrame(nodes)
