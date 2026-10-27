from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,re,datetime
import numpy as np,pandas as pd
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase115_incremental_function_20261023'
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
X=sparse.load_npz(O/'objects/log1p_cp10k_shared.npz');m=pd.read_csv(O/'tables/cell_registry.tsv.gz',sep='\t').fillna('');genes=pd.read_csv(O/'tables/shared_genes.tsv',sep='\t').symbol.to_numpy()
ex=set(json.loads((O/'resources/MAPPING_EXCLUDED_GENES.json').read_text()))
allowed=np.array([g not in ex and not re.match('^(MT-|RPL|RPS)',g) and g not in {'XIST','TSIX','MALAT1','DDX3Y','KDM5D','UTY','EIF1AY','RPS4Y1','RPS4Y2','ZFY','USP9Y'} for g in genes])
ref=m.cohort.eq('baseline').to_numpy()&~m.doublet.to_numpy(bool);added=m.cohort.eq('added').to_numpy()
split={}
for src,z in m[ref].groupby('source'):
 ds=sorted(z.donor_key.unique(),key=lambda x:hashlib.sha256(('20261023|'+x).encode()).hexdigest())
 if len(ds)<3:
  for d in ds:split[d]='train'
 else:
  n=max(1,round(len(ds)*.2));nt=len(ds)-2*n
  for i,d in enumerate(ds):split[d]='train' if i<nt else 'calibration' if i<nt+n else 'test'
m['split']=m.donor_key.map(split).fillna('query');save(m[['source','donor','donor_key','split']].drop_duplicates(),'donor_partitions.tsv')
def aggregate(mask,label):
 ix=np.flatnonzero(mask);meta=m.iloc[ix].reset_index(drop=True);group=meta.groupby(['source','donor_key',label],sort=True).indices
 out=[];inf=[]
 for (src,d,l),idx in group.items():
  out.append(np.asarray(X[ix[idx]].mean(0)).ravel());inf.append(dict(source=src,donor=d,label=l,cells=len(idx)))
 return np.asarray(out),pd.DataFrame(inf)
def unit(a):
 a=a-a.mean(1,keepdims=True);return a/np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-12)
def fit(mask,label):
 data,info=aggregate(mask,label);labels=sorted(info.label.unique());t=[]
 for l in labels:
  t.append(np.mean([data[(info.label.eq(l)&info.source.eq(s)).to_numpy()].mean(0) for s in info.loc[info.label.eq(l),'source'].unique()],axis=0))
 t=np.asarray(t);chosen=set()
 for i in range(len(labels)):
  diff=t[i]-(t.sum(0)-t[i])/max(1,len(labels)-1);ix=np.flatnonzero(allowed&(diff>0));chosen.update(ix[np.argsort(diff[ix])[-(80 if label=='class_label' else 40):]])
 ix=np.array(sorted(chosen));assert len(ix)>30
 mu=data[:,ix].mean(0);sd=np.maximum(data[:,ix].std(0),.1)
 return dict(genes=ix,mu=mu,sd=sd,centroids=unit((t[:,ix]-mu)/sd),profiles=unit((data[:,ix]-mu)/sd),info=info,labels=np.array(labels))
def pred(model,rows):
 a=X[rows][:,model['genes']].toarray();z=unit((a-model['mu'])/model['sd']);sim=z@model['centroids'].T;rank=np.argsort(-sim,axis=1);best=rank[:,0]
 ds=z@model['profiles'].T;info=model['info'];donors=info.donor.to_numpy();labs=info.label.to_numpy();votes=[];ties=[]
 for vals in ds:
  seen=set();picked=[]
  for j in np.argsort(-vals):
   if donors[j] in seen:continue
   seen.add(donors[j]);picked.append(j)
   if len(picked)==min(3,len(set(donors))):break
  count=pd.Series(labs[picked]).value_counts();w=count[count==count.max()].index.tolist();votes.append(w[0]);ties.append(len(w)>1)
 return pd.DataFrame(dict(nearest=model['labels'][best],similarity=sim[np.arange(len(rows)),best],margin=sim[np.arange(len(rows)),best]-sim[np.arange(len(rows)),rank[:,1]],neighbor=votes,tie=ties))
def accept(p,cut):return p.nearest.eq(p.neighbor)&~p.tie&p.similarity.ge(cut['similarity'])&p.margin.ge(cut['margin'])&p.nearest.isin(cut['calibrated_labels'])
results=m.copy();validation=[];diagnostics=[];unknown=[];calibrations={}
for level,col in [('class','class_label'),('state','state_label')]:
 train=ref&m.split.eq('train').to_numpy();calrows=np.flatnonzero(ref&m.split.eq('calibration').to_numpy());testrows=np.flatnonzero(ref&m.split.eq('test').to_numpy())
 model=fit(train,col);cp=pred(model,calrows);truth=m.iloc[calrows][col].to_numpy();correct=cp.nearest.eq(truth)&cp.neighbor.eq(truth)&~cp.tie
 assert correct.any()
 cut=dict(similarity=float(cp.loc[correct,'similarity'].quantile(.05)),margin=float(cp.loc[correct,'margin'].quantile(.05)),calibrated_labels=sorted(set(truth[correct])),calibration_cells=len(cp),correct_dual_calibration_cells=int(correct.sum()),mapping_genes=len(model['genes']))
 calibrations[level]=cut;save(pd.DataFrame({'gene':genes[model['genes']]}),level+'_mapping_genes.tsv');assert not set(genes[model['genes']])&ex
 cp['accepted']=accept(cp,cut);cp['true_label']=truth;cp['cell_id']=m.iloc[calrows].cell_id.to_numpy();save(cp,level+'_calibration_predictions.tsv.gz')
 tp=pred(model,testrows);tp['accepted']=accept(tp,cut);tp['true_label']=m.iloc[testrows][col].to_numpy();tp['source']=m.iloc[testrows].source.to_numpy();tp['donor']=m.iloc[testrows].donor.to_numpy();tp['cell_id']=m.iloc[testrows].cell_id.to_numpy();tp['level']=level;tp['correct']=tp.nearest.eq(tp.true_label)
 if level=='state':
  parent=validation[0].set_index('cell_id').loc[tp.cell_id]
  tp['accepted']=tp.accepted.to_numpy()&parent.accepted.to_numpy()&tp.nearest.str.split('_').str[0].eq(parent.nearest.to_numpy()).to_numpy()
 validation.append(tp)
 for label in sorted(m[ref][col].unique()):
  a=tp[tp.accepted&tp.nearest.eq(label)];t=tp[tp.true_label.eq(label)]
  diagnostics.append(dict(level=level,label=label,test_true_cells=len(t),accepted_cells=len(a),accepted_donors=a[['source','donor']].drop_duplicates().shape[0],correct_accepted=int(a.correct.sum()),accepted_agreement=float(a.correct.mean()) if len(a) else np.nan,true_label_recall=float((t.accepted&t.correct).mean()) if len(t) else np.nan,donor_balanced_accepted_agreement=float(a.groupby(['source','donor']).correct.mean().mean()) if len(a) else np.nan))
 qr=np.flatnonzero(added);p=pred(model,qr);ok=accept(p,cut);results[level+'_compatible']=results[col];results.loc[added,level+'_compatible']=np.where(ok,p.nearest,'Unresolved')
 for c in ['nearest','similarity','margin','neighbor','tie']:results.loc[added,level+'_'+c]=p[c].to_numpy()
 full=fit(ref,col);fp=pred(full,qr);results.loc[added,level+'_fullref_nearest']=fp.nearest.to_numpy();results.loc[added,level+'_fullref_dual_agree']=(fp.nearest.eq(fp.neighbor)&~fp.tie).to_numpy()
 # Closed-set false acceptance diagnostic: leave each parent class entirely out of training.
 if level=='class':
  for held in sorted(m[ref][col].unique()):
   mod=fit(train&~m[col].eq(held).to_numpy(),col);rr=np.flatnonzero(ref&m.split.eq('test').to_numpy()&m[col].eq(held).to_numpy())
   if not len(rr):continue
   pp=pred(mod,rr);ok=accept(pp,cut)
   unknown.append(dict(omitted_class=held,test_cells=len(rr),false_accepted_cells=int(ok.sum()),false_acceptance_fraction=float(ok.mean()),interpretation='Synthetic omitted reference class; not biological unknown truth'))
 print('MAPPED',level,'accepted',int(ok.sum()) if level!='class' else int(accept(p,cut).sum()),flush=True)
# A fine-state assignment is retained only if its parent agrees with independently fitted class assignment.
par=results.state_compatible.str.split('_').str[0];bad=added&par.ne(results.class_compatible).to_numpy();results.loc[bad,'state_compatible']='Unresolved'
save(results,'cell_mapping.tsv.gz');save(pd.concat(validation),'heldout_test_predictions.tsv.gz');save(pd.DataFrame(diagnostics),'heldout_label_diagnostics.tsv');save(pd.DataFrame(unknown),'unknown_class_diagnostic.tsv');dump('resources/CALIBRATION.json',calibrations)
counts=results[added].groupby(['source','class_compatible','state_compatible']).agg(cells=('cell_id','size'),source_patients=('donor','nunique'),doublets=('doublet','sum')).reset_index();save(counts,'incremental_mapping_counts.tsv')
labelnames=pd.read_csv(R/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',sep='\t');save(labelnames[['subgroup_id','working_name','expanded_class_id','class_name']].drop_duplicates(),'frozen_label_dictionary.tsv')
dump('MAPPING_COMPLETE.json',dict(status='PASS_EXPLORATORY',cells=len(results),query_cells=int(added.sum()),heldout_test_cells=int((ref&m.split.eq('test').to_numpy()).sum()),mapping_function_genes_disjoint=True,baseline_labels_unchanged=bool(results.loc[~added,'class_compatible'].equals(results.loc[~added,'class_label'])),state_parent_mismatch_rejected=int(bad.sum()),query_class_counts=results[added].groupby('source').class_compatible.value_counts().to_dict().__str__(),reliable_identity_claim=False))
print(counts.to_string(index=False))
