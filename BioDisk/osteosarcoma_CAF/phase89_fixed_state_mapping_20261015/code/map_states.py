from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='2'
import json,hashlib,datetime,gc,re
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
A=P/'phase41_full_transcriptome_setup_20261001'
Q=P/'phase88_expansion_identity_20261014'
O=P/'phase89_fixed_state_mapping_20261015'
for d in ['code','resources','tables','figures','logs']:(O/d).mkdir(parents=True,exist_ok=True)
def dump(name,obj):
 (O/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False,encoding='utf-8')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as g:
  for b in iter(lambda:g.read(8388608),b''):h.update(b)
 return h.hexdigest()
def var_ids(f):
 with h5py.File(f) as h:
  key=h['var'].attrs['_index']
  return pd.Index(h['var'][key].asstr()[:])
def normalize(x,den):
 x=x.astype(np.float32).multiply((10000/den)[:,None]).tocsr()
 x.data=np.log1p(x.data)
 return x
def aggregate(X,meta,cols):
 key=meta[cols].astype(str).agg('|'.join,axis=1)
 uniq=pd.Index(key.unique());codes=uniq.get_indexer(key)
 ct=np.bincount(codes,minlength=len(uniq))
 S=sparse.csr_matrix((1/ct[codes],(codes,np.arange(len(meta)))),shape=(len(uniq),len(meta)))
 vals=S@X
 if sparse.issparse(vals):vals=vals.toarray()
 info=meta.iloc[[np.flatnonzero(codes==i)[0] for i in range(len(uniq))]][cols].reset_index(drop=True)
 info['unit_cells']=ct
 return np.asarray(vals,np.float32),info
meta_path=P/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz'
query_path=P/'audit_expanded_admission_20261015/tables/final_cell_review_decisions.tsv.gz'
tracked=[meta_path,query_path,A/'objects/all_candidate_full_genes.h5ad',P/'CAF_FUNCTIONAL_BONE_LUNG_CURRENT.json']
hashes={str(f):sha(f) for f in tracked}
query=pd.read_csv(query_path,sep='\t',dtype={'cluster':str},low_memory=False)
query=query[query.expanded_review_retained].copy().reset_index(drop=True)
assert len(query)==1268 and query.cell_id.is_unique
ref=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad')
meta=pd.read_csv(meta_path,sep='\t',low_memory=False).set_index('cell_id',drop=False).loc[ref.obs_names].reset_index(drop=True)
assert len(meta)==8174 and meta.subgroup_id.nunique()==33
pres=sparse.load_npz(A/'objects/sample_feature_presence.npz')
pidx=pd.read_csv(A/'tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample')
samples=sorted(meta['sample'].unique());keep=np.asarray(pres[pidx.loc[samples,'presence_row'].to_numpy()].sum(0)).ravel()==len(samples)
shared=set(ref.var_names[keep])
for sm in query['sample'].unique():shared &= set(var_ids(Q/'preprocess/accepted'/sm/'analysis_ready.h5ad'))
shared=sorted(shared);ridx=ref.var_names.get_indexer(shared)
symbols=ref.var.iloc[ridx].symbol.astype(str).to_numpy()
ambig=ref.var.iloc[ridx].symbol_ambiguous.to_numpy(bool)
den=np.asarray(ref.X.sum(1)).ravel()
assert np.allclose(den,meta.full_library_counts)
R=normalize(ref.X[:,ridx],den)
reference_ids=meta.cell_id.to_numpy()
rv=ref.var.iloc[ridx].copy()
del ref;gc.collect()
query_parts=[];query_order=[]
for sm,g in query.groupby('sample',sort=False):
 a=ad.read_h5ad(Q/'preprocess/accepted'/sm/'analysis_ready.h5ad')
 ix=a.obs_names.get_indexer(g.cell_id);assert (ix>=0).all()
 x=a.layers['counts'][ix,:].tocsr()
 den=np.asarray(x.sum(1)).ravel();assert np.allclose(den,g.n_counts)
 qi=a.var_names.get_indexer(shared);assert (qi>=0).all()
 query_parts.append(normalize(x[:,qi],den));query_order.extend(g.index.tolist())
 print('QUERY_READ',sm,len(g),flush=True)
 del a,x;gc.collect()
X=sparse.vstack(query_parts,format='csr')
query=query.loc[query_order].reset_index(drop=True)
assert query.cell_id.is_unique and X.shape==(1268,len(shared))
# Marker selection is reference-only. No clinical group or new-donor independence enters mapping.
technical=np.array([bool(re.match(r'^(MT-|RPL|RPS)',g)) or g in {'XIST','TSIX','MALAT1','DDX3Y','KDM5D','UTY','EIF1AY','RPS4Y1','RPS4Y2','ZFY','USP9Y'} for g in symbols])
allowed=~technical & ~ambig & np.array([bool(g) and not g.startswith('ENSG') for g in symbols])
basekeep=~meta.doublet_positive_any.to_numpy(bool)
D,dm=aggregate(R[basekeep],meta.loc[basekeep].reset_index(drop=True),['dataset','donor_record','subgroup_id','class_id'])
states=sorted(meta.subgroup_id.unique());classes=sorted(meta.class_id.unique())
state_names=meta.groupby('subgroup_id').working_name.first().to_dict()
protocol={'created_utc':stamp(),'reference':'Frozen Phase30/58 classes and33 states; Phase66 metadata; all8174 labels unchanged','reference_training':'Exclude26 explicit model-positive doublets from templates; preserve original baseline and all query risk flags','query':'1268 reviewed expanded exploratory cells, including228 reference-failure cells; no patient-independence or anatomical-site gate','common_features':len(shared),'normalization':'Per-cell log1p(CP10k), denominator is the full original library before feature selection; shared measured Ensembl features only, no missing-feature zero imputation','templates':'Mean cells per source-scoped donor/state; mean donor profiles within source/state; equal source mean for state templates','mapping_markers':'Top40 positive state-versus-other-template differences per state, union, selected on reference only; exclude ambiguous symbols, MT/ribosomal and listed sex/technical genes; no new cell entry thresholds','method1':'Correlation of gene-standardized query profile to source-balanced state templates','method2':'Nearest donor/state profiles, retaining only best profile per source-scoped donor; plurality among5 nearest distinct donor records; tie not accepted','validation':'Leave one original dataset out, recompute features/templates; test heldout donor/state mean profiles; internal transfer diagnostic, not independent external validation or patient-level biological inference','calibration':'Global lower5th percentiles of similarity and margin among correctly state-mapped heldout profiles; state must have at least one correct heldout profile for supported state mapping; uncalibrated states retained as nearest-only','group_decision':'Original new sample/cluster mean profiles; two methods and calibrated score/margin for state support; same-class agreement for class-only; no automatic closed-set assignment','cell_decision':'State-compatible only if both single-cell methods agree with the calibrated group state; class-compatible when the group and both single-cell methods agree at parent class; otherwise uncertain','functions':'Existing Phase63 program definitions scored on common measured genes, descriptive only; functional differences not tested','order':'Complete mapping first, then inspect source patient independence','no_reclustering':True,'no_new_inferential_tests':True}
dump('resources/PROTOCOL.json',protocol)
save(pd.DataFrame({'feature':shared,'symbol':symbols,'ambiguous':ambig,'technical_excluded_from_mapping':technical}),'shared_feature_registry.tsv')
def templates(data,info):
 st=sorted(info.subgroup_id.unique());out=[]
 for s in st:
  source_means=[data[np.flatnonzero((info.subgroup_id.eq(s)&info.dataset.eq(ds)).to_numpy())].mean(0) for ds in info.loc[info.subgroup_id.eq(s),'dataset'].unique()]
  out.append(np.mean(source_means,axis=0))
 return np.asarray(out,np.float32),np.array(st,dtype=object)
def unit(z):
 z=z-z.mean(1,keepdims=True)
 return z/np.maximum(np.linalg.norm(z,axis=1,keepdims=True),1e-12)
def fit(data,info):
 t,st=templates(data,info);select=set()
 for i in range(len(st)):
  contrast=t[i]-(t.sum(0)-t[i])/max(1,len(st)-1)
  valid=np.flatnonzero(allowed & (contrast>0))
  select.update(valid[np.argsort(contrast[valid])[-40:]].tolist())
 ids=np.array(sorted(select),dtype=int)
 mu=data[:,ids].mean(0);sd=np.maximum(data[:,ids].std(0),.1)
 tt=unit((t[:,ids]-mu)/sd)
 dd=unit((data[:,ids]-mu)/sd)
 return {'ids':ids,'mu':mu,'sd':sd,'templates':tt,'states':st,'donor_profiles':dd,'meta':info.reset_index(drop=True)}
def predict(model,data):
 z=unit((data[:,model['ids']]-model['mu'])/model['sd'])
 sim=z@model['templates'].T
 order=np.argsort(-sim,axis=1)
 out=pd.DataFrame({'template_state':model['states'][order[:,0]],'second_state':model['states'][order[:,1]],'similarity':sim[np.arange(len(sim)),order[:,0]],'margin':sim[np.arange(len(sim)),order[:,0]]-sim[np.arange(len(sim)),order[:,1]]})
 ds=z@model['donor_profiles'].T
 mi=model['meta'];donors=mi.donor_record.astype(str).to_numpy();labels=mi.subgroup_id.astype(str).to_numpy()
 votes=[];details=[]
 for i in range(len(data)):
  seen=set();picked=[]
  for j in np.argsort(-ds[i]):
   if donors[j] in seen:continue
   seen.add(donors[j]);picked.append(j)
   if len(picked)==min(5,len(set(donors))):break
  labs=labels[picked];counts=pd.Series(labs).value_counts();mx=int(counts.max());ties=counts[counts==mx].index.tolist()
  winner=max(ties,key=lambda s:float(ds[i,np.array(picked)[labs==s]].mean()))
  parent=pd.Series([s.split('_')[0] for s in labs]).value_counts();pmax=int(parent.max());pties=parent[parent==pmax].index.tolist()
  pwin=max(pties,key=lambda s:float(ds[i,np.array(picked)[np.array([l.split('_')[0] for l in labs])==s]].mean()))
  votes.append({'neighbor_state':winner,'neighbor_state_votes':mx,'neighbor_state_tie':len(ties)>1,'neighbor_class':pwin,'neighbor_class_votes':pmax,'neighbor_class_tie':len(pties)>1,'neighbor_donors':'|'.join(donors[picked]),'neighbor_states':'|'.join(labs)})
 out=pd.concat([out,pd.DataFrame(votes)],axis=1);out['template_class']=out.template_state.str.split('_').str[0]
 return out,sim
valid=[]
for source in sorted(dm.dataset.unique()):
 train=dm.dataset.ne(source).to_numpy();test=~train
 model=fit(D[train],dm.loc[train].reset_index(drop=True))
 pp,_=predict(model,D[test])
 pp['heldout_source']=source;pp['true_state']=dm.loc[test,'subgroup_id'].to_numpy();pp['true_class']=dm.loc[test,'class_id'].to_numpy();pp['donor_record']=dm.loc[test,'donor_record'].to_numpy()
 pp['state_available_in_training']=pp.true_state.isin(model['states'])
 pp['correct_state']=pp.template_state.eq(pp.true_state);pp['correct_class']=pp.template_class.eq(pp.true_class)
 pp['dual_correct_state']=pp.correct_state & pp.neighbor_state.eq(pp.true_state) & ~pp.neighbor_state_tie
 valid.append(pp);print('VALIDATED_SOURCE',source,len(pp),int(pp.correct_state.sum()),flush=True)
cv=pd.concat(valid,ignore_index=True);save(cv,'leave_source_out_predictions.tsv')
correct=cv[cv.correct_state]
assert len(correct)>0
simcut=float(correct.similarity.quantile(.05));margincut=float(correct.margin.quantile(.05))
calstates=set(correct.true_state)
classcorrect=cv[cv.correct_class];classcut=float(classcorrect.similarity.quantile(.05))
valstate=cv.groupby('true_state').agg(validation_profiles=('true_state','size'),available_profiles=('state_available_in_training','sum'),correct_state=('correct_state','sum'),dual_correct_state=('dual_correct_state','sum'),correct_class=('correct_class','sum'),heldout_sources=('heldout_source','nunique')).reset_index().rename(columns={'true_state':'subgroup_id'})
valstate['state_recall']=valstate.correct_state/valstate.validation_profiles;valstate['working_name']=valstate.subgroup_id.map(state_names)
save(valstate,'state_transfer_validation.tsv')
cal={'state_similarity_lower05':simcut,'state_margin_lower05':margincut,'class_similarity_lower05':classcut,'validation_profiles':len(cv),'state_correct':int(cv.correct_state.sum()),'class_correct':int(cv.correct_class.sum()),'calibratable_states':sorted(calstates),'not_independent_external_validation':True}
dump('resources/CALIBRATION.json',cal)
model=fit(D,dm)
save(pd.DataFrame({'feature':np.array(shared)[model['ids']],'symbol':symbols[model['ids']]}),'reference_selected_mapping_features.tsv')
G,gm=aggregate(X,query,['sample','dataset','patient_id','cluster'])
gp,gsim=predict(model,G)
groups=pd.concat([gm,gp],axis=1)
stateok=groups.template_state.eq(groups.neighbor_state)&~groups.neighbor_state_tie & groups.similarity.ge(simcut)&groups.margin.ge(margincut)&groups.template_state.isin(calstates)
classok=groups.template_class.eq(groups.neighbor_class)&~groups.neighbor_class_tie & groups.similarity.ge(classcut)
groups['mapping_evidence']=np.select([stateok,classok],['state_compatible','class_only_compatible'],default='unresolved')
groups['working_name']=groups.template_state.map(state_names)
save(groups,'sample_cluster_state_mapping.tsv')
save(pd.concat([gm,pd.DataFrame(gsim,columns=model['states'])],axis=1),'sample_cluster_all_state_similarities.tsv')
cp,csim=predict(model,X[:,model['ids']].toarray() if False else X.toarray())
cell=pd.concat([query,cp],axis=1)
gmapping=groups[['sample','cluster','template_state','template_class','mapping_evidence']].rename(columns={'template_state':'group_state','template_class':'group_class','mapping_evidence':'group_mapping_evidence'})
cell=cell.merge(gmapping,on=['sample','cluster'],validate='many_to_one')
cellstate=(cell.group_mapping_evidence.eq('state_compatible')&cell.template_state.eq(cell.group_state)&cell.neighbor_state.eq(cell.group_state)&~cell.neighbor_state_tie)
cellclass=(cell.group_mapping_evidence.isin(['state_compatible','class_only_compatible'])&cell.template_class.eq(cell.group_class)&cell.neighbor_class.eq(cell.group_class)&~cell.neighbor_class_tie)
cell['mapping_evidence']=np.select([cellstate,cellclass],['state_compatible','class_only_compatible'],default='unresolved')
cell['mapped_state']=np.where(cellstate,cell.group_state,'')
cell['mapped_class']=np.where(cellstate|cellclass,cell.group_class,'')
cell['working_name']=cell.mapped_state.map(state_names).fillna('')
cell['formal_baseline_admission']=False
save(cell,'cell_state_mapping.tsv.gz')
# Query group function scores use the same available measured genes as the reference.
programs=json.loads((P/'phase63_matrix_resistance_20261005/resources/PROGRAMS.json').read_text(encoding='utf-8'))['sets']
funrows=[];funref=[];cover=[]
for name,genes in programs.items():
 ids=np.flatnonzero(np.isin(symbols,genes)&~ambig)
 cover.append({'program':name,'measured_symbols':len(set(symbols[ids])),'requested':len(set(genes))})
 if len(ids)==0:continue
 qs=G[:,ids].mean(1)
 rt,st=templates(D,dm)
 rs=rt[:,ids].mean(1)
 for i in range(len(gm)):funrows.append({'sample':gm.iloc[i]['sample'],'cluster':gm.iloc[i]['cluster'],'program':name,'mean_logCP10k':float(qs[i])})
 for i in range(len(st)):funref.append({'subgroup_id':st[i],'program':name,'mean_logCP10k':float(rs[i])})
save(pd.DataFrame(funrows),'query_group_function_profiles.tsv');save(pd.DataFrame(funref),'reference_state_function_profiles.tsv');save(pd.DataFrame(cover),'function_feature_coverage.tsv')
# Summaries keep biological identity, mapping evidence, and reference calibration separate.
state=cell[cell.mapping_evidence.eq('state_compatible')].groupby(['mapped_class','mapped_state','working_name'],dropna=False).agg(cells=('cell_id','size'),source_patient_labels=('patient_id','nunique'),samples=('sample','nunique'),reference_pass=('reference_QC_pass','sum')).reset_index()
save(state,'compatible_state_summary.tsv')
save(cell.groupby(['sample','patient_id','mapping_evidence','mapped_class','mapped_state'],dropna=False).size().reset_index(name='cells'),'sample_state_counts.tsv')
save(cell.groupby(['dataset','mapping_evidence']).size().reset_index(name='cells'),'mapping_evidence_by_source.tsv')
# Mapping-only stability check: do not confuse a best match with out-of-reference admission.
baseline_hashes_after={f:sha(f) for f in hashes};assert hashes==baseline_hashes_after
assert len(cell)==1268 and cell.cell_id.is_unique and set(cell.cell_id)==set(query.cell_id)
assert not cell.formal_baseline_admission.any()
assert set(cell.loc[cellstate.to_numpy() if False else cell.mapping_evidence.eq('state_compatible'),'mapped_state'])<=set(states)
dump('resources/INPUT_HASHES.json',{'before':hashes,'after':baseline_hashes_after,'unchanged':True})
plot=cell[cell.mapped_class.ne('')].groupby(['sample','mapped_class']).size().unstack(fill_value=0).reindex(index=query['sample'].unique(),columns=classes,fill_value=0).fillna(0)
fig,ax=plt.subplots(figsize=(11,6.5));bottom=np.zeros(len(plot));colors=['#3b8c81','#d99538','#6678a0','#b96d6d','#7c6e96']
for cl,col in zip(classes,colors):
 vals=plot[cl].to_numpy();ax.barh(np.arange(len(plot)),vals,left=bottom,label=cl,color=col);bottom+=vals
ax.set_yticks(np.arange(len(plot)),plot.index);ax.invert_yaxis();ax.set_xlabel('Cells with group- and cell-level compatible class mapping')
ax.set_title('New expanded candidates projected onto fixed F1-F5')
ax.legend(loc='lower right');fig.tight_layout();fig.savefig(O/'figures/class_mapping.png',dpi=170);plt.close(fig)
counts=cell.mapping_evidence.value_counts().to_dict()
result={'status':'MAPPING_COMPUTATION_COMPLETE_REVIEW_PENDING','completed_utc':stamp(),'reference_cells':8174,'reference_states':33,'reference_templates_exclude_doublets':int((~basekeep).sum()),'query_cells':1268,'query_groups':len(groups),'shared_measured_features':len(shared),'selected_mapping_features':len(model['ids']),'cell_mapping_counts':{k:int(v) for k,v in counts.items()},'group_mapping_counts':{k:int(v) for k,v in groups.mapping_evidence.value_counts().to_dict().items()},'mapped_states':int(cell.mapped_state.replace('',np.nan).nunique()),'validation':cal,'baseline_unchanged':True,'patient_independence_review':'not_started_until_mapping_review','new_inferential_tests':False}
dump('MAPPING_COMPLETE.json',result)
print('MAPPING_COMPLETE',json.dumps(result),flush=True)
