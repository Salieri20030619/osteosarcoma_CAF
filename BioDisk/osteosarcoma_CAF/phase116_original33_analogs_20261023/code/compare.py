from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd,json,hashlib,itertools,datetime
from scipy import sparse
from scipy.stats import rankdata,false_discovery_control
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=R/'phase115_incremental_function_20261023';O=R/'phase116_original33_analogs_20261023'
for d in ['code','resources','tables','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def dump(n,v):(O/n).write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
paths=[R/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv',R/'phase62_F1F5_internal_20261005/tables/marker_raw_corrected_evidence.tsv',R/'phase62_F1F5_internal_20261005/tables/pathway_display.tsv',R/'phase62_F1F5_internal_20261005/tables/four_axis_internal_overlay.tsv',B/'objects/log1p_cp10k_shared.npz',B/'tables/cell_registry.tsv.gz',B/'tables/shared_genes.tsv',B/'tables/cell_function_scores.tsv.gz',B/'tables/query_cluster_membership.tsv.gz',B/'tables/query_cluster_stability.tsv',B/'tables/functional_statistics.tsv',B/'tables/NCC_cluster_features.tsv',B/'tables/NCI_cluster_features.tsv',B/'tables/program_coverage.tsv']
paths.append(R/'phase62_F1F5_internal_20261005/resources/libraries_canonical.json')
hashes=[dict(path=str(p),sha256=sha(p)) for p in paths];dump('resources/INPUT_HASHES.json',hashes)
old=read(paths[0]).fillna('');assert len(old)==33;save(old,'original_33_subgroup_dictionary.tsv')
for f in paths[1:4]:save(read(f),'inherited_'+f.name)
m=read(B/'tables/cell_registry.tsv.gz');X=sparse.load_npz(B/'objects/log1p_cp10k_shared.npz');genes=read(B/'tables/shared_genes.tsv').symbol.tolist();gi={g:i for i,g in enumerate(genes)};scores=read(B/'tables/cell_function_scores.tsv.gz');assert m.cell_id.equals(scores.cell_id)
cl=read(B/'tables/query_cluster_membership.tsv.gz').set_index('cell_id');m['query_group']=m.cell_id.map(cl.cluster_compatible).fillna('reference')
signatures={};sigrows=[]
for _,r in old.iterrows():
 requested=[x for x in str(r.stable_markers).split(';') if x];mode='original_stable_markers'
 if not requested:requested=[x for x in str(r.display_markers).split(';') if x];mode='original_display_fallback'
 measured=[g for g in requested if g in gi];signatures[r.subgroup_id]=measured
 sigrows.append(dict(old_state=r.subgroup_id,working_name=r.working_name,signature_origin=mode,requested_genes=len(requested),measured_genes=len(measured),coverage=len(measured)/max(len(requested),1),genes=';'.join(measured),missing=';'.join(set(requested)-set(measured))))
siginfo=pd.DataFrame(sigrows).set_index('old_state');save(siginfo.reset_index(),'original_marker_coverage.tsv')
union=sorted(set().union(*map(set,signatures.values())));G=X[:,[gi[g] for g in union]].toarray();gix={g:i for i,g in enumerate(union)}
cov=read(B/'tables/program_coverage.tsv');programs=cov.loc[cov.coverage.ge(.8),'program'].tolist();F=scores[programs].to_numpy();oldstates=old.subgroup_id.tolist()
libraries=json.loads(paths[-1].read_text());archived=read(paths[2]);archived=archived[archived.scope.eq('all')].copy();terms=archived[['library','Term']].drop_duplicates().sort_values(['library','Term']);pathcov=[];pathvals=[];pathnames=[]
for _,t in terms.iterrows():
 gs=libraries[t.library][t.Term];measured=[g for g in gs if g in gi];ratio=len(measured)/len(gs);pathcov.append(dict(library=t.library,Term=t.Term,defined=len(gs),measured=len(measured),coverage=ratio))
 if ratio>=.8 and len(measured)>=3:
  pathnames.append(t.library+'::'+t.Term);pathvals.append(np.asarray(X[:,[gi[g] for g in measured]].mean(1)).ravel())
PF=np.column_stack(pathvals);F=np.column_stack([F,PF]);programs+=pathnames;save(pd.DataFrame(pathcov),'inherited_pathway_coverage.tsv');save(pd.concat([m[['cell_id','source','donor','query_group']],pd.DataFrame(PF,columns=pathnames)],axis=1),'inherited_pathway_cell_scores.tsv.gz')
fit=m.cohort.eq('baseline')&~m.doublet
def profiles(values,indices,label):
 z=m.loc[indices,['source','donor','site',label]].copy();z.index=np.arange(len(z));z=pd.concat([z,pd.DataFrame(values[indices],columns=[str(i) for i in range(values.shape[1])])],axis=1)
 feats=[str(i) for i in range(values.shape[1])]
 return z.groupby(['source','donor','site',label],dropna=False)[feats].mean().groupby(['source','donor',label]).mean().reset_index()
rp=profiles(G,np.flatnonzero(fit),'state_label');fp=profiles(F,np.flatnonzero(fit),'state_label')
def template(p,label):return p.groupby(['source',label]).mean(numeric_only=True).groupby(label).mean().reindex(oldstates).to_numpy()
T=template(rp,'state_label');TF=template(fp,'state_label');assert np.isfinite(T).all() and np.isfinite(TF).all()
def unit(a):
 a=a-a.mean(1,keepdims=True);return a/np.maximum(np.linalg.norm(a,axis=1,keepdims=True),1e-12)
mu=T.mean(0);sd=np.maximum(T.std(0),.1);tm=unit((T-mu)/sd)
fm=TF.mean(0);fs=np.maximum(TF.std(0),.05);tf=unit((TF-fm)/fs)
query=m.cohort.eq('added')&~m.doublet
qp=profiles(G,np.flatnonzero(query),'query_group');qfp=profiles(F,np.flatnonzero(query),'query_group')
Q=qp.groupby(['source','query_group']).mean(numeric_only=True);QF=qfp.groupby(['source','query_group']).mean(numeric_only=True);assert Q.index.equals(QF.index)
qsim=unit((Q.to_numpy()-mu)/sd)@tm.T;fsim=unit((QF.to_numpy()-fm)/fs)@tf.T
sr=unit(rankdata(Q.to_numpy(),axis=1))@unit(rankdata(T,axis=1)).T
save(pd.DataFrame(T,index=oldstates,columns=union).reset_index(names='old_state'),'old_marker_templates.tsv');save(pd.DataFrame(TF,index=oldstates,columns=programs).reset_index(names='old_state'),'old_function_templates.tsv')
# No parent-class prerequisite and no relabeling. Every query cell retains 33-way continuous resemblance.
allq=np.flatnonzero(m.cohort.eq('added'));cs=unit((G[allq]-mu)/sd)@tm.T;cr=np.argsort(-cs,axis=1)
cell=m.iloc[allq][['cell_id','source','donor','query_group','doublet']].copy();cell['nearest_original_state']=np.array(oldstates)[cr[:,0]];cell['similarity']=cs[np.arange(len(cs)),cr[:,0]];cell['runner_up_original_state']=np.array(oldstates)[cr[:,1]];cell['margin']=cs[np.arange(len(cs)),cr[:,0]]-cs[np.arange(len(cs)),cr[:,1]];save(cell,'cell_original33_retrieval.tsv.gz')
# Exact signature differences, conditional on existing query groups. Unused-feature arm is recorded separately.
tests=[];effects=[];details=[];signcache={}
def test(e):
 n=len(e)
 if n<2 or np.all(np.abs(e)<1e-12):return np.nan
 if n not in signcache:signcache[n]=np.array(list(itertools.product((-1.,1.),repeat=n)))
 return float(np.mean(np.abs(signcache[n]@e)>=abs(np.sum(e))-1e-10))
for src in ['NCC','NCI']:
 hvg=set(read(B/'tables'/(src+'_cluster_features.tsv')).symbol)
 base=m[m.source.eq(src)&~m.doublet]
 for group in sorted(base.query_group.unique()):
  for state,gs in signatures.items():
   for arm in ['all_original_markers','unused_in_query_clustering']:
    selected=gs if arm=='all_original_markers' else [g for g in gs if g not in hvg]
    for minimum in [1,5]:
     ee=[]
     for donor,d in base.groupby('donor'):
      a=d.index[d.query_group.eq(group)].to_numpy();b=d.index[d.query_group.ne(group)].to_numpy()
      if min(len(a),len(b))<minimum or not selected:continue
      dv=G[a][:,[gix[g] for g in selected]].mean(0)-G[b][:,[gix[g] for g in selected]].mean(0);delta=float(dv.mean());ee.append(delta)
      effects.append(dict(source=src,query_group=group,old_state=state,arm=arm,minimum=minimum,donor=donor,effect=delta,target_cells=len(a),rest_cells=len(b)))
      if arm=='all_original_markers' and minimum==1:
       for g,v in zip(selected,dv):details.append(dict(source=src,query_group=group,old_state=state,donor=donor,gene=g,effect=v))
     ee=np.array(ee);valid=len(selected)>=3 and siginfo.loc[state,'coverage']>=.8
     tests.append(dict(source=src,query_group=group,old_state=state,arm=arm,minimum=minimum,genes=len(selected),N=len(ee),mean_effect=float(ee.mean()) if len(ee) else np.nan,positive_patients=int((ee>1e-12).sum()),P=test(ee) if valid else np.nan,reason='low_gene_coverage' if not valid else 'insufficient_patients' if len(ee)<2 else 'all_zero' if np.all(abs(ee)<1e-12) else 'exact_patient_signflip'))
res=pd.DataFrame(tests);res['FDR']=np.nan
for _,idx in res.groupby(['source','arm','minimum']).groups.items():
 ii=res.loc[idx].index[res.loc[idx].P.notna()];res.loc[ii,'FDR']=false_discovery_control(res.loc[ii,'P'].to_numpy())
save(res,'all33_signature_tests.tsv');save(pd.DataFrame(effects),'signature_patient_effects.tsv.gz');gd=pd.DataFrame(details);save(gd,'original_marker_patient_effects.tsv.gz')
# Reevaluate the archived pathway portraits using the unchanged historical gene memberships.
pathstats=[]
for src in ['NCC','NCI']:
 base=m[m.source.eq(src)&~m.doublet]
 for group in sorted(base.query_group.unique()):
  for minimum in [1,5]:
   pe=[]
   for donor,d in base.groupby('donor'):
    a=d.index[d.query_group.eq(group)].to_numpy();b=d.index[d.query_group.ne(group)].to_numpy()
    if min(len(a),len(b))>=minimum:pe.append(PF[a].mean(0)-PF[b].mean(0))
   pe=np.array(pe) if pe else np.empty((0,len(pathnames)))
   for j,name in enumerate(pathnames):pathstats.append(dict(source=src,query_group=group,minimum=minimum,pathway=name,N=len(pe),mean_effect=float(pe[:,j].mean()) if len(pe) else np.nan,positive_patients=int((pe[:,j]>1e-12).sum()) if len(pe) else 0,P=test(pe[:,j])))
ps=pd.DataFrame(pathstats);ps['FDR']=np.nan
for _,idx in ps.groupby(['source','minimum']).groups.items():
 ii=ps.loc[idx].index[ps.loc[idx].P.notna()];ps.loc[ii,'FDR']=false_discovery_control(ps.loc[ii,'P'].to_numpy())
save(ps,'inherited_pathway_newgroup_statistics.tsv')
ge=gd.groupby(['source','query_group','old_state','gene']).agg(mean_effect=('effect','mean'),N=('donor','nunique'),positive_patients=('effect',lambda x:int((x>1e-12).sum()))).reset_index();save(ge,'original_marker_gene_evidence.tsv')
stability=read(B/'tables/query_cluster_stability.tsv').set_index('cluster');rows=[]
for i,(src,group) in enumerate(Q.index):
 ranks=rankdata(-qsim[i],method='min');franks=rankdata(-fsim[i],method='min');rranks=rankdata(-sr[i],method='min')
 for j,state in enumerate(oldstates):
  ev=ge[ge.source.eq(src)&ge.query_group.eq(group)&ge.old_state.eq(state)];z=res[(res.source==src)&(res.query_group==group)&(res.old_state==state)&(res.arm=='all_original_markers')&(res.minimum==1)].iloc[0]
  hu=res[(res.source==src)&(res.query_group==group)&(res.old_state==state)&(res.arm=='unused_in_query_clustering')&(res.minimum==1)].iloc[0]
  mi=res[(res.source==src)&(res.query_group==group)&(res.old_state==state)&(res.arm=='all_original_markers')&(res.minimum==5)].iloc[0]
  pos=ev[ev.mean_effect.gt(0)].gene.tolist();cons=ev[(ev.mean_effect.gt(0))&(ev.positive_patients>=np.ceil(ev.N*.6))].gene.tolist()
  rows.append(dict(source=src,query_group=group,old_state=state,working_name=old.set_index('subgroup_id').loc[state,'working_name'],marker_similarity=qsim[i,j],marker_rank_within_query=int(ranks[j]),raw_marker_spearman=sr[i,j],raw_marker_rank=int(rranks[j]),function_similarity=fsim[i,j],function_rank_within_query=int(franks[j]),signature_genes=len(signatures[state]),positive_marker_genes=len(pos),consistent_positive_marker_genes=len(cons),consistent_markers=';'.join(cons),N=int(z.N),signature_effect=z.mean_effect,P=z.P,FDR=z.FDR,unused_genes=int(hu.genes),unused_effect=hu.mean_effect,unused_P=hu.P,unused_FDR=hu.FDR,min5_N=int(mi.N),min5_effect=mi.mean_effect,min5_P=mi.P,min5_FDR=mi.FDR,cells=int(stability.loc[group,'singlet_cells']),group_patients=int(stability.loc[group,'source_patients']),dominant_patient_fraction=stability.loc[group,'dominant_patient_fraction'],seed_jaccard=stability.loc[group,'seed_best_jaccard']))
out=pd.DataFrame(rows)
# Rank each old state across new groups separately by source; do not force a one-to-one match.
out['marker_rank_for_old_state']=out.groupby(['source','old_state']).marker_similarity.rank(ascending=False,method='min').astype(int)
out['function_rank_for_old_state']=out.groupby(['source','old_state']).function_similarity.rank(ascending=False,method='min').astype(int)
out['retrieval_tier']=np.select([(out.marker_rank_within_query<=3)&(out.function_rank_within_query<=3)&(out.signature_effect>0)&(out.consistent_positive_marker_genes>=3),(out.marker_rank_within_query<=3)&(out.signature_effect>0),(out.function_rank_within_query<=3)],['marker_and_function_candidate','marker_candidate','function_only_candidate'],default='lower_rank_or_discordant')
save(out,'all_33_by_19_comparisons.tsv');top=out.sort_values(['source','old_state','marker_rank_for_old_state']).groupby(['source','old_state']).head(3);save(top,'per_original_subgroup_top3.tsv');save(out[out.retrieval_tier.eq('marker_and_function_candidate')],'concordant_retrieval_candidates.tsv')
# Retrieval diagnostic: leave a reference donor out of all33 templates and retrieve their subgroup profiles.
cv=[]
for donor in rp.donor.unique():
 train=rp[rp.donor.ne(donor)];testp=rp[rp.donor.eq(donor)];tt=train.groupby(['source','state_label']).mean(numeric_only=True).groupby('state_label').mean();labels=tt.index.to_numpy();a=tt.to_numpy();mm=a.mean(0);ss=np.maximum(a.std(0),.1);vv=unit((testp[[str(i) for i in range(len(union))]].to_numpy()-mm)/ss)@unit((a-mm)/ss).T;rank=np.argsort(-vv,axis=1)
 for k,(_,z) in enumerate(testp.iterrows()):cv.append(dict(donor=donor,source=z.source,true_state=z.state_label,nearest_state=labels[rank[k,0]],top3=';'.join(labels[rank[k,:3]]),top1_correct=labels[rank[k,0]]==z.state_label,top3_correct=z.state_label in labels[rank[k,:3]],state_in_training=z.state_label in labels))
cv=pd.DataFrame(cv);save(cv,'reference_leave_donor_out_retrieval.tsv');save(cv.groupby('true_state').agg(profiles=('donor','size'),top1_agreement=('top1_correct','mean'),top3_agreement=('top3_correct','mean'),donors=('donor','nunique')).reset_index(),'reference_retrieval_diagnostics.tsv')
dump('resources/PROTOCOL.json',dict(scope='All 33 original subgroups, original Phase62 marker signatures and functional portraits; no F1-F5 parent gate',reference='Frozen8174, excluding26model-positive cells for templates; every old subgroup retained',markers='Original stable markers; display fallback if no stable list; common assayed unique symbols; no newly selected signatures',template='Site then donor then source equal means; standardized correlation across union of original markers; raw Spearman sensitivity',functions=f'22 adequately covered inherited axes plus {len(pathnames)} archived Phase62 displayed pathway gene sets; full gene memberships scored from raw counts, not merely old pathway titles. Similarity is descriptive, not independent validation or a rerun of GSEA.',query='Reuse Phase115 within-source groups as retrieval containers only; all2336 query cells also get continuous33-way ranking; no new cell identity assigned',tests='Patient-level within new-source group vs rest signature mean difference, exact signflip, source/arm/minimum BH families across all33 signatures and query groups; pathway family source/minimum across all inherited pathways and query groups separate. Original signatures were within-parent contrasts; comparison universe differs in queries.',circularity='Full signature and inherited pathways can overlap query-clustering features; separate unused-marker-feature arm records residual support. Reference retrieval uses historically selected markers, not independent external validation.',ranking='All627 old-state/new-group pairs retained; each old state top3 per source; no forced one-to-one mapping or global similarity cutoff',candidate_rule='Navigation flag only: old state top3 within query in marker and functional similarity, signature effect positive, >=3 markers positive in >=60% paired donors. Rank4 or nonnominal candidates remain in full results; not a significance or identity gate.',significance='RawP<0.05 nominal; FDR separately; N<10 exploratory; no untested P values',figures='remain paused',original_results_modified=False))
assert len(out)==33*19 and top.groupby(['source','old_state']).size().eq(3).all();assert len(cell)==2336
for h in hashes:assert sha(h['path'])==h['sha256']
dump('COMPUTATION_COMPLETE.json',dict(status='PASS_RETRIEVAL_SCOPE',original_subgroups=33,new_groups=19,comparisons=len(out),query_cells=len(cell),marker_genes=len(union),sources=2,parent_class_gate=False,concordant_pairs=int(out.retrieval_tier.eq('marker_and_function_candidate').sum()),input_hashes_unchanged=len(hashes),labels_modified=False))
print(out[out.retrieval_tier.eq('marker_and_function_candidate')][['source','query_group','old_state','marker_similarity','function_similarity','consistent_positive_marker_genes','N','P','FDR','min5_N','seed_jaccard']].to_string(index=False));print('F2C3_TOP');print(top[top.old_state.eq('F2_C3')][['source','query_group','marker_similarity','marker_rank_within_query','function_similarity','function_rank_within_query','consistent_markers','N','signature_effect','P','FDR','min5_N','seed_jaccard']].to_string(index=False))
