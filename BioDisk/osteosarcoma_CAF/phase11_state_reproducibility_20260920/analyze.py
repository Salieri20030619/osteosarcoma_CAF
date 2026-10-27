"""Read-only analysis of Phase10 inputs; writes a new Phase11 directory.
SYNTAX_CHECK: python -m py_compile analyze.py
RUN_COMMAND: python -u analyze.py
EXPECTED_SUCCESS_MARKER: ANALYSIS_COMPLETE and exit_code=0
"""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='4'
import json, hashlib, datetime
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
B=P/'phase10_batch_gates_20260919'
O=P/'phase11_state_reproducibility_20260920'
EPS=1e-8  # Numeric-zero tolerance on log1p(CP10k) differences, not a biological effect threshold.
for sub in ['tables','resources','logs','delivery','figures']:
    (O/sub).mkdir(parents=True,exist_ok=True)
def save(df,name): df.to_csv(O/'tables'/name,sep='\t',index=False)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
labels=pd.read_csv(B/'tables/unified_primary_labels.tsv.gz',sep='\t',index_col='cell_id')
labels['cluster']=labels.unified_primary.astype(str)
assert len(labels)==3976 and labels.index.is_unique
programs=json.loads((B/'resources/programs.json').read_text())
programs={k:v for k,v in programs.items() if not k.endswith('_marker')}
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'input':'Phase10 frozen strict 3976 cells; no labels or gates changed',
 'programs':'pre-existing panels; exclude single-marker panel; no SR1004 work',
 'normalization':'gene log1p(CP10000) using full source library denominator, separately raw and DecontX',
 'program_score':'unweighted mean of measured panel gene log1p(CP10000); coverage explicit; >=0.8 for adequately covered program summaries',
 'effect':'within same donor record cluster minus other strict CAF cells; mean cell program score difference',
 'gates':'primary >=20 cells on each side; >=10 sensitivity only; no biological inference tests',
 'replication':'study-stratified summaries; pooled values are donor-record descriptions, not independent N; GSE250015 excluded from pooled primary descriptions',
 'marker_followup':'top30 existing cluster markers, descriptive post-selection patient consistency only',
 'transition':'existing refits only, common strict cells; best Jaccard mapping, purity and recall, full cross-tab retained',
 'refit_program_check':'>=20 remaining original cluster cells, then same donor >=20 cells in/out matched refit cluster; pair with original-label effect on same cell universe',
 'numeric_zero_tolerance':EPS,
 'claim':'exploratory internal reproducibility; no independent external subtype validation'}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
input_hashes={str(B/'objects/unified_strict_integrated.h5ad'):sha(B/'objects/unified_strict_integrated.h5ad')}
markers=pd.read_csv(B/'tables/unified_top30_markers.tsv',sep='\t')
marker_genes=sorted(set(markers.symbol.astype(str)))
genes=sorted(set(sum(programs.values(),[]))|set(marker_genes))
raw=np.full((len(labels),len(genes)),np.nan)
cor=raw.copy(); coverage=[]; seen=np.zeros(len(labels),int)
manifest=pd.read_csv(B/'resources/object_manifest.tsv',sep='\t')
for row in manifest.to_dict('records'):
    path=Path(row['path']); b=ad.read_h5ad(path)
    ids=labels.index.intersection(b.obs_names)
    if len(ids)==0:continue
    input_hashes[str(path)]=sha(path)
    ii=labels.index.get_indexer(ids);bi=b.obs_names.get_indexer(ids)
    sy=b.var.symbol.astype(str).to_numpy(); measured=set(sy)
    for layer,out,den_col in [('counts',raw,'full_library_counts'),('decontx_counts',cor,'corrected_library_counts')]:
        x=b.layers[layer][bi].tocsr();den=b.obs.iloc[bi][den_col].to_numpy(float)
        assert np.isfinite(den).all() and (den>0).all()
        for j,g in enumerate(genes):
            ix=np.flatnonzero(sy==g)
            if len(ix):
                counts=np.asarray(x[:,ix].sum(axis=1)).ravel()
                out[ii,j]=np.log1p(1e4*counts/den)
    for prog,gs in programs.items():
        mm=[g for g in gs if g in measured]
        coverage.append(dict(sample=row['sample'],dataset=row['dataset'],program=prog,measured=len(mm),requested=len(gs),coverage=len(mm)/len(gs),missing=';'.join(g for g in gs if g not in measured)))
    seen[ii]+=1
    print('SOURCE',row['sample'],len(ids),flush=True)
assert (seen==1).all()
save(pd.DataFrame(coverage),'program_measurement_coverage.tsv')
np.savez_compressed(O/'resources/target_gene_expression.npz',raw=raw,corrected=cor,genes=np.array(genes),cell_ids=labels.index.to_numpy(str))
score={};cover={}
for kind,mat in [('raw',raw),('corrected',cor)]:
    score[kind]=np.column_stack([np.nanmean(mat[:,[genes.index(g) for g in gs]],axis=1) for gs in programs.values()])
for name,gs in programs.items():cover[name]=np.isfinite(raw[:,[genes.index(g) for g in gs]]).sum(1)/len(gs)
pn=list(programs)
cellout=labels.copy()
for j,name in enumerate(pn):
    cellout['raw_'+name]=score['raw'][:,j];cellout['corrected_'+name]=score['corrected'][:,j];cellout['coverage_'+name]=cover[name]
save(cellout.reset_index(),'cell_program_scores.tsv.gz')
effects=[];gene_effects=[];eligibility=[]
donors=labels.donor_record.unique();cl=labels.cluster.to_numpy()
for donor in donors:
    dm=labels.donor_record.eq(donor).to_numpy();ds=str(labels.loc[dm,'dataset'].iloc[0])
    for c in sorted(set(cl),key=int):
        inside=dm&(cl==c);outside=dm&(cl!=c);ni=int(inside.sum());no=int(outside.sum())
        eligibility.append(dict(donor_record=donor,dataset=ds,cluster=c,in_cells=ni,out_cells=no,primary20=ni>=20 and no>=20,sensitivity10=ni>=10 and no>=10))
        if min(ni,no)<10:continue
        for j,name in enumerate(pn):
            cv=float(min(cover[name][dm]));r=score['raw'][inside,j].mean()-score['raw'][outside,j].mean();v=score['corrected'][inside,j].mean()-score['corrected'][outside,j].mean()
            effects.append(dict(cluster=c,program=name,donor_record=donor,dataset=ds,in_cells=ni,out_cells=no,primary20=min(ni,no)>=20,coverage=cv,raw_delta=r,corrected_delta=v,direction_same=bool((r>EPS and v>EPS) or (r < -EPS and v < -EPS))))
        wanted=set(markers.loc[markers.cluster.astype(str)==c,'symbol'])|set(sum(programs.values(),[]))
        for gene in sorted(wanted):
            j=genes.index(gene)
            if not np.isfinite(raw[dm,j]).all():continue
            rr=raw[inside,j].mean()-raw[outside,j].mean();cc=cor[inside,j].mean()-cor[outside,j].mean()
            gene_effects.append(dict(cluster=c,gene=gene,donor_record=donor,dataset=ds,primary20=min(ni,no)>=20,raw_delta=rr,corrected_delta=cc,raw_detection_in=float((raw[inside,j]>0).mean()),raw_detection_out=float((raw[outside,j]>0).mean()),existing_top30=gene in set(markers.loc[markers.cluster.astype(str)==c,'symbol'])))
e=pd.DataFrame(effects);save(e,'patient_program_effects.tsv');save(pd.DataFrame(gene_effects),'patient_gene_effects.tsv.gz');save(pd.DataFrame(eligibility),'patient_comparison_eligibility.tsv')
summ=[]
for gate in [20,10]:
    sub=e[(e.in_cells>=gate)&(e.out_cells>=gate)&(e.coverage>=.8)]
    for scope in ['all_records_descriptive','exclude_GSE250015']+sorted(labels.dataset.unique()):
        q=sub if scope=='all_records_descriptive' else sub[sub.dataset.ne('GSE250015')] if scope=='exclude_GSE250015' else sub[sub.dataset.eq(scope)]
        for (c,p),g in q.groupby(['cluster','program']):
            summ.append(dict(scope=scope,min_cells_each_side=gate,cluster=c,program=p,donor_records=len(g),studies=g.dataset.nunique(),positive_raw=int((g.raw_delta>EPS).sum()),positive_corrected=int((g.corrected_delta>EPS).sum()),same_direction=int(g.direction_same.sum()),median_raw_delta=g.raw_delta.median(),min_raw_delta=g.raw_delta.min(),max_raw_delta=g.raw_delta.max(),median_corrected_delta=g.corrected_delta.median()))
save(pd.DataFrame(summ),'program_consistency_summary.tsv')

# Reconstruct the entire transition table, not just the best overlapping cluster.
trans=[];mapping=[];refit=[]
for path in sorted(list((B/'tables').glob('*_labels.tsv.gz'))+list((O/'tables').glob('exclude_*_labels.tsv.gz'))):
    model=path.name.removesuffix('_labels.tsv.gz')
    if model in ['unified_primary','unified_uncorrected']:continue
    t=pd.read_csv(path,sep='\t',index_col='cell_id')
    assert t.index.is_unique
    ids=labels.index.intersection(t.index)
    if not len(ids):continue
    input_hashes[str(path)]=sha(path)
    ii=labels.index.get_indexer(ids);base=labels.loc[ids];new=t.loc[ids,model].astype(str).to_numpy();old=base.cluster.to_numpy()
    for c in sorted(set(cl),key=int):
        cm=old==c;n=int(cm.sum())
        if not n:
            mapping.append(dict(model=model,cluster=c,remaining_primary_cells=0,evaluable=False,best_cluster='',Jaccard=np.nan,recall=np.nan,purity=np.nan,target_cells=0));continue
        matches=[]
        for target in sorted(set(new),key=int):
            tm=new==target;nij=int((cm&tm).sum());nt=int(tm.sum())
            if nij:
                jac=nij/(n+nt-nij);matches.append((jac,target,nij,nt))
                trans.append(dict(model=model,cluster=c,target=target,cells=nij,remaining_primary_cells=n,target_cells=nt,recall=nij/n,purity=nij/nt,Jaccard=jac))
        jac,target,nij,nt=max(matches)
        mapping.append(dict(model=model,cluster=c,remaining_primary_cells=n,evaluable=n>=20,best_cluster=target,Jaccard=jac,recall=nij/n,purity=nij/nt,target_cells=nt))
        if n<20:continue
        for donor in base.donor_record.unique():
            dm=base.donor_record.eq(donor).to_numpy();om=dm&cm;ox=dm&~cm;nm=dm&(new==target);nx=dm&(new!=target)
            if min(om.sum(),ox.sum(),nm.sum(),nx.sum())<20:continue
            ds=str(base.loc[dm,'dataset'].iloc[0])
            for j,p in enumerate(pn):
                cv=float(cover[p][ii[dm]].min())
                for kind in ['raw','corrected']:
                    v=score[kind][ii,j];od=float(v[om].mean()-v[ox].mean());nd=float(v[nm].mean()-v[nx].mean())
                    refit.append(dict(model=model,cluster=c,best_cluster=target,donor_record=donor,dataset=ds,program=p,expression=kind,coverage=cv,original_delta=od,refit_delta=nd,direction_same=(od>EPS and nd>EPS) or (od < -EPS and nd < -EPS),original_direction_evaluable=abs(od)>EPS,original_in=int(om.sum()),refit_in=int(nm.sum())))
    print('TRANSITION',model,flush=True)
save(pd.DataFrame(trans),'all_cluster_transitions.tsv');save(pd.DataFrame(mapping),'best_cluster_transitions.tsv');save(pd.DataFrame(refit),'refit_patient_program_effects.tsv.gz')
rr=pd.DataFrame(refit);summary=[]
for (c,p,kind),g in rr[(rr.coverage>=.8)&rr.dataset.ne('GSE250015')].groupby(['cluster','program','expression']):
    # Repeated models are diagnostic conditions, never biological replicates.
    for donor,d in g.groupby('donor_record'):
        valid=d[d.original_direction_evaluable]
        summary.append(dict(cluster=c,program=p,expression=kind,donor_record=donor,models=len(d),nonzero_baseline_models=len(valid),same_direction_models=int(valid.direction_same.sum()),direction_retention=valid.direction_same.mean(),median_refit_delta=d.refit_delta.median(),min_refit_delta=d.refit_delta.min(),max_refit_delta=d.refit_delta.max()))
save(pd.DataFrame(summary),'refit_program_direction_by_donor.tsv')
assert sha(B/'objects/unified_strict_integrated.h5ad')==input_hashes[str(B/'objects/unified_strict_integrated.h5ad')]
(O/'resources/input_sha256.json').write_text(json.dumps(input_hashes,indent=2))
status={'status':'ANALYSIS_COMPLETE_REVIEW_PENDING','cells':len(labels),'clusters':len(set(cl)),'donor_records':len(donors),'programs':len(programs),'patient_program_rows':len(e),'refit_effect_rows':len(rr),'original_object_unchanged':True,'inferential_tests':False,'independent_validation':False}
(O/'STATUS.json').write_text(json.dumps(status,indent=2));(O/'ANALYSIS_COMPLETE').write_text('ANALYSIS_COMPLETE\n')
print(json.dumps(status),flush=True)
