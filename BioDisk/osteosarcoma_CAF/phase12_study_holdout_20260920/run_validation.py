"""SYNTAX_CHECK: python -m py_compile run_validation.py
RUN_COMMAND: python -u run_validation.py
EXPECTED_SUCCESS_MARKER: ANALYSIS_COMPLETE, exit_code=0
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase10_batch_gates_20260919';V=P/'phase11_state_reproducibility_20260920';O=P/'phase12_study_holdout_20260920'
for d in ['tables','resources','logs','delivery']: (O/d).mkdir(parents=True,exist_ok=True)
SEED=20260920;EPS=1e-8
programs=['Stress_core','Cell_cycle','Contractile','Inflammatory','Matrix_remodeling','Antigen_presentation']
protocol={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'exploratory leave-study-out reanalysis of previously seen cohorts, not prospective external validation',
 'excluded_study':'GSE250015 unresolved cross-study overlap','fold':'one of four remaining studies held out in full; no pooled cross-study patient inference',
 'features':'fixed measured common-gene universe, exclude ALL prior 17 program-panel genes, stress512, S43/G2M54, mitochondrial/ribosomal genes; min3 cells and HVG2000 selected only in training',
 'discovery':'train-only scaling, PCA30, Harmony study_technology theta2 dynamic lambda, k15 Leiden0.5',
 'prediction':'C1 L2 logistic trained on uncorrected train PCA coordinates against train-only discovered labels; training-only preprocessing applied to heldout cells',
 'weight':'inverse donor count times inverse class count, normalized mean1',
 'rejection':'primary max_probability>=0.70 and distance to assigned training-class centroid <= class training-distance 95th percentile; descriptive working rules, not calibrated clinical probabilities',
 'sensitivity':'0.50 and 0.90 probability with same distance rule; no outcome-driven threshold choice',
 'evaluation_genes':'not used as predictor features; raw and DecontX donor-internal mean program differences, coverage>=0.8; at least20 accepted cells each side',
 'candidate_rule':'training cluster program raw/corrected both positive in >=75% of >=3 evaluable donor records and at least2 studies; both median effects > numeric epsilon; negative results retained',
 'statistics':'no P values, no cell-level inference; folds and models not biological replicates',
 'programs':programs,'seed':SEED}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=B/'objects/unified_strict_integrated.h5ad';before=sha(source)
a=ad.read_h5ad(source);a=a[a.obs.dataset.ne('GSE250015')].copy();assert a.n_obs==3577
scores=pd.read_csv(V/'tables/cell_program_scores.tsv.gz',sep='\t',index_col='cell_id').loc[a.obs_names]
assert scores.index.equals(a.obs_names)
panels=json.loads((B/'resources/programs.json').read_text());excluded_symbols=set(sum(panels.values(),[]))
R=P/'phase7_review_states_20260917/remediation_20260917/resources'
stress=pd.read_csv(R/'stress_core_author.csv');cycle=pd.read_csv(R/'seurat_cc_updated2019.tsv',sep='\t')
excluded_symbols.update(cycle.gene)
eligible=~a.var.symbol.isin(excluded_symbols).to_numpy()&~a.var_names.isin(stress.ensembl_gene_id)&~a.var.symbol.str.startswith(('MT-','RPL','RPS')).to_numpy()
all_effects=[];train_summary=[];folds=[];assignments=[];candidate_results=[];features=[]
def effect_rows(ids,lab,accepted,fold,split,cutoff):
    ob=a.obs.loc[ids];z=scores.loc[ids];rows=[]
    for donor in ob.donor_record.unique():
        dm=ob.donor_record.eq(donor).to_numpy();ds=str(ob.loc[dm,'dataset'].iloc[0])
        for c in sorted(set(lab)):
            inn=dm&accepted&(lab==c);out=dm&accepted&(lab!=c);ni=int(inn.sum());no=int(out.sum())
            for p in programs:
                cv=float(z.loc[dm,'coverage_'+p].min());ok=min(ni,no)>=20 and cv>=.8
                rows.append(dict(fold=fold,split=split,threshold=cutoff,cluster=str(c),donor_record=donor,dataset=ds,program=p,in_cells=ni,out_cells=no,coverage=cv,evaluable=ok,
                 raw_delta=float(z.loc[inn,'raw_'+p].mean()-z.loc[out,'raw_'+p].mean()) if min(ni,no)>0 else np.nan,
                 corrected_delta=float(z.loc[inn,'corrected_'+p].mean()-z.loc[out,'corrected_'+p].mean()) if min(ni,no)>0 else np.nan))
    return pd.DataFrame(rows)
for study in sorted(a.obs.dataset.unique()):
    train=a[a.obs.dataset.ne(study),eligible].copy();test=a[a.obs.dataset.eq(study)].copy()
    assert not set(train.obs_names)&set(test.obs_names)
    assert not set(train.obs.donor_record)&set(test.obs.donor_record)
    sc.pp.filter_genes(train,min_cells=3)
    sc.pp.highly_variable_genes(train,n_top_genes=min(2000,train.n_vars),flavor='seurat',batch_key='study_technology')
    hv=train.var_names[train.var.highly_variable].tolist();train=train[:,hv].copy()
    assert not train.var.symbol.isin(excluded_symbols).any() and not train.var_names.isin(stress.ensembl_gene_id).any()
    features.extend(dict(fold=study,ensembl_id=g,symbol=str(train.var.loc[g,'symbol'])) for g in hv)
    xx=train.X.toarray();xt=test[:,hv].X.toarray();scaler=StandardScaler().fit(xx)
    xx=np.clip(scaler.transform(xx),-10,10);xt=np.clip(scaler.transform(xt),-10,10)
    pc=PCA(n_components=30,svd_solver='arpack',random_state=SEED).fit(xx);z=pc.transform(xx);zt=pc.transform(xt)
    h=hm.run_harmony(z,train.obs,'study_technology',theta=2,lamb=None,random_state=SEED,ncores=4,max_iter_harmony=30,verbose=False)
    hz=np.asarray(h.Z_corr);hz=hz if hz.shape==z.shape else hz.T
    assert np.isfinite(hz).all();train.obsm['train_harmony']=hz
    sc.pp.neighbors(train,use_rep='train_harmony',n_neighbors=15,random_state=SEED)
    sc.tl.leiden(train,resolution=.5,key_added='discovered',random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
    y=train.obs.discovered.astype(str).to_numpy();don=train.obs.donor_record.astype(str)
    weights=1/don.map(don.value_counts()).to_numpy()/pd.Series(y).map(pd.Series(y).value_counts()).to_numpy();weights/=weights.mean()
    model=LogisticRegression(C=1.0,max_iter=3000,solver='lbfgs',random_state=SEED).fit(z,y,sample_weight=weights)
    assert model.n_iter_.max()<3000
    prob=model.predict_proba(zt);pred=model.classes_[prob.argmax(1)];confidence=prob.max(1)
    distances=np.empty(len(test));limits=np.empty(len(test));class_rows=[]
    for c in model.classes_:
        center=z[y==c].mean(0);r=np.linalg.norm(z[y==c]-center,axis=1);limit=float(np.quantile(r,.95));tm=pred==c
        distances[tm]=np.linalg.norm(zt[tm]-center,axis=1);limits[tm]=limit
        class_rows.append(dict(cluster=c,training_cells=int((y==c).sum()),radius95=limit))
    pd.DataFrame(class_rows).to_json(O/'resources'/(study+'_training_class_radii.json'),orient='records',indent=2)
    train_rows=effect_rows(train.obs_names,y,np.ones(len(train),bool),study,'train',1.0);all_effects.append(train_rows)
    chosen=[]
    for (c,p),g in train_rows[train_rows.evaluable].groupby(['cluster','program']):
        positive=(g.raw_delta>EPS)&(g.corrected_delta>EPS)
        candidate=len(g)>=3 and g.dataset.nunique()>=2 and positive.mean()>=.75 and g.raw_delta.median()>EPS and g.corrected_delta.median()>EPS
        rec=dict(fold=study,cluster=c,program=p,donor_records=len(g),studies=g.dataset.nunique(),both_positive=int(positive.sum()),median_raw=g.raw_delta.median(),median_corrected=g.corrected_delta.median(),candidate=bool(candidate))
        train_summary.append(rec)
        if candidate:chosen.append((c,p))
    for cut in [.5,.7,.9]:
        accepted=(confidence>=cut)&(distances<=limits)
        ef=effect_rows(test.obs_names,pred,accepted,study,'heldout',cut);all_effects.append(ef)
        for c,p in chosen:
            g=ef[(ef.cluster==c)&(ef.program==p)&ef.evaluable]
            candidate_results.append(dict(fold=study,threshold=cut,cluster=c,program=p,heldout_evaluable_records=len(g),positive_raw=int((g.raw_delta>EPS).sum()),positive_corrected=int((g.corrected_delta>EPS).sum()),both_positive=int(((g.raw_delta>EPS)&(g.corrected_delta>EPS)).sum()),median_raw=g.raw_delta.median(),median_corrected=g.corrected_delta.median(),assessment='not_evaluable' if not len(g) else 'all_direction_preserved' if ((g.raw_delta>EPS)&(g.corrected_delta>EPS)).all() else 'direction_not_universal'))
        for donor in test.obs.donor_record.unique():
            dm=test.obs.donor_record.eq(donor).to_numpy()
            assignments.append(dict(fold=study,donor_record=donor,threshold=cut,total=int(dm.sum()),accepted=int((accepted&dm).sum()),rejected=int((~accepted&dm).sum()),low_probability=int(((confidence<cut)&dm).sum()),out_of_distribution=int(((distances>limits)&dm).sum())))
    output=test.obs[['dataset','sample','donor_record','technology','unified_cluster']].copy();output['predicted_training_cluster']=pred;output['probability']=confidence;output['distance']=distances;output['radius95']=limits;output['accepted70']=(confidence>=.7)&(distances<=limits)
    output.to_csv(O/'tables'/(study+'_heldout_assignments.tsv.gz'),sep='\t',index_label='cell_id')
    train.obs[['dataset','donor_record','discovered']].to_csv(O/'tables'/(study+'_training_labels.tsv.gz'),sep='\t',index_label='cell_id')
    rec=dict(heldout_study=study,train_cells=len(train),test_cells=len(test),train_clusters=len(set(y)),training_label_fit_balanced_accuracy=balanced_accuracy_score(y,model.predict(z)),accepted70=int(output.accepted70.sum()),accepted_fraction70=float(output.accepted70.mean()),candidate_program_pairs=len(chosen),excluded_prediction_genes=int((~eligible).sum()))
    folds.append(rec);print('FOLD_COMPLETE',json.dumps(rec),flush=True)
    pd.DataFrame(folds).to_csv(O/'tables/fold_summary.tsv',sep='\t',index=False)
    # Store transform/model coefficients to make predictions reconstructable without pickle.
    np.savez_compressed(O/'resources'/(study+'_model.npz'),genes=np.array(hv),scale_mean=scaler.mean_,scale_scale=scaler.scale_,pca_mean=pc.mean_,pca_components=pc.components_,classes=model.classes_,coef=model.coef_,intercept=model.intercept_)
save(pd.concat(all_effects),'all_patient_program_effects.tsv.gz');save(pd.DataFrame(train_summary),'training_program_candidates.tsv');save(pd.DataFrame(candidate_results),'heldout_candidate_results.tsv');save(pd.DataFrame(assignments),'assignment_attrition.tsv');save(pd.DataFrame(features),'training_selected_features.tsv')
assert sha(source)==before
qa={'source_unchanged':True,'train_test_cells_disjoint':True,'donor_records_disjoint':True,'panel_features_excluded':True,'test_expression_not_used_for_training':True,'logistic_converged':True,'folds':len(folds),'cells':a.n_obs}
(O/'COMPUTATION_QA.json').write_text(json.dumps(qa,indent=2))
(O/'resources/input_sha256.json').write_text(json.dumps({str(p):sha(p) for p in [source,V/'tables/cell_program_scores.tsv.gz',B/'resources/programs.json',R/'stress_core_author.csv',R/'seurat_cc_updated2019.tsv']},indent=2))
(O/'STATUS.json').write_text(json.dumps({'status':'STUDY_HOLDOUT_COMPUTED_REVIEW_PENDING','folds':len(folds),'cells':a.n_obs,'prospective_external_validation':False,'biological_inferential_tests':False},indent=2));(O/'ANALYSIS_COMPLETE').write_text('ANALYSIS_COMPLETE\n')
