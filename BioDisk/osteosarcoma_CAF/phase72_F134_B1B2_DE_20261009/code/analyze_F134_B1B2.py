from pathlib import Path
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '2'
import json, hashlib, datetime, itertools, importlib.util, traceback, shutil, zipfile, warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests

P = Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O = P / 'phase72_F134_B1B2_DE_20261009'
A = P / 'phase41_full_transcriptome_setup_20261001'
CLASSES = ['F1', 'F3', 'F4']

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(f, **kw): return pd.read_csv(f, sep='\t', low_memory=False, **kw)
def dump(f, x): f.write_text(json.dumps(x, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
def save(d, name): d.to_csv(O / 'tables' / name, sep='\t', index=False, na_rep='NA')
def sha(f):
    h = hashlib.sha256()
    with f.open('rb') as s:
        for b in iter(lambda: s.read(8388608), b''): h.update(b)
    return h.hexdigest()
def valid_mean(x):
    n = np.isfinite(x).sum(axis=0)
    out = np.full(x.shape[1], np.nan)
    np.divide(np.nansum(x, axis=0), n, out=out, where=n>0)
    return out
def truth(x): return x.astype(str).str.lower().isin(['true', '1'])

def run():
    assert not (O / 'COMPLETE.json').exists()
    for d in ['code', 'tables', 'objects', 'resources', 'QA', 'delivery']:
        (O / d).mkdir(parents=True, exist_ok=True)
    inputs = {
        'cells': P / 'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',
        'grouping': P / 'phase71_allcell_AB1B2_grouping_20261009/tables/sample_group_manifest.tsv',
        'counts': A / 'objects/all_candidate_full_genes.h5ad',
        'presence': A / 'objects/sample_feature_presence.npz',
        'sample_index': A / 'tables/feature_presence_sample_index.tsv',
        'feature_index': A / 'tables/feature_presence_gene_index.tsv.gz',
    }
    dump(O / 'resources/INPUT_MANIFEST.json', {k:dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for k,f in inputs.items()})
    protocol = {
        'frozen_at': now(), 'scope': 'Existing expanded F1/F3/F4 labels in Phase66 8174-cell set, assigned B1/B2 using Phase71 sample metadata',
        'primary_contrasts': ['F1:B2-B1', 'F3:B2-B1', 'F4:B2-B1'],
        'no_reclustering_or_relabeling': True, 'no_F2_F5_analysis': True,
        'primary_population': 'Retain all selected cells and original lineage/CNV/confidence annotations; doublet flags preserved',
        'unit': 'Independent patient within source, never a cell or sample from the same patient',
        'primary_design': 'Within-source same-patient B2 minus B1; require measured cells in both groups for that class',
        'unmatched_data': 'All unmatched B2-only donors retained in descriptive expression; not mixed into a paired hypothesis test',
        'pseudobulk': 'Sum raw counts across cells and all same-group specimens within each patient; normalize by full original library UMI sum; log2(CP10k+0.1)',
        'effect': 'Equal-patient mean and median paired differences; positive means higher B2 expression',
        'gene_universe': 'All features measured in at least one relevant source sample; validity per patient-group requires measurement in every contributing sample',
        'missing_features': 'NA, not biological zero; per-gene effective N reported',
        'test': 'Exact two-sided patient sign-flip test; N<2 or all-zero paired effects => no P value',
        'FDR_family': 'Each class and source B2-vs-B1 full feature family separately; NA gets 1 only for BH family accounting, reported FDR stays NA',
        'candidate_gate': 'Target raw P<0.05, FDR reported, not mandatory for nomination',
        'small_N': 'Exploratory; report effect and donor consistency',
        'doublet_sensitivity': 'Check whether exclusions change any primary paired pseudobulk inputs; do not rerun identical inputs',
        'clinical_scope': 'B1 history-supported bone; B2 documented lung metastasis; contemporaneity and treatment matching not established',
        'cross_source_pooling': False, 'downstream_pathways_or_flow_selection': 'Not part of this step',
    }
    dump(O / 'resources/PROTOCOL_FROZEN.json', protocol)
    spec = importlib.util.spec_from_file_location('current_metadata', P / 'phase66_doublet_repair_DSC1_20261005/code/current_metadata.py')
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    cells = mod.load_current_cells().set_index('cell_id', drop=False)
    clinical = read(inputs['grouping'], keep_default_na=False).set_index('sample')
    cells['AB_group'] = cells['sample'].map(clinical.AB_group)
    cells['AB_evidence_status'] = cells['sample'].map(clinical.AB_evidence_status)
    assert cells.AB_group.notna().all()
    selected = cells[cells.class_id.isin(CLASSES) & cells.AB_group.isin(['B1','B2'])].copy()
    assert len(selected)==3241
    save(selected.reset_index(drop=True), 'F134_B1B2_cell_metadata.tsv.gz')
    a = ad.read_h5ad(inputs['counts'])
    assert a.shape==(8174,63398) and set(a.obs_names)==set(cells.index)
    assert (a.obs.class_id.astype(str)==cells.loc[a.obs_names,'class_id']).all()
    assert np.allclose(np.asarray(a.X.sum(axis=1)).ravel(), cells.loc[a.obs_names,'full_library_counts'].to_numpy(float))
    assert np.all(a.X.data>=0) and np.allclose(a.X.data,np.round(a.X.data))
    presence = sparse.load_npz(inputs['presence']).tocsr()
    sidx = read(inputs['sample_index']).set_index('sample').presence_row.to_dict()
    features = a.var_names.to_numpy(str); symbols = a.var.symbol.astype(str).to_numpy()
    gidx = read(inputs['feature_index'], keep_default_na=False)
    assert np.array_equal(gidx.feature.astype(str).to_numpy(),features)
    shutil.copy2(inputs['presence'], O / 'objects/sample_feature_presence.npz')
    shutil.copy2(inputs['sample_index'], O / 'resources/feature_presence_sample_index.tsv')
    shutil.copy2(inputs['feature_index'], O / 'resources/feature_presence_gene_index.tsv.gz')
    sample_rows=[]; overview=[]; roster=[]; availability=[]; all_stats=[]; primary_checks=[]; paired_effects=[]
    datasets=clinical.loc[clinical.AB_group.isin(['B1','B2']),'dataset'].unique()
    for cls in CLASSES:
        cs = selected[selected.class_id==cls]
        ix = a.obs_names.get_indexer(cs.index)
        ob=cs.drop(columns='cell_id').copy();ob.index.name='cell_id'
        for col in ob.select_dtypes(include='object').columns: ob[col]=ob[col].fillna('').astype(str)
        obj=ad.AnnData(X=a.X[ix].copy(),obs=ob,var=a.var.copy())
        obj.uns['counts']='raw full-feature counts; missing features require sample presence matrix'
        obj.uns['contrast']='B2 minus B1; patient unit, same-source paired primary'
        obj.uns['scope']='historical F label subset, not all-cell reannotation'
        obj.write_h5ad(O/'objects'/f'{cls}_B1B2_raw_counts.h5ad',compression='gzip')
        for grp in ['B1','B2']:
            z=cs[cs.AB_group==grp]
            save(z[['cell_id','sample','dataset','donor_record','class_id','AB_group','doublet_positive_any']].reset_index(drop=True),f'{cls}_{grp}_cell_index.tsv.gz')
            overview.append(dict(class_id=cls,AB_group=grp,cells=len(z),samples=z['sample'].nunique(),source_scoped_donor_records=z.donor_record.nunique(),doublet_positive_cells=int(truth(z.doublet_positive_any).sum())))
        for (ds,sm,grp),z in cs.groupby(['dataset','sample','AB_group']):
            sample_rows.append(dict(class_id=cls,dataset=ds,sample=sm,AB_group=grp,donor_record=z.donor_record.iloc[0],cells=len(z),doublet_positive_cells=int(truth(z.doublet_positive_any).sum())))
        for ds in datasets:
            z=cs[cs.dataset==ds]
            donors={grp:sorted(z.loc[z.AB_group==grp,'donor_record'].unique()) for grp in ['B1','B2']}
            pairs=sorted(set(donors['B1'])&set(donors['B2']))
            available='paired_exploratory' if len(pairs)>=2 else ('paired_direction_only' if len(pairs)==1 else ('B1_absent' if not donors['B1'] else 'no_complete_pair'))
            availability.append(dict(class_id=cls,dataset=ds,B1_patients=len(donors['B1']),B2_patients=len(donors['B2']),B1_cells=int((z.AB_group=='B1').sum()),B2_cells=int((z.AB_group=='B2').sum()),paired_patients=len(pairs),status=available,min_attainable_two_sided_P=2**(1-len(pairs)) if len(pairs)>=2 else np.nan))
            if z.empty: continue
            samples=z['sample'].unique().tolist()
            gi=np.flatnonzero(np.asarray(presence[[sidx[s] for s in samples]].sum(axis=0)).ravel()>0)
            fe=features[gi];sy=symbols[gi];m=len(gi)
            per={}; counts_by_patient={}; valid_by_patient={}
            for (donor,grp),q in z.groupby(['donor_record','AB_group']):
                idx=a.obs_names.get_indexer(q.index)
                raw=np.asarray(a.X[idx][:,gi].sum(axis=0)).ravel()
                lib=float(q.full_library_counts.sum())
                assert lib>0
                sm=q['sample'].unique().tolist()
                valid=np.asarray(presence[[sidx[s] for s in sm]][:,gi].sum(axis=0)).ravel()==len(sm)
                values=np.log2(raw/lib*1e4+0.1);values[~valid]=np.nan
                per[(donor,grp)]=values;counts_by_patient[(donor,grp)]=raw;valid_by_patient[(donor,grp)]=valid
                roster.append(dict(class_id=cls,dataset=ds,donor_record=donor,AB_group=grp,samples=';'.join(sm),cells=len(q),full_library_UMIs=lib,paired=donor in pairs,doublet_positive_cells=int(truth(q.doublet_positive_any).sum())))
            expr={g:np.vstack([per[(d,g)] for d in donors[g]]) if donors[g] else np.empty((0,m)) for g in ['B1','B2']}
            e1=valid_mean(expr['B1']);e2=valid_mean(expr['B2'])
            delta=np.vstack([per[(d,'B2')]-per[(d,'B1')] for d in pairs]) if pairs else np.empty((0,m))
            n=np.isfinite(delta).sum(axis=0)
            effect=valid_mean(delta)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                median=np.nanmedian(delta,axis=0) if pairs else np.full(m,np.nan)
            pv=np.full(m,np.nan);status=np.full(m,'not_tested_no_complete_pair',dtype=object)
            status[n==1]='descriptive_single_patient';status[n>=2]='not_tested_all_zero_effect'
            for j in range(m):
                d=delta[:,j];d=d[np.isfinite(d)]
                if len(d)<2 or np.all(d==0):continue
                signs=np.array(list(itertools.product([-1.,1.],repeat=len(d))))
                null=np.abs(signs@d/len(d));pv[j]=np.mean(null>=abs(d.mean())-1e-12)
                status[j]='exact_two_sided_patient_signflip'
            qv=multipletests(np.where(np.isfinite(pv),pv,1.),method='fdr_bh')[1]
            qv[~np.isfinite(pv)]=np.nan
            family=f'{ds}|{cls}|B2_vs_B1|paired_all_source_features'
            frame=pd.DataFrame(dict(class_id=cls,dataset=ds,feature=fe,symbol=sy,
                B1_all_donors_N=np.isfinite(expr['B1']).sum(axis=0),B2_all_donors_N=np.isfinite(expr['B2']).sum(axis=0),
                B1_all_donors_mean_log2CP10k=e1,B2_all_donors_mean_log2CP10k=e2,
                unpaired_group_mean_difference_descriptive_only=e2-e1,
                paired_N=n,paired_mean_B2_minus_B1=effect,paired_median_B2_minus_B1=median,
                paired_positive_patients=np.sum(delta>0,axis=0),paired_negative_patients=np.sum(delta<0,axis=0),
                P=pv,FDR=qv,test_status=status,family=family,nominal_P_lt_0_05=pv<.05))
            all_stats.append(frame)
            pair_names=np.array(pairs,dtype=str)
            np.savez_compressed(O/'objects'/f'{cls}__{ds}__patient_expression_effects.npz',features=fe,symbols=sy,B1_patients=np.array(donors['B1'],dtype=str),B2_patients=np.array(donors['B2'],dtype=str),B1_log2CP10k=expr['B1'],B2_log2CP10k=expr['B2'],paired_patients=pair_names,paired_effects=delta)
            for d,v in zip(pairs,delta):
                paired_effects.append(pd.DataFrame(dict(class_id=cls,dataset=ds,donor_record=d,feature=fe,symbol=sy,B2_minus_B1=v)))
            if pairs:
                paired_cells=z[z.donor_record.isin(pairs)]
                positives=int(truth(paired_cells.doublet_positive_any).sum())
                primary_checks.append(dict(class_id=cls,dataset=ds,paired_patients=len(pairs),paired_cells=len(paired_cells),doublet_positive_cells=positives,exclusion_changes_inputs=positives>0))
                assert positives==0, 'If paired input changes, add actual doublet sensitivity before completing'
            print('CONTEXT',cls,ds,'B1',len(donors['B1']),'B2',len(donors['B2']),'paired',len(pairs),'features',m,flush=True)
    save(pd.DataFrame(overview),'F134_group_summary.tsv')
    save(pd.DataFrame(sample_rows),'F134_sample_group_counts.tsv')
    save(pd.DataFrame(roster),'patient_pseudobulk_roster.tsv')
    save(pd.DataFrame(availability),'comparison_availability.tsv')
    save(pd.DataFrame(primary_checks),'doublet_sensitivity_input_check.tsv')
    stats=pd.concat(all_stats,ignore_index=True)
    for cls in CLASSES:
        f=stats[stats.class_id==cls]
        save(f,f'{cls}_B2_vs_B1_all_gene_statistics.tsv.gz')
        save(f[f.nominal_P_lt_0_05],f'{cls}_nominal_genes.tsv')
    save(pd.concat(paired_effects,ignore_index=True),'paired_patient_gene_effects.tsv.gz')
    # Independently reconstruct every source feature for samples in the primary paired comparisons.
    paired_roster=pd.DataFrame(roster);paired_roster=paired_roster[paired_roster.paired]
    rawqa=[]
    used_samples=sorted(set(s for v in paired_roster.samples for s in v.split(';')))
    for sm in used_samples:
        src=ad.read_h5ad(P/'phase26_client_reanalysis_20260922/inputs'/sm/'full_called_counts.h5ad',backed='r')
        for cls in paired_roster.loc[paired_roster.samples.str.split(';').map(lambda ss:sm in ss),'class_id'].unique():
            q=selected[(selected['sample']==sm)&(selected.class_id==cls)]
            oldix=src.obs_names.get_indexer(q.index);assert (oldix>=0).all()
            raw=src.X[np.sort(oldix),:].tocsr()
            newix=a.obs_names.get_indexer(q.index);vi=a.var_names.get_indexer(src.var_names);assert (vi>=0).all()
            actual=np.asarray(raw.sum(axis=0)).ravel();frozen=np.asarray(a.X[newix][:,vi].sum(axis=0)).ravel()
            err=float(np.max(np.abs(actual-frozen)));assert err==0
            assert float(raw.sum())==float(q.full_library_counts.sum())
            rawqa.append(dict(sample=sm,class_id=cls,cells=len(q),features=src.n_vars,max_raw_pseudobulk_error=err))
        src.file.close()
    save(pd.DataFrame(rawqa),'raw_source_pseudobulk_QA.tsv')
    # Explicit closed-form check of the N=2 primary test, independent of the generic sign generator.
    f4=stats[(stats.class_id=='F4')&(stats.dataset=='SCPCP000017')]
    e=np.load(O/'objects/F4__SCPCP000017__patient_expression_effects.npz')['paired_effects']
    ref=(np.abs(np.vstack([e[0]+e[1],e[0]-e[1],-e[0]+e[1],-e[0]-e[1]]))>=np.abs(e[0]+e[1])[None,:]-2e-12).mean(axis=0)
    tested=f4.P.notna().to_numpy();assert np.allclose(ref[tested],f4.P.to_numpy()[tested])
    assert stats.loc[stats.paired_N<2,'P'].isna().all()
    assert not stats.nominal_P_lt_0_05.any()
    for cls in CLASSES:
        sub=ad.read_h5ad(O/'objects'/f'{cls}_B1B2_raw_counts.h5ad',backed='r')
        assert set(sub.obs_names)==set(selected.index[selected.class_id==cls])
        sub.file.close()
    summary=[]
    for cls in CLASSES:
        f=stats[stats.class_id==cls]
        summary.append(dict(class_id=cls,statistic_rows=len(f),paired_effect_rows=int(f.paired_mean_B2_minus_B1.notna().sum()),valid_P_rows=int(f.P.notna().sum()),nominal_rows=int(f.nominal_P_lt_0_05.sum()),max_paired_N=int(f.paired_N.max())))
    save(pd.DataFrame(summary),'differential_analysis_summary.tsv')
    qa={'checked_at':now(),'status':'PASS','scope':'Existing F1 F3 F4; B2 versus B1 within class', 'selected_cells':len(selected),'group_summary':overview,'differential_summary':summary,'raw_source_pseudobulk_checks':len(rawqa),'raw_count_error':0,'F4_exact_P_closed_form_verified':True,'F1_no_inferential_P':True,'F3_no_B1_comparison':True,'unmatched_B2_not_added_as_paired_replicates':True,'no_nominal_genes':True,'doublet_exclusion_primary_inputs_unchanged':True,'missing_features_not_zero':True,'class_labels_preserved':True}
    dump(O/'QA/ACCEPTANCE.json',qa)
    text='F1/F3/F4 的 B1/B2 分组与差异分析\n\n范围：沿用Phase66的现有F1/F3/F4标签，映射Phase71临床分组。不是对全细胞重新定义F类。B1为肺转移病史支持的骨样本；B2为明确肺转移灶；同一患者多份肺样本先汇总。\n\n结果：F1，B1 3个细胞，B2 442个细胞；只有1对患者，输出全基因方向及效应，不计算推断性P值。F3，B1 0个细胞，B2 124个细胞；保留B2表达与分组结果，不把缺失B1当零表达，不计算差异P。F4，B1 102个细胞，B2 2570个细胞；主来源SCPCP000017有2对患者。按患者伪批量计算完整基因差异和精确双侧符号翻转检验，最小可达P为0.5，无原始P<0.05的条目。无名义差异不能解释为无生物学差异。\n\n统计：按来源分别分析。患者组内汇总完整原始计数，以完整文库归一化为CP10k，再取log2(CP10k+0.1)。主效应为同患者B2减B1；正值代表肺转移灶侧更高。不同来源及仅有B2的患者保留描述性数据，不混入配对检验。每个F类/来源的完整基因比较为一个BH家族，原始P和FDR均保留。N<2或全零效应不产生P值。\n\n组均值差列仅作描述，包含未配对患者；paired_mean_B2_minus_B1才是主配对效应。未测feature保留NA。主配对细胞均未标记为双细胞，排除双细胞不会改变主比较。其他B2细胞的标记保留在清单中。\n\nobjects/F*_B1B2_raw_counts.h5ad为3份分组后的原始计数对象，obs保留AB_group与原标签。其feature是跨来源并集，必须结合sample_feature_presence.npz及resources下索引区分未测基因。tables/F*_B2_vs_B1_all_gene_statistics.tsv.gz含各来源全基因统计与不可检验原因；paired_patient_gene_effects.tsv.gz保留逐患者效应。comparison_availability.tsv、patient_pseudobulk_roster.tsv和raw_source_pseudobulk_QA.tsv分别给出可评估性、患者清单和原始矩阵复核。\n\n本轮未新增聚类、富集、流式marker筛选。尚不能将F1/F3/F4差异解释为因果性转移效应，且B1/B2取材时序与治疗并未完全匹配。\n'
    (O/'delivery/RESULTS_README_ZH.txt').write_text(text,encoding='utf-8')
    files=[f for sub in ['code','tables','objects','resources','QA','delivery'] for f in (O/sub).rglob('*') if f.is_file() and f.suffix!='.zip']
    manifest=[dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=sha(f)) for f in files]
    dump(O/'delivery/FILE_MANIFEST.json',manifest)
    archive=O/'delivery/CAF_F134_B1B2_DE_20261009.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for f in files+[O/'delivery/FILE_MANIFEST.json']:z.write(f,str(f.relative_to(O)))
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    completed=dict(completed_at=now(),status='completed_verified',phase=O.name,archive=str(archive),archive_sha256=sha(archive),acceptance=str(O/'QA/ACCEPTANCE.json'),summary=summary,scope=protocol['scope'])
    dump(O/'COMPLETE.json',completed)
    print(json.dumps(completed,ensure_ascii=False),flush=True)

if __name__=='__main__':
    try:run()
    except BaseException:
        if O.exists():(O/'exit_code').write_text('1\n')
        traceback.print_exc();raise
    else:(O/'exit_code').write_text('0\n')
