"""Versioned NCC preprocessing; does not admit provisional cells to the CAF baseline."""
from pathlib import Path
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[key]='4'
import sys, json, datetime, hashlib, subprocess, concurrent.futures, gzip, traceback, math
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import h5py, anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase87_NCC_preprocessing_20261014'
RAW=P/'CAF_access_recheck_20261013/downloads/NCC'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    if isinstance(x,float) and not math.isfinite(x):return None
    return x
def dump(p,x):p.write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def read_counts(path):
    with h5py.File(path,'r') as h:
        m=h['matrix']; ft=m['features']
        ids=ft['id'].asstr()[:];syms=ft['name'].asstr()[:];bar=m['barcodes'].asstr()[:]
        x=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
        keep=ft['feature_type'].asstr()[:]=='Gene Expression'
        x=x[:,keep].copy();ids=ids[keep];syms=syms[keep]
    x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    assert len(set(ids))==len(ids) and len(set(bar))==len(bar)
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
    assert (np.asarray(x.sum(1)).ravel()>0).all()
    a=ad.AnnData(x,obs=pd.DataFrame(index=bar),var=pd.DataFrame({'symbol':syms},index=ids))
    a.obs_names.name='barcode';a.var_names.name='gene_id'
    return a
def worker(sm):
    print('START',sm,stamp(),flush=True)
    log=O/'logs'/(sm+'.log')
    with log.open('w',encoding='utf-8') as f:
        rc=subprocess.call([PY,str(O/'code/geo_qc_snapshot.py'),'--manifest',str(O/'manifest.tsv'),'--raw',str(O/'inputs'),'--out',str(O/'results'),'--sample',sm],stdout=f,stderr=subprocess.STDOUT)
    (O/'logs'/(sm+'.qc.exit_code')).write_text(str(rc)+'\n',encoding='utf-8')
    if rc!=0:return {'sample':sm,'status':'qc_failed','exit_code':rc}
    out=O/'results'/sm
    summary=json.loads((out/'complete.json').read_text(encoding='utf-8'))
    a=ad.read_h5ad(out/'annotated.h5ad')
    extra=0
    a.obs['doublet_final_method']='scrublet_automatic_threshold_pass'
    a.obs['doublet_assessed_any']=True
    a.obs['doublet_positive_any']=False
    if summary['doublet_status']!='automatic_threshold_pass':
        d=O/'doublet_review'/sm;d.mkdir(exist_ok=True)
        with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,a.layers['counts'].T,field='integer')
        (d/'genes.tsv').write_text('\n'.join(a.var_names)+'\n',encoding='utf-8')
        (d/'cells.tsv').write_text('\n'.join(a.obs_names)+'\n',encoding='utf-8')
        with (O/'logs'/(sm+'.scdblfinder.log')).open('w',encoding='utf-8') as f:
            rc=subprocess.call([R,'--vanilla',str(O/'code/scdblfinder.R'),str(d)],stdout=f,stderr=subprocess.STDOUT)
        (O/'logs'/(sm+'.scdblfinder.exit_code')).write_text(str(rc)+'\n',encoding='utf-8')
        if rc!=0:return {'sample':sm,'status':'alternate_doublet_failed','exit_code':rc}
        calls=pd.read_csv(d/'calls.tsv',sep='\t').set_index('cell_id')
        assert calls.index.is_unique and set(calls.index)==set(a.obs_names)
        calls=calls.loc[a.obs_names]
        assert calls.scDblFinder_class.isin(['singlet','doublet']).all()
        a.obs['scDblFinder_score']=calls.scDblFinder_score.to_numpy()
        a.obs['scDblFinder_class']=calls.scDblFinder_class.to_numpy()
        a.obs['doublet_final_method']='scDblFinder_alternate_wide_prior'
        a.obs['doublet_positive_any']=a.obs.predicted_doublet.fillna(False).to_numpy(bool)|(calls.scDblFinder_class.to_numpy()=='doublet')
        extra=int((calls.scDblFinder_class=='doublet').sum())
        a=a[calls.scDblFinder_class.to_numpy()=='singlet'].copy()
    a.obs['CAF_admission_status']='provisional_CNV_and_ambient_review_pending'
    a.obs['site_verified_for_Bone_Lung']=False
    a.uns['preprocessing_stage']=O.name
    a.uns['post_doublet_embedding_policy']='Initial QC coordinates retained; not a recomputed post-doublet embedding.'
    d=O/'accepted'/sm;d.mkdir(exist_ok=True)
    a.write_h5ad(d/'analysis_ready.h5ad',compression='gzip')
    a.obs.to_csv(d/'cell_metadata.tsv.gz',sep='\t',index_label='cell_id',encoding='utf-8')
    b=ad.read_h5ad(d/'analysis_ready.h5ad',backed='r')
    assert b.shape==a.shape and b.obs_names.tolist()==a.obs_names.tolist() and b.var_names.tolist()==a.var_names.tolist()
    b.file.close()
    counts=a.layers['counts']
    assert np.isfinite(counts.data).all() and np.equal(counts.data,np.floor(counts.data)).all()
    summary.update(status='preprocessing_complete_identity_pending',extra_doublets_excluded=extra,analysis_ready_cells=a.n_obs,
        fibroblast_like_cells=int((a.obs.lineage_prelim=='FibLike').sum()),
        mesenchymal_ambiguous_cells=int((a.obs.lineage_prelim=='MesAmb').sum()),
        doublet_positive_any_retained=int(a.obs.doublet_positive_any.sum()),
        confirmed_new_CAF_cells=None,accepted_sha256=sha(d/'analysis_ready.h5ad'))
    assert summary['input_cells']==summary['input_cells']-summary['qc_pass']+summary['doublets']+extra+a.n_obs
    dump(d/'ACCEPTANCE.json',summary)
    print('DONE',sm,a.n_obs,stamp(),flush=True)
    return summary
def main():
    for d in ['code','inputs','results','logs','accepted','doublet_review','tables']:(O/d).mkdir(parents=True,exist_ok=True)
    dump(O/'STATUS.json',{'status':'preparing_inputs','started_utc':stamp(),'pid':os.getpid(),'baseline_changed':False})
    accept=json.loads((P/'CAF_access_recheck_20261013/SERVER_ACCEPTANCE.json').read_text(encoding='utf-8'))
    expected={x['file']:x for x in accept['files']}
    meta=pd.read_csv(RAW/'Patient metadata',sep='\t')
    rows=[];integrity=[]
    for r in meta.to_dict('records'):
        bs=r['Sample']
        if bs=='BS6':continue
        source=RAW/(bs+'.h5');digest=sha(source)
        assert digest==expected[source.name]['sha256']
        a=read_counts(source);sm='NCC__'+bs
        d=O/'inputs'/sm;d.mkdir(exist_ok=True)
        a.write_h5ad(d/'counts.h5ad',compression='gzip')
        integrity.append({'sample':sm,'input_cells':a.n_obs,'features':a.n_vars,'total_counts':int(a.X.sum()),'raw_sha256':digest,'input_sha256':sha(d/'counts.h5ad')})
        row=dict(sample=sm,short=bs,patient_id='NCC__'+r['ID'],matrix_dir=str(d),dataset='NCC_Zheng2024',site='not_explicitly_verified',group='unassigned',disease_stage_raw=r['Tumor status'],treatment_raw=r['Treatment'],library_mapping_status='one_author_H5_per_sample_capture_details_unverified',source_evidence='author Patient metadata; fixed commit e5741ca70544508da3a09d6a99de47bb591ced9a',technology='10x_Gene_Expression_H5',patient_independence_status='within_NCC_ID_verified_cross_study_independence_unresolved',source_type='human_tumor',clinical_pairing='NCC__Patient4_pre_post' if r['ID']=='Patient4' else 'unpaired',cross_cohort_patient_pooling='not_allowed',Gender=r['Gender'],Age=r['Age'])
        rows.append(row)
        print('INPUT_VERIFIED',sm,a.shape,flush=True)
        del a
    assert len(rows)==9 and sum(x['input_cells'] for x in integrity)==79071
    assert len(set(x['patient_id'] for x in rows))==8
    pd.DataFrame(rows).to_csv(O/'manifest.tsv',sep='\t',index=False,encoding='utf-8')
    dump(O/'INPUT_ACCEPTANCE.json',{'status':'PASS','samples':9,'patient_labels':8,'input_cells':79071,'files':integrity,'BS6':'excluded_missing_gene_and_barcode_labels'})
    dump(O/'STATUS.json',{'status':'QC_doublet_lineage_running','pid':os.getpid(),'utc':stamp(),'samples':9,'workers':2,'baseline_changed':False,'formal_new_CAF_patients':0})
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        fs={ex.submit(worker,r['sample']):r['sample'] for r in rows}
        for fut in concurrent.futures.as_completed(fs):
            try:res=fut.result()
            except Exception:
                res={'sample':fs[fut],'status':'failed','traceback':traceback.format_exc()}
                print(res['traceback'],flush=True)
            results.append(res);dump(O/'SAMPLE_PROGRESS.json',results)
    pd.DataFrame(results).to_csv(O/'tables/sample_preprocessing_summary.tsv',sep='\t',index=False,encoding='utf-8')
    assert all(r['status']=='preprocessing_complete_identity_pending' for r in results),results
    # Rebuild global tables after isolated sample jobs; do not trust their shared summaries.
    tables=[pd.read_csv(O/'results'/r['sample']/'cluster_identity_evidence.tsv',sep='\t') for r in rows]
    pd.concat(tables,ignore_index=True).to_csv(O/'tables/cluster_identity_evidence.tsv',sep='\t',index=False,encoding='utf-8')
    final={'status':'preprocessing_complete_CNV_ambient_identity_pending','completed_utc':stamp(),'samples':9,'patient_labels':8,'input_cells':79071,'analysis_ready_cells':sum(r['analysis_ready_cells'] for r in results),'fibroblast_like_cells':sum(r['fibroblast_like_cells'] for r in results),'baseline_changed':False,'formal_new_CAF_patients':0,'next_step':'CNV with calibrated normal reference and holdouts, ambient sensitivity, lineage adjudication and cross-study patient deduplication','all_sample_exit_codes_zero':True,'all_output_readbacks_passed':True}
    dump(O/'PREPROCESSING_COMPLETE.json',final);dump(O/'STATUS.json',final)
    print('NCC_PREPROCESSING_COMPLETE',json.dumps(final),flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        dump(O/'FAILED.json',{'utc':stamp(),'traceback':traceback.format_exc()})
        dump(O/'STATUS.json',{'status':'failed','utc':stamp(),'see':'FAILED.json'})
        raise
