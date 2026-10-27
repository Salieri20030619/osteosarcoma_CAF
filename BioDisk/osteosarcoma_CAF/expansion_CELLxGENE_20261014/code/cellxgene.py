from pathlib import Path
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='2'
import json,subprocess,datetime,hashlib,traceback
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/expansion_CELLxGENE_20261014')
URL='https://datasets.cellxgene.cziscience.com/36d1ebc5-9ad3-4b6d-b321-eb1ceaa550cf.h5ad'
EXPECTED=2190278862
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def main():
    for d in ['downloads','objects','tables','logs']:(O/d).mkdir(parents=True,exist_ok=True)
    path=O/'downloads/extension_cohort.h5ad';tmp=path.with_suffix('.h5ad.partial')
    dump(O/'STATUS.json',{'status':'downloading','utc':stamp(),'pid':os.getpid(),'url':URL,'expected_bytes':EXPECTED})
    if not path.exists():
        with (O/'logs/curl.log').open('w',encoding='utf-8') as f:
            r=subprocess.run(['curl','--fail','--location','--show-error','--silent','--connect-timeout','20','--retry','3','--retry-delay','10','--speed-time','180','--speed-limit','1024','--continue-at','-','--output',str(tmp),URL],stdout=f,stderr=subprocess.STDOUT)
        (O/'logs/download.exit_code').write_text(str(r.returncode)+'\n',encoding='utf-8')
        assert r.returncode==0
        assert tmp.stat().st_size==EXPECTED
        tmp.rename(path)
    assert path.stat().st_size==EXPECTED
    digest=sha(path)
    dump(O/'DOWNLOAD_COMPLETE.json',{'url':URL,'bytes':EXPECTED,'sha256':digest,'size_matches_official_asset':True,'utc':stamp()})
    dump(O/'STATUS.json',{'status':'extracting_osteosarcoma','utc':stamp(),'pid':os.getpid()})
    a=ad.read_h5ad(path,backed='r')
    assert a.n_obs==172370
    obs=a.obs
    mask=obs['disease_ontology_term_id'].astype(str).eq('MONDO:0009807')
    assert mask.any()
    b=a[np.flatnonzero(mask)].to_memory();a.file.close()
    assert set(b.obs.donor_id.astype(str))=={'EXT_OS_1','EXT_OS_2'}
    b.write_h5ad(O/'objects/osteosarcoma_author_object.h5ad',compression='gzip')
    b.obs.to_csv(O/'tables/osteosarcoma_author_metadata.tsv.gz',sep='\t',index_label='cell_id',encoding='utf-8')
    assert b.raw is not None,'Declared raw.X missing'
    x=sp.csr_matrix(b.raw.X);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
    var=b.raw.var.copy()
    assert 'feature_name' in var.columns
    var['symbol']=var['feature_name'].astype(str)
    counts=ad.AnnData(X=x,obs=b.obs.copy(),var=var)
    counts.obs['source']='CELLxGENE_Whitfield2026_10x_multiome_nuclei'
    counts.obs['CAF_status']='author_stromal_label_is_not_CAF_confirmation'
    counts.obs['n_counts']=np.asarray(x.sum(1)).ravel()
    counts.obs['n_genes']=x.getnnz(axis=1)
    for name,genes in {'FibLike':['DCN','LUM','PDGFRA','DPT','COL3A1','COL6A1'],'OsteoLike':['RUNX2','SP7','ALPL','IBSP','BGLAP','DMP1'],'Mural':['RGS5','CSPG4','NOTCH3','MCAM','MYH11'],'Immune':['PTPRC','CD3D','CD79A','LYZ','LST1'],'Endo':['PECAM1','VWF','CLDN5','EMCN','KDR']}.items():
        detected=np.zeros(counts.n_obs,dtype=int)
        for gene in genes:
            ix=np.flatnonzero(var.symbol.to_numpy()==gene)
            if len(ix):detected+=(np.asarray(x[:,ix].sum(1)).ravel()>0)
        counts.obs['det_'+name]=detected
    counts.write_h5ad(O/'objects/osteosarcoma_counts.h5ad',compression='gzip')
    counts.obs.to_csv(O/'tables/osteosarcoma_qc_and_lineage_evidence.tsv.gz',sep='\t',index_label='cell_id',encoding='utf-8')
    groups=[k for k in ['donor_id','tissue','cell_type'] if k in counts.obs]
    c=counts.obs.groupby(groups,observed=True).size().reset_index(name='cells')
    c.to_csv(O/'tables/author_cell_counts.tsv',sep='\t',index=False,encoding='utf-8')
    summary={'status':'OS_counts_extracted_identity_review_pending','utc':stamp(),'input_cohort_cells':172370,'OS_cells':b.n_obs,'donor_labels':sorted(set(b.obs.donor_id.astype(str))),'raw_count_features':counts.n_vars,'author_cell_counts':c.to_dict('records'),'raw_integer_counts_verified':True,'counts_object_sha256':sha(O/'objects/osteosarcoma_counts.h5ad'),'formal_new_CAF_patients':0,'clinical_grouping':'pending_per_donor_anatomical_and_treatment_metadata','modality':'single_nucleus_10x_multiome','baseline_changed':False}
    dump(O/'EXTRACTION_COMPLETE.json',summary);dump(O/'STATUS.json',summary)
    print(json.dumps(summary),flush=True)
if __name__=='__main__':
    rc=1
    try:main();rc=0
    except Exception:dump(O/'FAILED.json',{'utc':stamp(),'traceback':traceback.format_exc()});traceback.print_exc()
    (O/'logs/pipeline.exit_code').write_text(str(rc)+'\n',encoding='utf-8')
    dump(O/'JOB_FINISHED.json',{'exit_code':rc,'success':rc==0 and (O/'EXTRACTION_COMPLETE.json').exists(),'utc':stamp()})
    raise SystemExit(rc)
