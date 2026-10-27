"""Validate all 16 GEO count matrices and prepare a versioned QC adapter."""
from pathlib import Path
import os, json, gzip, hashlib, csv, datetime, shutil
os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
os.environ.setdefault('OMP_NUM_THREADS', '4')
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import anndata as ad
from threadpoolctl import threadpool_limits
threadpool_limits(4)
ROOT = Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
STAGE = ROOT/'geo_preprocess_20260918'
SOURCE = ROOT/'phase7_review_states_20260917'
RAW = SOURCE/'download_expansion_20260917/files'
def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
rows=[r for r in json.loads((SOURCE/'candidate_sources_20260917/GEO_candidate_sample_manifest.json').read_text()) if r['dataset'] in ['GSE318840','GSE162454','GSE250015']]
assert len(rows)==16
# Prioritize the two new lung tissues, then other FLEX tissues and bone cohorts.
rows.sort(key=lambda r:(0 if r['sample'] in ['GSM9504303','GSM9504304'] else 1 if r['dataset']=='GSE318840' else 2 if r['dataset']=='GSE162454' else 3,r['sample']))
checks=json.loads((SOURCE/'download_expansion_20260917/geo_download_completion.json').read_text())
expected={r['url'].split('/')[-1]:r for r in checks['results']}
file_checks=[]; manifest=[]; audits=[]
for r in rows:
    gsm=r['sample']; ds=r['dataset']; sample=ds+'__'+gsm
    files={u.split('/')[-1]:RAW/ds/u.split('/')[-1] for u in r['supplementary_files']}
    for name,f in files.items():
        digest=sha(f); assert digest==expected[name]['sha256'],f
        assert f.stat().st_size==expected[name]['bytes'],f
        file_checks.append(dict(path=str(f),bytes=f.stat().st_size,sha256=digest,verified=True))
    if ds!='GSE250015':
        f=next(f for n,f in files.items() if 'features.tsv' in n)
        b=next(f for n,f in files.items() if 'barcodes.tsv' in n)
        m=next(f for n,f in files.items() if 'matrix.mtx' in n)
        with gzip.open(f,'rt') as h: features=[l.rstrip('\n').split('\t') for l in h]
        with gzip.open(b,'rt') as h: barcodes=[l.strip() for l in h]
        counts=mmread(str(m)).tocsr(); counts.sum_duplicates(); counts.eliminate_zeros(); counts.sort_indices()
        assert counts.shape==(len(features),len(barcodes))
        keep=np.array([len(f)<3 or f[2]=='Gene Expression' for f in features])
        excluded=int((~keep).sum()); features=[f for f,k in zip(features,keep) if k]
        ids=[f[0] for f in features]; symbols=[f[1] for f in features]
        counts=counts[keep,:]
        id_type='Ensembl'
    else:
        f=next(iter(files.values()))
        with gzip.open(f,'rt') as h: barcodes=next(csv.reader(h))[1:]
        pieces=[]; symbols=[]
        for frame in pd.read_csv(f,index_col=0,chunksize=256):
            assert list(frame.columns)==barcodes
            v=frame.to_numpy(dtype=np.float64)
            assert np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all()
            symbols.extend(frame.index.astype(str)); pieces.append(sp.csr_matrix(v))
        counts=sp.vstack(pieces,format='csr'); counts.eliminate_zeros(); counts.sort_indices()
        ids=symbols.copy(); excluded=0; id_type='gene_symbol'
    assert len(ids)==len(set(ids)),(gsm,'duplicate gene identifiers')
    assert len(barcodes)==len(set(barcodes)),(gsm,'duplicate cell barcodes')
    assert np.isfinite(counts.data).all() and (counts.data>=0).all() and np.equal(counts.data,np.floor(counts.data)).all()
    assert counts.data.max()<2**24,'float32 count precision would be lost'
    x=counts.T.tocsr().astype(np.float32)
    assert np.array_equal(np.asarray(x.sum(axis=1,dtype=np.float64)).ravel(),np.asarray(counts.sum(axis=0)).ravel())
    a=ad.AnnData(x,obs=pd.DataFrame(index=barcodes),var=pd.DataFrame({'symbol':symbols,'source_gene_id':ids},index=ids))
    a.uns['source_files']=[str(f) for f in files.values()]
    a.uns['count_validation']='nonnegative_integer_finite_dimensions_and_identifiers_passed'
    a.uns['source_identifier_type']=id_type
    path=STAGE/'inputs'/sample; path.mkdir(exist_ok=True)
    a.write_h5ad(path/'counts.h5ad',compression='gzip')
    mt=np.array([g.startswith('MT-') for g in symbols]); total=np.asarray(x.sum(1)).ravel()
    mt_counts=np.asarray(x[:,mt].sum(1)).ravel()
    pct=np.divide(100*mt_counts,total,out=np.zeros_like(total),where=total>0)
    chars=dict(v.split(': ',1) for v in r['characteristics'] if ': ' in v)
    if ds=='GSE318840':
        site={'Bone':'bone','Lung':'lung','Soft tissue':'soft_tissue'}[r['source']]
        group=r['title'].split('_')[1]; tech='10x_FLEX_singleplex_probe_based'
        patient=ds+'__'+r['title']; mapping='series_reports_eight_distinct_patients_singleplex_sample'
        treatment=chars['treatment']; independence='within_cohort_8_distinct_patients_reported_cross_cohort_not_verified'
    elif ds=='GSE162454':
        site='bone';group='primary_tumor';tech='10x_3prime_v3';patient=ds+'__Patient_'+r['title'].split()[1]
        mapping='series_reports_six_distinct_patients_individual_capture_mapping_unavailable'
        treatment='not_reported_in_GEO';independence='within_cohort_6_distinct_patients_reported_cross_cohort_not_verified'
    else:
        site='bone';group=chars['phenotype'];tech='Singleron_GEXSCOPE';patient=ds+'__'+gsm
        mapping='one_reported_primary_lesion_per_group_individual_capture_mapping_unavailable'
        treatment='not_reported_in_GEO';independence='cross_cohort_overlap_with_existing_Shanghai_cohorts_unresolved'
    manifest.append(dict(sample=sample,short=gsm,dataset=ds,patient_id=patient,site=site,group=group,
        disease_stage_raw=r['title'],treatment_raw=treatment,library_mapping_status=mapping,
        source_evidence='GEO_SOFT:'+ds,technology=tech,matrix_dir=str(path),
        patient_independence_status=independence,source_type='direct_patient_tumor',cnv_status='not_run',
        clinical_pairing='not_established',cross_cohort_patient_pooling='not_authorized_by_identity_evidence'))
    audits.append(dict(sample=sample,dataset=ds,cells=a.n_obs,features=a.n_vars,nnz=x.nnz,total_counts=int(total.astype(np.float64).sum()),
        identifier_type=id_type,duplicate_symbols=int(pd.Series(symbols).duplicated().sum()),nonRNA_removed=excluded,
        mitochondrial_features=int(mt.sum()),mitochondrial_detected_features=int((np.asarray(x[:,mt].sum(0)).ravel()>0).sum()),
        pct_mt_median=float(np.median(pct)),pct_mt_p99=float(np.quantile(pct,.99)),
        cells_zero_counts=int((total==0).sum()),matrix_status='PASS',site=site,technology=tech))
    dump(STAGE/'review/matrix_audit.json',audits)
    print('MATRIX_ACCEPTED',gsm,a.shape,'mt_median',float(np.median(pct)),flush=True)
pd.DataFrame(manifest).to_csv(STAGE/'manifest.tsv',sep='\t',index=False)
pd.DataFrame(audits).to_csv(STAGE/'review/matrix_audit.tsv',sep='\t',index=False)
dump(STAGE/'review/file_integrity.json',file_checks)
assert len(file_checks)==44
source=ROOT/'phase5_nonpdx_qc_20260916/code/phase5.py'
script=source.read_text()
start=script.index('def read_matrix(path):');end=script.index('\ndef mad(x):',start)
script=script[:start]+'''def read_matrix(path):
 a=ad.read_h5ad(path/'counts.h5ad')
 assert sp.issparse(a.X) and a.obs_names.is_unique and a.var_names.is_unique
 return a,0
'''+script[end:]
# Retain the original sample-level thresholds and embeddings, but suppress
# cell-level inferential marker P values: use descriptive mean/detection tables.
start=script.index(' markers=[]\n');end=script.index(' rows=[];labelmap={}',start)
script=script[:start]+script[end:]
old="a.obs['technology_reported']=row['technology'];a.obs['contains_peaks']=peaks>0"
new=old+"\n for key in ['patient_independence_status','source_type','clinical_pairing','cross_cohort_patient_pooling']:a.obs[key]=row[key]"
assert old in script;script=script.replace(old,new)
old="'capture_mapping':'author tissue matrix; individual capture mapping unresolved'"
assert old in script
script=script.replace(old,"'capture_mapping':row['library_mapping_status'],'doublet_prior_policy':'legacy 0.06 prior, exploratory sample-level assessment; not a measured loading rate','technology':row['technology'],'mt_policy':'sample-specific empirical thresholds; platform percentages not compared as biological effects'")
script=script.replace(" a.obs['qc_pass']=~(a.obs.qc_low_genes|a.obs.qc_low_counts|a.obs.qc_high_mt)"," a.obs['qc_pass_without_mt']=~(a.obs.qc_low_genes|a.obs.qc_low_counts)\n a.obs['qc_pass']=~(a.obs.qc_low_genes|a.obs.qc_low_counts|a.obs.qc_high_mt)")
script=script.replace('PHASE5_QC_CLUSTERING_COMPLETE','GEO_QC_CLUSTERING_COMPLETE')
(STAGE/'code/geo_qc.py').write_text(script)
dump(STAGE/'review/pipeline_provenance.json',{'prepared_utc':now(),'source_script':str(source),'source_sha256':sha(source),'adapter_sha256':sha(STAGE/'code/geo_qc.py'),'changes':['validated h5ad input adapter','metadata preservation','no cell-level inferential marker tests','record MT-free QC sensitivity and platform-specific interpretation'],'seed':270231,'threads':4,'max_concurrent_samples':1,'doublet_prior':0.06,'ScPCA':'client downloading; not processed'})
(STAGE/'inputs/COMPLETE').write_text('16_MATRICES_44_FILES_VALIDATED\n')
print('GEO_INPUT_ACCEPTANCE_COMPLETE',flush=True)
