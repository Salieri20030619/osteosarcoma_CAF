"""Check saved counts, generate descriptive markers and lineage review packets."""
from pathlib import Path
import os,json,hashlib,subprocess,sys,datetime,gzip
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918')
manifest=pd.read_csv(P/'manifest.tsv',sep='\t').fillna('')
assert len(manifest)==16
assert (P/'results/SUCCESS.txt').read_text().strip()=='GEO_QC_CLUSTERING_COMPLETE'
audits=[];coverage=[];summary=[]
gene_order=pd.read_csv(P.parent/'phase6_nonpdx_cnv_20260916/resources/infercnv_gene_order_GRCh38_98.tsv',sep='\t',header=None)
ordered_genes=set(gene_order.iloc[:,0].astype(str))
feature_file=next((P.parent/'phase7_review_states_20260917/download_expansion_20260917/files/GSE162454').glob('*features.tsv.gz'))
symbol_ids={}
with gzip.open(feature_file,'rt') as f:
    for line in f:
        fields=line.rstrip().split('\t');symbol_ids.setdefault(fields[1],set()).add(fields[0])
symbol_to_ens={g:next(iter(ids)) for g,ids in symbol_ids.items() if len(ids)==1}

for r in manifest.to_dict('records'):
    sample=r['sample'];out=P/'results'/sample
    comp=json.loads((out/'complete.json').read_text())
    a=ad.read_h5ad(out/'annotated.h5ad');x=a.layers['counts'].tocsr()
    raw=ad.read_h5ad(P/'inputs'/sample/'counts.h5ad')
    assert a.var_names.equals(raw.var_names)
    pos=raw.obs_names.get_indexer(a.obs.barcode.astype(str));assert (pos>=0).all()
    assert (x!=raw.X[pos,:]).nnz==0,'saved raw counts changed'
    allobs=pd.read_csv(out/'all_cells_qc.tsv.gz',sep='\t')
    assert len(allobs)==raw.n_obs==comp['input_cells']
    assert int(allobs.analysis_keep.sum())==a.n_obs==comp['retained']
    assert a.obs_names.is_unique and a.var_names.is_unique and x.nnz>0
    assert np.isfinite(a.X.data).all() and np.isfinite(x.data).all()
    assert np.equal(x.data,np.floor(x.data)).all() and (x.data>=0).all()
    assert np.isfinite(a.obsm['X_umap']).all()
    assert a.obs.patient_id.nunique()==1 and str(a.obs.patient_id.iloc[0])==r['patient_id']
    norm=sp.diags(1e4/np.asarray(x.sum(1)).ravel()).dot(x).tocsr()
    sums=np.asarray(norm.sum(0)).ravel();dets=np.asarray((x>0).sum(0)).ravel()
    records=[]
    for cl in a.obs.cluster.astype(str).unique():
        mask=(a.obs.cluster.astype(str)==cl).to_numpy();n=int(mask.sum());nr=a.n_obs-n
        if not nr:continue
        sm=np.asarray(norm[mask].sum(0)).ravel();dt=np.asarray((x[mask]>0).sum(0)).ravel()
        mean=sm/n;other=(sums-sm)/nr;pct=dt/n;pct_other=(dets-dt)/nr
        lfc=np.log2((mean+0.01)/(np.maximum(other,0)+0.01))
        score=lfc*np.maximum(pct-pct_other,0)
        eligible=np.flatnonzero((pct>=.1)&(lfc>0)&(pct>pct_other))
        top=eligible[np.argsort(-score[eligible],kind='stable')[:50]]
        for j in top:
            records.append(dict(group=cl,names=str(a.var_names[j]),symbol=str(a.var.symbol.iloc[j]),scores=float(score[j]),
                log2FC=float(lfc[j]),pct_in=float(pct[j]),pct_out=float(pct_other[j]),
                mean_normalized_in=float(mean[j]),mean_normalized_out=float(other[j]),
                marker_method='descriptive_log2FC_times_detection_difference_no_P_values'))
    pd.DataFrame(records).to_csv(out/'cluster_markers.tsv',sep='\t',index=False)
    symbols=a.var.symbol.astype(str)
    coordinate_ids=[g if g.startswith('ENSG') else symbol_to_ens.get(s,'') for g,s in zip(a.var_names,symbols)]
    match=np.array([g in ordered_genes for g in coordinate_ids])
    detected=np.asarray((x>0).sum(0)).ravel()>0
    thresholds=json.loads((out/'qc_thresholds.json').read_text())
    entry=dict(comp)
    entry.pop('fingerprint',None)
    entry.update(site=r['site'],technology=r['technology'],patient_independence_status=r['patient_independence_status'],
        min_genes=thresholds['min_genes'],min_counts=thresholds['min_counts'],max_pct_mt=thresholds['max_pct_mt'],
        qc_pass_without_mt=int(allobs.qc_pass_without_mt.sum()),mt_only_exclusions=int((allobs.qc_pass_without_mt & ~allobs.qc_pass).sum()),
        low_genes=int(allobs.qc_low_genes.sum()),low_counts=int(allobs.qc_low_counts.sum()),high_mt=int(allobs.qc_high_mt.sum()),
        retained_fraction=a.n_obs/raw.n_obs,
        gene_order_covered=int(match.sum()),gene_order_detected=int((match&detected).sum()),
        source_count_conservation=True,CAF_identity='NOT_CONFIRMED_CNV_PENDING')
    summary.append(entry)
    audits.append(dict(sample=sample,raw_counts_identical=True,cell_accounting_pass=True,identifiers_pass=True,finite_embedding_pass=True,retained=a.n_obs))
    del raw,a,x,norm
    print('OUTPUT_VERIFIED',sample,comp['retained'],flush=True)
pd.DataFrame(summary).to_csv(P/'review/geo_processing_summary.tsv',sep='\t',index=False)
(P/'review/output_integrity.json').write_text(json.dumps(audits,indent=2)+'\n')
# Reuse broad lineage marker panels and rules; top markers are descriptive.
source=P.parent/'phase5_nonpdx_qc_20260916/code/review_lineages.py'
script=source.read_text()
script=script.replace("r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase5_nonpdx_qc_20260916')", "r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918')")
script=script.replace('assert len(manifest)==17','assert len(manifest)==16')
script=script.replace('PHASE5_LINEAGE_REVIEW_READY','GEO_LINEAGE_REVIEW_READY')
(P/'code/geo_lineage_review.py').write_text(script)
subprocess.run([sys.executable,'-m','py_compile',str(P/'code/geo_lineage_review.py')],check=True)
subprocess.run([sys.executable,'-u',str(P/'code/geo_lineage_review.py')],check=True)
(P/'review/marker_method.json').write_text(json.dumps({'ranking':'log2((mean normalized expression + 0.01)/(rest mean + 0.01)) * positive detection difference','minimum_within_cluster_detection':0.1,'maximum_markers_per_cluster':50,'P_values_computed':False,'source_lineage_script':str(source),'source_lineage_script_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'interpretation':'descriptive cluster annotation and normal-reference feasibility; not patient-level differential expression or confirmed CAF identity'},indent=2))
(P/'review/OUTPUTS_VERIFIED').write_text('GEO_COUNTS_QC_LINEAGE_OUTPUTS_VERIFIED\n')
print('GEO_REVIEW_COMPLETE',flush=True)
