"""SYNTAX_CHECK: python -m py_compile takeover_audit.py
RUN_COMMAND: python takeover_audit.py --root <OS_CAF_phase2> --out <audit>
EXPECTED_SUCCESS_MARKER: TAKEOVER_AUDIT_COMPLETE
"""
import argparse,json,csv,hashlib,platform
from pathlib import Path
import h5py,numpy as np,pandas as pd
from scipy.sparse import csr_matrix

def arr(x):
    if isinstance(x,h5py.Group):
        if 'values' in x:
            a=arr(x['values']).astype(object)
            if 'mask' in x:a[arr(x['mask'])]=None
            return a
        if 'codes' in x:
            c=arr(x['categories']);return np.array([c[i] if i>=0 else None for i in arr(x['codes'])])
        raise ValueError(list(x))
    a=x[()]
    if getattr(a,'dtype',None) is not None and a.dtype.kind in 'SO':a=np.array([v.decode() if isinstance(v,bytes) else v for v in a])
    return a
def frame(g):
    idx=g.attrs.get('_index','_index');idx=idx.decode() if isinstance(idx,bytes) else idx
    d={}
    for k in g:
        if k==idx:continue
        x=g[k]
        if isinstance(x,h5py.Dataset):d[k]=arr(x)
        elif 'codes' in x and 'categories' in x:
            cat=arr(x['categories']);code=arr(x['codes']);d[k]=np.array([cat[c] if c>=0 else None for c in code])
        elif 'values' in x:
            v=arr(x['values']).astype(object);v[arr(x['mask'])]=None;d[k]=v
        else:raise ValueError((k,list(x)))
    return pd.DataFrame(d,index=arr(g[idx]))
def sparse(g):return csr_matrix((g['data'][()],g['indices'][()],g['indptr'][()]),shape=tuple(g.attrs['shape']))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    s=pd.read_csv(args.root/'results/sample_summary.tsv',sep='\t');details=[];cl=[];db=[]
    expected=[3086,6334,3567,4434,7952,925,2244]
    allcols=[]
    for i,r in s.iterrows():
        p=args.root/'results'/r['sample'];f=p/'annotated.h5ad'
        with h5py.File(f) as h:
            o=frame(h['obs']);v=frame(h['var']);x=sparse(h['layers/counts'])
            assert len(o)==r.retained==expected[i]
            assert len(v)==33694 and x.shape==(len(o),len(v))
            assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all()
            assert o.index.is_unique and v.index.is_unique
            assert (o['sample']==r['sample']).all()
            assert np.allclose(np.asarray(x.sum(1)).ravel(),o.n_counts)
            assert 'X_umap' in h['obsm'] and np.isfinite(h['obsm/X_umap'][()]).all()
            req=['lineage_prelim','lineage_reviewed','cnv_ref_immune_reviewed','cnv_ref_endo_reviewed','cnv_target_reviewed']
            assert set(req)<=set(o.columns)
            q=pd.read_csv(p/'cell_annotations_reviewed.tsv.gz',sep='\t',index_col=0)
            assert o.index.tolist()==q.index.tolist()
            for k in req:assert o[k].astype(str).tolist()==q[k].astype(str).tolist(),(r['short'],k)
            X=sparse(h['X']);assert np.isfinite(X.data).all()
            norm=x.copy().astype(float);tot=np.asarray(x.sum(1)).ravel();norm.data=np.log1p(norm.data*np.repeat(1e4/tot,np.diff(norm.indptr)))
            diff=X-norm;maxerr=float(np.abs(diff.data).max()) if diff.nnz else 0
            assert maxerr<1e-4
            details.append({'sample':r['sample'],'short':r['short'],'cells':len(o),'genes':len(v),'counts_nnz':x.nnz,'counts_integer':True,'counts_obs_sum_match':True,'review_tsv_match':True,'normalization_max_error':maxerr,'ensembl_fraction':float(v.index.str.startswith('ENSG').mean()),'duplicate_symbol_rows':int(v.symbol.duplicated(keep=False).sum()),'patient_ids':sorted(o.patient_id.astype(str).unique().tolist()),'cnv_status':o.cnv_status.unique().tolist(),'immune_ref':int(o.cnv_ref_immune_reviewed.sum()),'endo_ref':int(o.cnv_ref_endo_reviewed.sum()),'targets':int(o.cnv_target_reviewed.sum()),'lineage_counts':o.lineage_reviewed.value_counts().to_dict(),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'tmp_files':[a.name for a in p.glob('*tmp*')],'obs_columns':o.columns.tolist(),'var_columns':v.columns.tolist(),'uns_keys':list(h['uns'])})
            for c,g in o.groupby('cluster'):
                cl.append({'sample':r['sample'],'short':r['short'],'cluster':c,'cells':len(g),'label':g.lineage_reviewed.iloc[0],'immune_ref':int(g.cnv_ref_immune_reviewed.sum()),'endo_ref':int(g.cnv_ref_endo_reviewed.sum()),'targets':int(g.cnv_target_reviewed.sum()),'median_umi':g.n_counts.median(),'median_genes':g.n_genes.median(),'median_db':g.doublet_score.median(),'high_complexity':int(g.high_complexity_flag.sum())})
        qc=pd.read_csv(p/'all_cells_qc.tsv.gz',sep='\t');qt=qc[qc.qc_pass]
        assert len(qc)==r.input_cells and len(qt)==r.qc_pass and qc.analysis_keep.sum()==r.retained
        db.append({'sample':r['sample'],'short':r['short'],'threshold':r.scrublet_threshold,'qc_pass':len(qt),'doublets':r.doublets,**{f'q{int(k*100)}':float(qt.doublet_score.quantile(k)) for k in [0,.5,.9,.95,.99,1]}})
        print('AUDIT_OK',r['short'],len(o),flush=True)
    pd.DataFrame(cl).to_csv(args.out/'cluster_audit.tsv',sep='\t',index=False)
    pd.DataFrame(db).to_csv(args.out/'doublet_audit.tsv',sep='\t',index=False)
    (args.out/'object_audit.json').write_text(json.dumps(details,ensure_ascii=False,indent=2),encoding='utf-8')
    print('TAKEOVER_AUDIT_COMPLETE')
if __name__=='__main__':main()
