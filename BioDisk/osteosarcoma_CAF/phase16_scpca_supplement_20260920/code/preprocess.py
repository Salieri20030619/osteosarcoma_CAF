"""ScPCA supplementary QC; no CAF confirmation before ambient and dual CNV review."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']: os.environ[key]='4'
import json, datetime, hashlib, traceback, gzip
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase16_scpca_supplement_20260920'
for folder in ['code','logs','tables','objects','inputs','resources']: (O/folder).mkdir(parents=True,exist_ok=True)
panels=json.loads((P/'phase10_batch_gates_20260919/resources/panels.json').read_text())
customer=json.loads((P/'phase14_client_objectives_20260920/resources/customer_panels.json').read_text())
surfaces=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
def dump(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str));tmp.replace(path)
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''): h.update(block)
    return h.hexdigest()
state={'status':'INPUT_ACCEPTANCE_AND_QC_RUNNING','started_utc':now(),'pid':os.getpid(),'completed_samples':[],'strict_CAF_confirmed':None}
dump(O/'STATUS.json',state)
protocol={'source':'filtered X verified integer counts; official processed barcode membership retained as mandatory primary QC gate',
 'additional_qc':'same existing 3 unscaled MAD rule on log10 genes/counts; min200 genes,min300 UMI,MT max(min20,max5,median+3MAD)',
 'doublets':'retain official scDblFinder singlets only; do not rerun or override official rejection',
 'embedding':'per sample seurat HVG3000, exclude MT/ribosomal/sex drivers, PCA40, Leiden .4/.8/1.2; no integrated expression used for gene measurement',
 'lineage':'descriptive panel evidence and initial cluster labels; full marker review, ambient and dual-CNV acceptance required before strict CAF',
 'customer':'exact client panels, raw-count descriptive support >=2 genes; not final myCAF identity',
 'statistics':'no cell-level inferential tests; same-patient samples are paired records, never independent patient replicates',
 'clinical':'do not equate shared patient with synchronous collection; chemotherapy No only as defined in source; cross-study independence unresolved',
 'frozen_at':now()}
dump(O/'PROTOCOL.json',protocol);dump(O/'resources/customer_panels.json',customer)
manifest=[]; meta=[]
for root in [P/'phase15_scpca_download_20260920',P/'phase15_scpca_download_20260920/additional23']:
    assert json.loads((root/'DOWNLOAD_QA.json').read_text())['all_checks_passed']
    records=json.loads((root/'verified_files.json').read_text())
    m=pd.read_csv(root/'sample_download_acceptance.tsv',sep='\t',keep_default_na=False)
    for row in m.to_dict('records'):
        sm=row['scpca_sample_id']; rr={}
        for kind in ['filtered','processed']:
            r=next(r for r in records if sm in r['archive_member'] and '_'+kind+'_rna.h5ad' in r['archive_member'])
            rr[kind]=r['path'];rr[kind+'_sha256']=r['sha256']
        manifest.append({**row,**rr});meta.append(row)
assert len(manifest)==12 and len({r['scpca_sample_id'] for r in manifest})==12
pd.DataFrame(manifest).to_csv(O/'manifest.tsv',sep='\t',index=False)
md=pd.DataFrame(meta)
md.groupby(['scpca_project_id','participant_id']).agg(samples=('scpca_sample_id',';'.join),sites=('tissue_location',';'.join),records=('scpca_sample_id','size')).reset_index().to_csv(O/'tables/patient_sample_links.tsv',sep='\t',index=False)
summ=[];panelrows=[];clusterrows=[]
try:
 for row in manifest:
    sm=row['scpca_sample_id'];state.update(current_sample=sm);dump(O/'STATUS.json',state);print('SAMPLE_START',sm,now(),flush=True)
    for kind in ['filtered','processed']: assert sha(row[kind])==row[kind+'_sha256'], (sm,kind,'source hash mismatch')
    src=ad.read_h5ad(row['filtered']);proc=ad.read_h5ad(row['processed'],backed='r');accepted=set(proc.obs_names);proc.file.close()
    x=sparse.csr_matrix(src.X);x.eliminate_zeros()
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data),atol=1e-6)
    assert src.uns.get('X_name','').lower() in ['counts','count','raw_counts'], str(src.uns.get('X_name'))
    assert src.obs_names.is_unique and src.var_names.is_unique and accepted<=set(src.obs_names)
    sy=src.var.gene_symbol.astype(str).to_numpy();lib=np.asarray(x.sum(1)).ravel();ng=x.getnnz(1)
    assert np.allclose(lib,src.obs['sum'].to_numpy(),atol=1e-5)
    mt=np.char.startswith(sy.astype(str),'MT-');assert mt.any()
    pct=np.asarray(x[:,mt].sum(1)).ravel()/np.maximum(lib,1)*100
    mad=lambda v:np.median(np.abs(v-np.median(v)))
    lg=np.log10(ng+1);lc=np.log10(lib+1)
    ming=max(200,int(np.floor(10**(np.median(lg)-3*mad(lg))-1)))
    minc=max(300,int(np.floor(10**(np.median(lc)-3*mad(lc))-1)))
    maxmt=min(20.,max(5.,float(np.median(pct)+3*mad(pct))))
    obs=src.obs.copy();obs['barcode']=src.obs_names;obs.index=pd.Index([sm+'_'+b for b in src.obs_names],name='cell_id')
    obs['official_processed_member']=[b in accepted for b in src.obs_names]
    obs['n_counts']=lib;obs['n_genes']=ng;obs['pct_mt']=pct
    obs['official_singlet']=obs.scDblFinder_class.astype(str).str.lower().eq('singlet')
    obs['adaptive_qc_pass']=(ng>=ming)&(lib>=minc)&(pct<=maxmt)
    obs['analysis_keep']=obs.official_processed_member&obs.official_singlet&obs.adaptive_qc_pass
    obs['sample']=sm;obs['patient_id']=row['participant_id'];obs['dataset']=row['scpca_project_id'];obs['site']=row['tissue_location']
    obs.to_csv(O/'tables'/f'{sm}_all_cell_qc.tsv.gz',sep='\t')
    var=pd.DataFrame({'symbol':sy},index=src.var_names.copy())
    keep=obs.analysis_keep.to_numpy();a=ad.AnnData(x[keep].astype(np.float32),obs=obs.loc[keep].copy(),var=var)
    assert a.n_obs>0
    a.layers['counts']=a.X.copy()
    sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
    # Gene-level support computed from original counts, not normalized or integrated values.
    xc=a.layers['counts'].tocsc();den=a.obs.n_counts.to_numpy();sy=a.var.symbol.to_numpy()
    evidence=pd.DataFrame(index=a.obs_names)
    for name,genes in {**panels,**customer}.items():
        vals=[np.asarray(xc[:,np.flatnonzero(sy==g)].sum(1)).ravel() for g in genes if (sy==g).any()]
        vals=np.column_stack(vals) if vals else np.zeros((a.n_obs,0))
        evidence['det_'+name]=(vals>=1).sum(1);evidence['pct_'+name]=vals.sum(1)/np.maximum(den,1)*100
        if name in customer:
            panelrows.append({'sample':sm,'patient_id':row['participant_id'],'site':row['tissue_location'],'panel':name,'all_QC_cells':a.n_obs,'covered_genes':vals.shape[1],'requested_genes':len(genes),'raw_supported2_all_lineages':int(((vals>=1).sum(1)>=2).sum()),'interpretation':'all-lineage QC expression; not CAF counts'})
    if a.n_obs>=30:
        sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=3000)
        excluded=a.var.symbol.str.match(r'^(MT-|RPL|RPS)')|a.var.symbol.isin(['XIST','RPS4Y1','DDX3Y','EIF1AY'])
        a.var.loc[excluded,'highly_variable']=False
        b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10)
        npc=min(40,b.n_obs-1,b.n_vars-1);sc.tl.pca(b,n_comps=npc,svd_solver='arpack',random_state=270231)
        a.obsm['X_pca']=b.obsm['X_pca'];del b
        sc.pp.neighbors(a,n_neighbors=min(15,a.n_obs-1),n_pcs=npc,random_state=270231)
        for res in [.4,.8,1.2]:sc.tl.leiden(a,resolution=res,key_added='leiden_'+str(res),flavor='igraph',n_iterations=2,directed=False,random_state=270231)
        a.obs['cluster']=a.obs['leiden_0.8'].astype(str)
    else:a.obs['cluster']='0'
    lineages=['T_NK','B_Plasma','Myeloid','DC','Endo','Epi','Fib','Osteo','Mural','OC','Myogenic','Mast','Erythroid']
    labels={}
    for cl in a.obs.cluster.unique():
        mask=a.obs.cluster.eq(cl).to_numpy();r={'sample':sm,'cluster':cl,'cells':int(mask.sum())};active=[]
        for name in lineages:
            frac=float((evidence.loc[mask,'det_'+name]>=2).mean());p=float(evidence.loc[mask,'pct_'+name].median())
            r['frac2_'+name]=frac;r['median_pct_'+name]=p
            if frac>=.2:active.append((name,p))
        active.sort(key=lambda t:t[1],reverse=True)
        label=active[0][0] if active else 'Unresolved'
        if len(active)>1 and active[1][1]>=.8*active[0][1]:label='Mixed'
        if label in ['Fib','Osteo','Mural'] and sum(n in ['Fib','Osteo','Mural'] for n,v in active)>1:label='MesAmb'
        label={'Fib':'FibLike','Osteo':'OsteoLike'}.get(label,label);labels[cl]=label;r['lineage_candidate']=label;clusterrows.append(r)
    a.obs['lineage_candidate']=a.obs.cluster.map(labels).astype(str)
    a.obs['cnv_status']='not_run';a.obs['ambient_status']='not_run';a.uns['identity_status']='provisional_not_CAF'
    for col in evidence:a.obs[col]=evidence[col]
    a.write_h5ad(O/'objects'/f'{sm}_qc.h5ad',compression='gzip')
    a.obs.to_csv(O/'tables'/f'{sm}_cell_evidence.tsv.gz',sep='\t',index_label='cell_id')
    inp=O/'inputs'/sm;inp.mkdir(exist_ok=True)
    a.obs.reset_index(names='cell_id').to_csv(inp/'cells.tsv',sep='\t',index=False)
    pd.DataFrame({'gene_id':a.var_names,'ensembl_id':a.var_names,'symbol':a.var.symbol.to_numpy()}).to_csv(inp/'genes.tsv',sep='\t',index=False)
    with gzip.open(inp/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,a.layers['counts'].T,field='integer')
    (inp/'COMPLETE').write_text('COUNTS_EXPORT_COMPLETE_NOT_CNV_READY\n')
    summ.append({'sample':sm,'project':row['scpca_project_id'],'patient_id':row['participant_id'],'site':row['tissue_location'],'filtered_cells':len(obs),'official_processed':len(accepted),'retained':a.n_obs,'min_genes':ming,'min_counts':minc,'max_mt':maxmt,'clusters':a.obs.cluster.nunique(),'raw_Fib2':int((evidence.det_Fib>=2).sum()),'strict_CAF':None,'status':'QC_COMPLETE_IDENTITY_PENDING'})
    pd.DataFrame(summ).to_csv(O/'tables/qc_summary.tsv',sep='\t',index=False)
    pd.DataFrame(panelrows).to_csv(O/'tables/customer_panels_all_lineages.tsv',sep='\t',index=False)
    pd.DataFrame(clusterrows).to_csv(O/'tables/cluster_lineage_review.tsv',sep='\t',index=False)
    state['completed_samples'].append(sm);dump(O/'STATUS.json',state);print('SAMPLE_QC_COMPLETE',json.dumps(summ[-1]),flush=True)
    del a,src,x,xc,evidence
 state.update(status='QC_COMPLETE_LINEAGE_AMBIENT_CNV_REVIEW_PENDING',finished_utc=now(),samples=len(summ),retained_cells=sum(r['retained'] for r in summ))
 dump(O/'STATUS.json',state);(O/'QC_COMPLETE').write_text('QC_COMPLETE\n');print('QC_COMPLETE',flush=True)
except Exception:
 state.update(status='QC_FAILED_REVIEW_REQUIRED',error=traceback.format_exc());dump(O/'STATUS.json',state);raise
