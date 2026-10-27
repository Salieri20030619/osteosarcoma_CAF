from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import sys,json,gzip,hashlib,subprocess,datetime,importlib.util,concurrent.futures,shutil
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.io import mmwrite,mmread
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase88_expansion_identity_20261014';Q=O/'preprocess'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
def dump(p,o):p.write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def symcount(x,symbols,g):return np.asarray(x[:,np.flatnonzero(symbols==g)].sum(1)).ravel()
def main():
    ncc=P/'phase87_NCC_preprocessing_20261014'
    assert json.loads((ncc/'JOB_FINISHED.json').read_text(encoding='utf-8'))['success']
    manifest=pd.read_csv(ncc/'manifest.tsv',sep='\t').fillna('');rows=manifest.to_dict('records');sources=[]
    for r in rows:
        sm=r['sample'];src=ncc/'accepted'/sm/'analysis_ready.h5ad'
        old=json.loads((ncc/'accepted'/sm/'ACCEPTANCE.json').read_text(encoding='utf-8'))
        digest=sha(src);assert digest==old['accepted_sha256']
        sources.append({'sample':sm,'path':str(src),'sha256':digest})
        d=Q/'accepted'/sm;d.mkdir(exist_ok=True)
        if not (d/'analysis_ready.h5ad').exists():(d/'analysis_ready.h5ad').symlink_to(src)
    cx=P/'expansion_CELLxGENE_20261014';info=json.loads((cx/'EXTRACTION_COMPLETE.json').read_text(encoding='utf-8'))
    src=cx/'objects/osteosarcoma_counts.h5ad';assert sha(src)==info['counts_object_sha256']
    sources.append({'sample':'CELLxGENE_OS','path':str(src),'sha256':info['counts_object_sha256']})
    a=ad.read_h5ad(src);assert a.n_obs==17465
    mapping={'EXT_OS_1-1':'CXG__EXT_OS_1_L1','EXT_OS_1-2':'CXG__EXT_OS_1_L2','EXT_OS_2':'CXG__EXT_OS_2'}
    assert set(a.obs.Sample_ID.astype(str))==set(mapping)
    for source_library,sm in mapping.items():
        b=a[a.obs.Sample_ID.astype(str)==source_library].copy()
        assert b.obs.donor_id.nunique()==1 and b.obs.Cellgen_ID.nunique()==1
        donor=str(b.obs.donor_id.iloc[0])
        # Preserve source cell IDs verbatim as barcodes and all author metadata columns.
        d=Q/'inputs'/sm;d.mkdir(exist_ok=True);b.write_h5ad(d/'counts.h5ad',compression='gzip')
        r=dict(sample=sm,short=donor,patient_id='CXG__'+donor,matrix_dir=str(d),dataset='Whitfield2026_CELLxGENE',site='limb',group='unassigned',disease_stage_raw='osteosarcoma',treatment_raw='unknown',library_mapping_status='donor-level object; capture mapping unresolved',source_evidence='CELLxGENE dataset_version 36d1ebc5-9ad3-4b6d-b321-eb1ceaa550cf',technology='10x_multiome_single_nucleus',patient_independence_status='donor_labels_only_cross_study_independence_unresolved',source_type='human_tumor',clinical_pairing='unknown',cross_cohort_patient_pooling='not_allowed',Gender='unknown',Age='unknown')
        rows.append(r)
        r.update(short=source_library,source_library=source_library,source_Cellgen_ID=str(b.obs.Cellgen_ID.iloc[0]),library_mapping_status='one observed author Sample_ID and Cellgen_ID per processing unit',clinical_pairing='two_libraries_same_donor' if donor=='EXT_OS_1' else 'unpaired',author_patient_status=';'.join(sorted(b.obs.Patient_Status.astype(str).unique())),author_histological_subtype=';'.join(sorted(b.obs.Subtype.astype(str).unique())))
    del a,b
    pd.DataFrame(rows).to_csv(Q/'manifest.tsv',sep='\t',index=False,encoding='utf-8')
    dump(O/'resources/SOURCE_HASHES.json',sources)
    spec=importlib.util.spec_from_file_location('preprocess_worker',Q/'code/process_ncc.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.O=Q
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        results=list(ex.map(mod.worker,list(mapping.values())))
    dump(Q/'CELLxGENE_QC_RESULTS.json',results)
    assert all(r['status']=='preprocessing_complete_identity_pending' for r in results),results
    # Recompute descriptive cluster markers on final post-doublet cell sets for all samples.
    for r in rows:
        sm=r['sample'];a=ad.read_h5ad(Q/'accepted'/sm/'analysis_ready.h5ad');x=a.layers['counts'].tocsr()
        assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
        assert np.array_equal(np.asarray(x.sum(1)).ravel(),a.obs.n_counts.to_numpy())
        d=Q/'results'/sm;d.mkdir(exist_ok=True)
        target=d/'annotated.h5ad'
        if target.exists():
            # CELLxGENE raw QC object exists here: preserve it before replacing with the final accepted view.
            target.rename(d/'before_alternate_doublet_review.h5ad')
        target.symlink_to(Q/'accepted'/sm/'analysis_ready.h5ad')
        if (d/'complete.json').exists():shutil.copy2(d/'complete.json',d/'QC_complete_before_identity.json')
        dump(d/'complete.json',{'retained':a.n_obs,'identity_input':'final accepted post-doublet cells'})
        norm=sp.diags(1e4/np.asarray(x.sum(1)).ravel()).dot(x).tocsr();sums=np.asarray(norm.sum(0)).ravel();dets=np.asarray((x>0).sum(0)).ravel()
        markers=[]
        for cl in a.obs.cluster.astype(str).unique():
            mask=(a.obs.cluster.astype(str)==cl).to_numpy();n=int(mask.sum());nr=a.n_obs-n
            if nr==0:continue
            smm=np.asarray(norm[mask].sum(0)).ravel();dt=np.asarray((x[mask]>0).sum(0)).ravel()
            mean=smm/n;other=np.maximum((sums-smm)/nr,0);pct=dt/n;pct_other=(dets-dt)/nr
            lfc=np.log2((mean+.01)/(other+.01));score=lfc*np.maximum(pct-pct_other,0)
            ix=np.flatnonzero((pct>=.1)&(lfc>0)&(pct>pct_other));top=ix[np.argsort(-score[ix],kind='stable')[:50]]
            markers.extend({'group':cl,'names':str(a.var_names[j]),'symbol':str(a.var.symbol.iloc[j]),'scores':float(score[j]),'log2FC':float(lfc[j]),'pct_in':float(pct[j]),'pct_out':float(pct_other[j])} for j in top)
        pd.DataFrame(markers).to_csv(d/'cluster_markers.tsv',sep='\t',index=False,encoding='utf-8')
        print('MARKERS_READY',r['sample'],a.n_obs,flush=True);del a,x,norm
    subprocess.run([PY,'-u',str(O/'code/review_lineages.py')],check=True)
    pack=pd.read_csv(Q/'review/cluster_review_packets.tsv',sep='\t',dtype={'cluster':str})
    go=pd.read_csv(O/'resources/infercnv_gene_order_GRCh38_98.tsv',sep='\t',header=None,names=['id','chr','start','end'],dtype={'chr':str}).set_index('id')
    summary=[];checks=[]
    for r in rows:
        sm=r['sample'];a=ad.read_h5ad(Q/'accepted'/sm/'analysis_ready.h5ad');x=a.layers['counts'].tocsr();symbols=a.var.symbol.astype(str).to_numpy()
        o=pd.read_csv(Q/'review'/(sm+'_cell_evidence.tsv.gz'),sep='\t',dtype={'cluster':str}).set_index('cell_id',drop=False).loc[a.obs_names]
        assert o.index.tolist()==a.obs_names.tolist()
        pp=pack[pack['sample']==sm].set_index('cluster');o['CTA_context']=o.cluster.map(pp.CTA_context).astype(bool)
        o['pctUMI_OC_specific']=sum(symcount(x,symbols,g) for g in ['OCSTAMP','DCSTAMP','ATP6V0D2','SIGLEC15','SUCNR1','CALCR'])/o.n_counts.to_numpy()*100
        o['pctUMI_competing_mes']=o.pctUMI_Fib+o.pctUMI_Osteo+o.pctUMI_Mural
        ds=o.scDblFinder_score if 'scDblFinder_score' in o and o.scDblFinder_score.notna().all() else o.doublet_score
        base=(o.n_genes>=500)&(o.n_counts>=750)&ds.notna()&(ds<=ds.quantile(.9))&(o.pctUMI_competing_mes<=.2)&((o.pctUMI_Epi+o.pctUMI_Myogenic)<=.2)&(o.pctUMI_OC_specific<=.1)&~o.CTA_context&~o.doublet_positive_any.astype(bool)
        o['ref_immune']=False
        for typ in ['T_NK','B_Plasma','Myeloid','DC']:
            o['ref_immune']|=(o.lineage_phase3==typ)&(o['det_'+typ+'_review']>=3)&(o['pctUMI_'+typ]>=.5)&base&(o.pctUMI_Endo<=.2)
        ipct=sum(o['pctUMI_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC'])
        o['ref_endo']=(o.lineage_phase3=='Endo')&(o.det_Endo_review>=3)&(o.pctUMI_Endo>=.5)&base&(ipct<=.2)
        o['cnv_role']='target';rng=np.random.default_rng(270231)
        for typ,cap in [('immune',400),('endo',160)]:
            ids=np.flatnonzero(o['ref_'+typ].to_numpy());rng.shuffle(ids);ids=ids[:cap]
            if len(ids)<20:continue
            nh=max(10,int(len(ids)*.2));o.iloc[ids[:nh],o.columns.get_loc('cnv_role')]='holdout_'+typ;o.iloc[ids[nh:],o.columns.get_loc('cnv_role')]='ref_'+typ
        ref=o.cnv_role.str.startswith('ref_');hold=o.cnv_role.str.startswith('holdout_')
        if hold.sum()<20 and ref.sum()+hold.sum()>=50:
            ix=np.flatnonzero(ref.to_numpy());rng.shuffle(ix)
            for j in ix[:20-int(hold.sum())]:o.iloc[j,o.columns.get_loc('cnv_role')]=str(o.iloc[j].cnv_role).replace('ref_','holdout_')
        ref=o.cnv_role.str.startswith('ref_');hold=o.cnv_role.str.startswith('holdout_')
        ids=a.var_names.astype(str).tolist();assert all(g.startswith('ENSG') for g in ids) and len(ids)==len(set(ids))
        expressed=np.asarray((x>0).sum(0)).ravel()>=3
        autocov={str(ch):sum(bool(e) and g in go.index and go.at[g,'chr']==str(ch) for g,e in zip(ids,expressed)) for ch in range(1,23)}
        eligible=bool(ref.sum()>=30 and hold.sum()>=20 and min(autocov.values())>=30)
        assert not (ref&hold).any();assert (o.loc[ref|hold,'pctUMI_competing_mes']<=.2).all()
        assert not o.loc[ref|hold,'lineage_phase3'].isin(['FibLike','OsteoLike','Mural','MesAmb','Mixed','Unresolved','OC_like','Myogenic_pending']).any()
        o['tool_id']=[sm+'_c'+str(i).zfill(6) for i in range(len(o))]
        d=O/'inputs'/sm;d.mkdir(exist_ok=True)
        o.to_csv(d/'cells.tsv',sep='\t',index=False,encoding='utf-8')
        pd.DataFrame({'ensembl_id':ids,'gene_id':ids,'symbol':symbols,'source_gene_id':ids}).to_csv(d/'genes.tsv',sep='\t',index=False,encoding='utf-8')
        with gzip.open(d/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
        # Exact matrix roundtrip and identifiers before either CNV or ambient computation.
        m=mmread(str(d/'counts.mtx.gz')).T.tocsr();assert m.shape==x.shape and (m!=x).nnz==0
        pd.DataFrame({'tool_id':o.tool_id,'annotation':np.where(ref,o.cnv_role,'C'+o.cluster.astype(str)+'_'+o.lineage_phase3)}).to_csv(d/'annotations.tsv',sep='\t',index=False,header=False,encoding='utf-8')
        (d/'INPUT_COMPLETE').write_text('CNV_INPUT_EXPORT_COMPLETE\n',encoding='utf-8');(d/'COMPLETE').write_text('COUNTS_AND_CELL_ORDER_VERIFIED\n',encoding='utf-8')
        dump(d/'gene_coverage.json',autocov)
        z=dict(sample=sm,dataset=r['dataset'],patient_id=r['patient_id'],site=r['site'],technology=r['technology'],cells=len(o),FibLike=int(o.lineage_phase3.eq('FibLike').sum()),MesAmb=int(o.lineage_phase3.eq('MesAmb').sum()),ref_cells=int(ref.sum()),holdout_cells=int(hold.sum()),cnv_eligible=eligible,min_expressed_genes_per_autosome=min(autocov.values()),ambient_eligible=o.lineage_phase3.nunique()>=2,reason='PASS' if eligible else 'insufficient_reference_holdout_or_autosome_coverage')
        summary.append(z);checks.append({'sample':sm,'counts_exact_roundtrip':True,'cell_order_verified':True,'genes_unique':True,'source':str(Q/'accepted'/sm/'analysis_ready.h5ad')})
        pd.DataFrame(summary).to_csv(O/'review/input_summary.tsv',sep='\t',index=False,encoding='utf-8');print('INPUT_GATE',json.dumps(z),flush=True)
        del a,x,m,o
    dump(O/'resources/EXPORT_ACCEPTANCE.json',checks)
    (O/'INPUTS_READY').write_text('ALL12_IDENTITY_INPUTS_READY\n',encoding='utf-8')
    print('PREPARE_COMPLETE',flush=True)
if __name__=='__main__':main()
