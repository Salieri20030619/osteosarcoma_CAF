from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,importlib.util,urllib.request,hashlib,sys,traceback,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase23_loose_first_20260921'
for d in ['code','inputs','raw','tables','objects','figures','logs','delivery','resources']:(O/d).mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def truth(s):return s.fillna(False).astype(str).str.lower().isin(['true','1','1.0'])
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str))
def status(s,**kw):dump(O/'STATUS.json',dict(status=s,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw));print(s,kw,flush=True)
def module(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'user_direction':'start with the loosest provider-called complete cell input, then add evidence; no project thresholds at entry','scope':'52 source sample records; original analyses immutable','entry':'all provider-called cells: legacy source count matrices and ScPCA filtered matrices, NOT raw empty-droplet matrices','excluded_at_entry':'only invalid matrices/identifiers; zero-count cells retained in ledger but technically unembeddable','not_entry_gates':['UMI','gene count','mitochondrial fraction','doublet score quantile','official processed membership','Fib/ECM marker count','competing lineage percentage','CNV high or low','normal reference count','heldout reference acceptance','ambient correction','cluster CNV fraction','original FibLike label'],'official_filters':'software-specific, applied only when running that software; no universal official CAF threshold exists','ScPCA':'official processed membership, singlet/doublet annotations, consensus lineage and official inferCNV stored as evidence; all filtered cells remain at entry','legacy_CNV':'previous calls preserved as historical evidence with original reference selection; project inferCNV high/low calls never called official defaults','new_model':'all input cells; expression-only unsupervised joint model; old C7/C9 used only for overlap after clustering','numerical_choices':'logCP10k full-library denominator; consensus batch HVG2000; PCA30; Harmony study/technology; neighbors15; Leiden resolution0.5; these are modeling choices, NOT claimed official acceptance thresholds','statistics':'independent patient unit; no cells as replicates; exploratory post-selection inference with family-specific FDR','sources':{'ScPCA':'https://scpca.readthedocs.io/en/stable/processing_information.html','CopyKAT':'https://github.com/navinlabcode/copykat/blob/master/man/copykat.Rd','inferCNV':'https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV'}}
dump(O/'PROTOCOL.json',protocol)
old=read(P/'phase17_all_cohort_expanded_20260921/tables/all_cell_layer_evidence.tsv.gz').set_index('cell_id')
oldmeta=read(P/'phase10_batch_gates_20260919/tables/all_cell_unified_evidence.tsv.gz').set_index('cell_id')
base=read(P/'phase17_all_cohort_expanded_20260921/tables/sample_layer_counts.tsv').set_index('sample')
legacy=read(P/'phase22_no_UMI_sensitivity_20260921/legacy_UMI_audit.tsv').set_index('sample')
m0=read(P/'OS_CAF_phase2/code/manifest.tsv').set_index('sample');m5=read(P/'phase5_nonpdx_qc_20260916/manifest.tsv').set_index('sample');mg=read(P/'geo_preprocess_20260918/manifest.tsv').set_index('sample');ms=read(P/'phase16_scpca_supplement_20260920/manifest.tsv').set_index('scpca_sample_id')
panels=json.loads((P/'phase10_batch_gates_20260919/resources/panels.json').read_text());customer=json.loads((P/'phase16_scpca_supplement_20260920/resources/customer_panels.json').read_text());panels.update(customer)
dump(O/'resources/panels.json',panels)
sources=[]
for sm in base.index:
 if sm in m0.index:sources.append((sm,'tar',m0.loc[sm]))
 elif sm in m5.index:sources.append((sm,'mtx',m5.loc[sm]))
 elif sm in mg.index:sources.append((sm,'h5ad',mg.loc[sm]))
 else:assert sm in ms.index;sources.append((sm,'ScPCA',ms.loc[sm]))
assert len(sources)==52
reader0=module(P/'OS_CAF_phase2/code/phase2.py','old_reader0');reader5=module(P/'phase5_nonpdx_qc_20260916/code/phase5.py','old_reader5')
summary=[];hvgs=[];official=[];ledgers=[]
status('BUILDING_UNFILTERED_PROJECT_ENTRY',samples=52)
for sm,kind,r in sources:
 dest=O/'inputs'/f'{sm}.h5ad';done=O/'inputs'/f'{sm}.json'
 if dest.exists() and done.exists():
  summary.append(json.loads(done.read_text()));hvgs.append(read(O/'tables'/f'{sm}_hvg.tsv'));ledgers.append(read(O/'tables'/f'{sm}_entry.tsv.gz'));print('CACHED',sm,flush=True);continue
 if kind=='ScPCA':
  a=ad.read_h5ad(r.filtered);a.var['symbol']=a.var.gene_symbol.astype(str)
  q=read(P/'phase16_scpca_supplement_20260920/tables'/f'{sm}_all_cell_qc.tsv.gz').set_index('cell_id')
  a.obs_names=[sm+'_'+str(z) for z in a.obs_names]
  proc=ad.read_h5ad(r.processed,backed='r');po=proc.obs.copy();po.index=[sm+'_'+str(z) for z in po.index]
  spec={'sample':sm,'official_infercnv_status':str(proc.uns.get('infercnv_status','missing')),'official_reference_count':str(proc.uns.get('infercnv_num_reference_cells','missing')),'official_reference_celltypes':str(proc.uns.get('infercnv_reference_celltypes','missing'))}
  official.append(spec);dump(O/'inputs'/f'{sm}_official_metadata.json',spec)
  for c in ['consensus_celltype_annotation','singler_celltype_annotation','cellassign_celltype_annotation','scimilarity_celltype_annotation','is_infercnv_reference','infercnv_total_cnv']:
   if c in po:q['official_'+c]=po[c].reindex(q.index)
  proc.file.close();source_path=str(r.filtered)
 else:
  q=read(legacy.loc[sm,'QC_table']).set_index('cell_id')
  if kind=='tar':
   source_path=str(O/'raw'/r.source_url.rsplit('/',1)[1]);pp=Path(source_path)
   if not pp.exists():
    print('DOWNLOADING_SOURCE',sm,flush=True)
    for attempt in range(3):
     try:
      with urllib.request.urlopen(r.source_url,timeout=180) as inp,pp.with_suffix('.part').open('wb') as out:
       while True:
        chunk=inp.read(4*1024*1024)
        if not chunk:break
        out.write(chunk)
      pp.with_suffix('.part').replace(pp);break
     except Exception:
      if attempt==2:raise
   a,_=reader0.read_matrix(pp)
  elif kind=='mtx':source_path=str(r.matrix_dir);a,_=reader5.read_matrix(Path(r.matrix_dir))
  else:source_path=str(Path(r.matrix_dir)/'counts.h5ad');a=ad.read_h5ad(source_path)
  a.obs_names=[sm+'_'+str(z) for z in a.obs_names]
 assert a.obs_names.is_unique and set(a.obs_names)==set(q.index),(sm,a.shape,len(q))
 q=q.loc[a.obs_names].copy();x=a.X.tocsr().astype(np.float32);x.eliminate_zeros();den=np.asarray(x.sum(1)).ravel();ng=x.getnnz(1)
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data))
 assert np.allclose(den,q.n_counts) and np.array_equal(ng,q.n_genes),(sm,'raw count integrity')
 q['sample']=sm;q['dataset']=base.loc[sm,'dataset'];q['entry']=True;q['embeddable']=den>0;q['historical_QC_keep']=q.index.isin(old.index)
 for c in ['strict','supported','expanded','CNV_state','reference_QC_pass','ambient_checked','mural_flag','osteo_flag','doublet','original_lineage']:
  q['historical_'+c]=old[c].reindex(q.index)
 if kind!='ScPCA':
  for c in ['infer_call','cnv2_call','technology']:q['historical_'+c]=oldmeta[c].reindex(q.index)
 else:
  ev=read(P/'phase16_scpca_supplement_20260920/tables'/f'{sm}_identity_evidence.tsv.gz').set_index('cell_id')
  q['historical_infer_call']=ev.infer_call.reindex(q.index);q['historical_cnv2_call']=ev.copykat_call.reindex(q.index)
  q['historical_technology']='ScPCA_'+q.suspension_type.astype(str)
 if sm=='SCPCS000611':
  e=read(P/'phase22_no_UMI_sensitivity_20260921/CNV611/delivery/identity_evidence.tsv.gz').set_index('cell_id')
  q['supplement_infer_call']=e.infer_call.reindex(q.index);q['supplement_copykat_call']=e.copykat_call.reindex(q.index)
 q['technology']=q.get('technology_reported',q.get('historical_technology',pd.Series('unknown',index=q.index))).fillna('unknown').astype(str)
 q['study_technology']=q.dataset.astype(str)+'|'+q.technology
 q['donor_record']=q.dataset.astype(str)+'|'+q.patient_id.astype(str)
 sy=a.var.symbol.astype(str).to_numpy();xc=x.tocsc()
 for name,gs in panels.items():
  indices=[np.flatnonzero(sy==gene) for gene in gs];found=[ix for ix in indices if len(ix)]
  if found:
   v=np.column_stack([np.asarray(xc[:,ix].sum(1)).ravel() for ix in found]);q['det_'+name]=(v>0).sum(1);q['pct_'+name]=100*v.sum(1)/np.maximum(den,1);q['score_'+name]=np.log1p(v/np.maximum(den,1)[:,None]*1e4).mean(1)
  else:q['det_'+name]=0;q['pct_'+name]=0.;q['score_'+name]=np.nan
  q['coverage_'+name]=len(found)/len(gs)
 # Expression is recorded continuously; no panel count or coverage threshold gates entry.
 for c in q:
  if q[c].dtype==object or isinstance(q[c].dtype,pd.StringDtype) or str(q[c].dtype)=='boolean':q[c]=q[c].fillna('missing').astype(str).astype(object)
 q.index=pd.Index(q.index.to_numpy(dtype=object),dtype=object)
 ids=np.array([str(z).split('.')[0] for z in a.var_names])
 if q.dataset.iloc[0]=='GSE250015':
  gm=read(P/'geo_cnv_20260918/inputs'/sm/'genes.tsv').set_index('source_gene_id');ids=gm.ensembl_id.reindex(a.var_names).fillna('').astype(str).to_numpy()
 a=ad.AnnData(x,obs=q,var=pd.DataFrame({'symbol':sy,'ensembl_id':ids},index=pd.Index(a.var_names.to_numpy(dtype=object),dtype=object)))
 a.write_h5ad(dest,compression='gzip')
 w=a[q.embeddable].copy();sc.pp.normalize_total(w,target_sum=1e4);sc.pp.log1p(w);sc.pp.highly_variable_genes(w,flavor='seurat',n_top_genes=min(2000,w.n_vars))
 hv=w.var[['ensembl_id','symbol','dispersions_norm','highly_variable']].copy();hv['sample']=sm;hv['rank']=hv.dispersions_norm.rank(ascending=False,method='min');hv.to_csv(O/'tables'/f'{sm}_hvg.tsv',sep='\t',index=False);hvgs.append(hv)
 cols=['sample','dataset','patient_id','donor_record','site','n_counts','n_genes','pct_mt','entry','embeddable','historical_QC_keep']+[c for c in q if c.startswith(('historical_','official_','supplement_','score_','det_','pct_','coverage_'))]+['study_technology']
 cols=list(dict.fromkeys(cols));ledger=q[cols].reset_index(names='cell_id');ledger.to_csv(O/'tables'/f'{sm}_entry.tsv.gz',sep='\t',index=False);ledgers.append(ledger)
 rec={'sample':sm,'dataset':q.dataset.iloc[0],'source_kind':kind,'source_path':source_path,'entry_cells':len(q),'embeddable_cells':int((den>0).sum()),'old_QC_cells':int(q.historical_QC_keep.sum()),'restored_cells':int((~q.historical_QC_keep).sum()),'raw_count_integrity':True,'official_processed':int(truth(q.official_processed_member).sum()) if kind=='ScPCA' else None,'official_singlets':int(truth(q.official_singlet).sum()) if kind=='ScPCA' else None}
 dump(done,rec);summary.append(rec);print('INPUT_COMPLETE',rec,flush=True)
 del a,w,x,xc,q;gc.collect()
ss=pd.DataFrame(summary);ss.to_csv(O/'tables/input_recovery.tsv',sep='\t',index=False)
allobs=pd.concat(ledgers,ignore_index=True);assert allobs.cell_id.is_unique
allobs.to_csv(O/'tables/all_entry_evidence.tsv.gz',sep='\t',index=False)
hv=pd.concat(hvgs,ignore_index=True);hv=hv[hv.ensembl_id.astype(str).str.startswith('ENSG')].drop_duplicates(['sample','ensembl_id'])
cons=hv.groupby('ensembl_id').agg(samples=('sample','nunique'),HVG_samples=('highly_variable','sum'),median_rank=('rank','median'));cons=cons[cons.samples==52].sort_values(['HVG_samples','median_rank'],ascending=[False,True]);assert len(cons)>=2000
cons.to_csv(O/'tables/consensus_HVG_rank.tsv',sep='\t');(O/'resources/HVG2000.txt').write_text('\n'.join(cons.head(2000).index)+'\n')
dump(O/'delivery/input_QA.json',{'samples':len(ss),'entry_cells':int(ss.entry_cells.sum()),'old_QC_cells':int(ss.old_QC_cells.sum()),'restored_cells':int(ss.restored_cells.sum()),'no_project_filters_at_entry':True,'raw_count_integrity_all_pass':bool(ss.raw_count_integrity.all()),'shared_ensembl_genes':len(cons)})
(O/'INPUTS_COMPLETE').write_text('ALL_52_COMPLETE_INPUTS_VERIFIED\n');status('INPUTS_COMPLETE',cells=len(allobs),samples=len(ss))
