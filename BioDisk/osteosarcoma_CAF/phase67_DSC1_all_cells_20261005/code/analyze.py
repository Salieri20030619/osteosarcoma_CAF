from pathlib import Path
import os,sys,json,datetime,hashlib,traceback,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=B/'phase67_DSC1_all_cells_20261005'
for s in ['code','tables','resources','QA','delivery','logs']:(O/s).mkdir(parents=True,exist_ok=True)
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
 annpath=B/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz'
 cols=['original_cell_id','sample','dataset','original_donor_record','site','technology','full_library','review_expression_lineage','review_lineage_code','lineage_confidence','malignancy_assessment','candidate_tier','cluster_id']
 ann=pd.read_csv(annpath,sep='\t',usecols=cols,keep_default_na=False).set_index('original_cell_id')
 assert len(ann)==876648 and ann.index.is_unique
 reg=pd.read_csv(B/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
 inv=pd.read_csv(B/'phase53_de_novo_reannotation_20261002/resources/SAMPLE_INPUT_INVENTORY.tsv',sep='\t').set_index('sample')
 cur=json.loads((B/'CAF_CURRENT_CELL_EVIDENCE.json').read_text())
 cm=pd.read_csv(cur['cell_table'],sep='\t',low_memory=False).set_index('cell_id')
 assert len(cm)==8174 and cm.index.is_unique and cm.index.isin(ann.index).all()
 prior=pd.read_csv(B/'phase66_doublet_repair_DSC1_20261005/tables/desmosome_cell_counts.tsv.gz',sep='\t')
 prior=prior[prior.matrix.eq('raw')].set_index('cell_id')
 protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'request':'Examine SR1004 target DSC1 expression across all cells following reported high wet-lab signal.','gene':'DSC1','feature':'ENSG00000134765','scope':'All 876648 source-called input cells, 107 samples, 9 sources; includes cells failing provider QC; no new cell filtering. This is not the 8174 candidate subset and not a post-QC singlet atlas.','matrix':'Per-sample original full_called_counts.h5ad X; original count layer, not normalized/integrated/HVG matrix.','mapping':'Unique exact Ensembl base ID ENSG00000134765 or exact DSC1 symbol; no family-member substitution; absence is missing, not zero.','metrics':'Raw count, raw>=1/2/5 detection, total-library CP10k, sample pseudobulk CP10k, within-source donor-site and donor equal-site descriptive summaries. Thresholds 2/5 are descriptive, not biological high-expression classifications.','statistics':'Descriptive only; cells not independent biological replicates; source-scoped donor records not pooled as globally independent patients. No P values.','annotation':str(annpath),'annotation_sha256':sha(annpath),'candidate_manifest':cur,'comparability':'Customer antibody catalog, assay, sample/cell compartment still unconfirmed; RNA measurements do not establish protein level.','doublets':'Do not inherit missing phase55 doublet fields as negatives; native caller/QC metadata retained for positive cells; phase66 current flags used only inside the validated 8174 subset.'}
 dump(O/'resources/PROTOCOL.json',protocol)
 parts=[];coverage=[];positives=[];maxliberr=0
 for n,row in enumerate(reg.itertuples(index=False),1):
  f=B/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad'
  stat=f.stat()
  a=ad.read_h5ad(f,backed='r')
  ids=a.obs_names.astype(str);assert ids.is_unique and ids.isin(ann.index).all()
  t=ann.loc[ids].copy();t.index.name='original_cell_id'
  assert t['sample'].eq(row.sample).all() and t.dataset.eq(row.dataset).all()
  assert len(t)==int(inv.loc[row.sample,'cells'])
  features=a.var_names.astype(str)
  hits=np.flatnonzero(np.array([v.split('.')[0]=='ENSG00000134765' or v=='DSC1' for v in features]))
  assert len(hits)<=1,(row.sample,features[hits].tolist())
  j=int(hits[0]) if len(hits) else None
  if j is not None:
   for col in ['symbol','gene_symbol']:
    if col in a.var:
     sym=str(a.var.iloc[j][col]);assert sym in ['DSC1','','nan'],(row.sample,sym)
  count=np.full(len(t),np.nan);lib=np.zeros(len(t));nz=np.zeros(len(t),dtype=int)
  with h5py.File(f,'r') as h:
   x=h['X'];enc=x.attrs.get('encoding-type','')
   if enc=='csr_matrix':
    ip=x['indptr'][:]
    for lo in range(0,len(t),4096):
     hi=min(lo+4096,len(t));st,en=int(ip[lo]),int(ip[hi])
     block=sparse.csr_matrix((x['data'][st:en].astype(np.float64),x['indices'][st:en],ip[lo:hi+1]-st),shape=(hi-lo,a.n_vars))
     assert np.isfinite(block.data).all() and (block.data>=0).all() and np.allclose(block.data,np.round(block.data),rtol=0,atol=1e-6)
     lib[lo:hi]=np.asarray(block.sum(axis=1)).ravel();nz[lo:hi]=np.asarray((block>0).sum(axis=1)).ravel()
     if j is not None:count[lo:hi]=block[:,j].toarray().ravel()
   else:
    raise RuntimeError('Unexpected matrix encoding '+str(enc))
  maxliberr=max(maxliberr,float(np.max(np.abs(lib-pd.to_numeric(t.full_library).to_numpy()))))
  assert np.allclose(lib,pd.to_numeric(t.full_library),rtol=1e-6,atol=1e-6)
  t=t.drop(columns='full_library');t['full_library']=lib;t['detected_features']=nz
  t['DSC1_assayed']=j is not None;t['DSC1_raw']=count
  t['DSC1_CP10k']=np.divide(count*1e4,lib,out=np.full(len(t),np.nan),where=lib>0)
  t['original_candidate8174']=t.index.isin(cm.index)
  t['current_candidate_doublet_positive']=t.index.map(cm.doublet_positive_any.to_dict())
  t['candidate_pool']=np.where(t.original_candidate8174,'inside_8174','outside_8174')
  pos=t.loc[t.DSC1_raw.ge(1)].copy()
  qcols=[c for c in a.obs if any(k in c.lower() for k in ['doublet','dblfinder','scrublet','miqc','scpca_filter','subsets_mito_percent','prob_compromised','provider_consensus_celltype'])]
  for c in qcols:pos['native_'+c]=a.obs.loc[pos.index,c].astype(str)
  positives.append(pos.reset_index())
  parts.append(t.reset_index())
  coverage.append({'sample':row.sample,'dataset':row.dataset,'cells':len(t),'DSC1_assayed':j is not None,'feature':features[j] if j is not None else '', 'positive_cells':int((count>=1).sum()),'raw_total':float(np.nansum(count)) if j is not None else None,'input_path':str(f),'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns})
  a.file.close()
  print(f'{n}/107 {row.sample} cells={len(t)} positive={int((count>=1).sum())}',flush=True)
 d=pd.concat(parts,ignore_index=True)
 assert len(d)==len(ann) and d.original_cell_id.is_unique and set(d.original_cell_id)==set(ann.index)
 assert d['sample'].nunique()==107 and d.dataset.nunique()==9
 cand=d[d.original_candidate8174].set_index('original_cell_id').loc[prior.index]
 assert np.allclose(cand.DSC1_raw,prior.DSC1,equal_nan=True)
 assert int(cand.DSC1_assayed.sum())==8170 and int(cand.DSC1_raw.ge(1).sum())==3
 d.to_csv(O/'tables/DSC1_all_cells.tsv.gz',sep='\t',index=False)
 pd.concat(positives,ignore_index=True).sort_values('DSC1_raw',ascending=False).to_csv(O/'tables/DSC1_positive_cells_with_native_QC.tsv.gz',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(O/'tables/sample_feature_coverage.tsv',sep='\t',index=False)
 def summ(g):
  v=g[g.DSC1_assayed];p=v[v.DSC1_raw.ge(1)]
  return {'cells':len(g),'assayed_cells':len(v),'missing_cells':len(g)-len(v),'positive_cells':len(p),'positive_pct':100*len(p)/len(v) if len(v) else np.nan,'ge2_cells':int(v.DSC1_raw.ge(2).sum()),'ge5_cells':int(v.DSC1_raw.ge(5).sum()),'raw_total':float(v.DSC1_raw.sum()) if len(v) else np.nan,'mean_raw_all_assayed':float(v.DSC1_raw.mean()),'max_raw':float(v.DSC1_raw.max()),'mean_CP10k_all_assayed':float(v.DSC1_CP10k.mean()),'pseudobulk_CP10k':float(v.DSC1_raw.sum()/v.full_library.sum()*1e4) if len(v) and v.full_library.sum()>0 else np.nan,'positive_median_raw':float(p.DSC1_raw.median()),'positive_median_CP10k':float(p.DSC1_CP10k.median()),'samples':g['sample'].nunique(),'positive_samples':p['sample'].nunique(),'source_scoped_donor_records':g.original_donor_record.nunique(),'positive_donor_records':p.original_donor_record.nunique(),'zero_library_cells':int(v.full_library.eq(0).sum())}
 total=summ(d)
 for name,keys in [('by_source',['dataset']),('by_sample',['dataset','sample','original_donor_record','site','technology']),('by_lineage',['review_lineage_code','review_expression_lineage']),('by_source_lineage',['dataset','review_lineage_code','review_expression_lineage']),('by_cluster',['dataset','cluster_id','review_lineage_code']),('by_candidate_pool',['candidate_pool']),('by_fibro_candidate_tier',['candidate_tier']),('by_donor_site',['dataset','original_donor_record','site'])]:
  out=[]
  for k,g in d.groupby(keys,dropna=False,sort=True):
   if not isinstance(k,tuple):k=(k,)
   out.append(dict(zip(keys,k),**summ(g)))
  pd.DataFrame(out).to_csv(O/'tables'/('DSC1_'+name+'.tsv'),sep='\t',index=False)
 ds=pd.read_csv(O/'tables/DSC1_by_donor_site.tsv',sep='\t')
 donor=ds.groupby(['dataset','original_donor_record'],as_index=False).agg(sites=('site','size'),assayed_cells=('assayed_cells','sum'),positive_cells=('positive_cells','sum'),site_equal_positive_pct=('positive_pct','mean'),site_equal_pseudobulk_CP10k=('pseudobulk_CP10k','mean'))
 donor.to_csv(O/'tables/DSC1_by_donor_equal_site.tsv',sep='\t',index=False)
 # Independent aggregate reconciliation against saved sample and lineage tables.
 for name in ['by_source','by_sample','by_lineage','by_cluster','by_candidate_pool']:
  z=pd.read_csv(O/'tables'/('DSC1_'+name+'.tsv'),sep='\t')
  assert int(z.cells.sum())==len(d) and int(z.positive_cells.sum())==total['positive_cells']
  assert np.isclose(z.raw_total.sum(),total['raw_total'])
 positive=d[d.DSC1_raw.ge(1)]
 raw_dist=positive.DSC1_raw.value_counts().sort_index()
 raw_dist.rename_axis('DSC1_raw').reset_index(name='cells').to_csv(O/'tables/DSC1_positive_count_distribution.tsv',sep='\t',index=False)
 result={'overall':total,'candidate8174':summ(d[d.original_candidate8174]),'outside_candidate8174':summ(d[~d.original_candidate8174]),'max_library_reproduction_error':maxliberr,'candidate_phase66_raw_counts_reproduced':True,'input_scope':'all source-called cells, includes provider QC failures','inference':'descriptive_only','protein_level_evaluated':False}
 dump(O/'resources/SUMMARY.json',result)
 dump(O/'QA/ACCEPTANCE.json',{'passed':True,'cells':len(d),'samples':107,'sources':9,'all_original_ids_matched':True,'counts_nonnegative_integer':True,'libraries_recomputed':True,'candidate_phase66_reproduced':True,'aggregate_totals_reconciled':True,'missing_not_zero':True,'no_cell_level_inference':True})
 (O/'ANALYSIS_COMPLETE').write_text('SUCCESS\n')
 print('RESULT '+json.dumps(result,ensure_ascii=False),flush=True)
if __name__=='__main__':
 status=1
 try:run();status=0
 except Exception:traceback.print_exc()
 finally:
  (O/'logs/analysis.exit_code').write_text(str(status)+'\n')
  dump(O/'QA/EXECUTION_STATUS.json',{'exit_code':status,'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 sys.exit(status)
