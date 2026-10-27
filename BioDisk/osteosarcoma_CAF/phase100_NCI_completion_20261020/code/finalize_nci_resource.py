import pathlib,json,hashlib,datetime,re,concurrent.futures
import numpy as np,pandas as pd,anndata as ad
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase100_NCI_completion_20261020';q=root/'phase95_NCI_expanded_reanalysis_20261016'
review=pd.read_csv(p/'tables/cluster_review_decisions.tsv',sep='\t',encoding='utf-8')
run=json.loads((q/'COMPUTATION_FINISHED.json').read_text(encoding='utf-8'))['samples']
genes=json.loads((root/'phase96_all14_literature_20261016/project_registry.json').read_text(encoding='utf-8'))
aliases=json.loads((root/'phase97_normal_background_20261020/resources/GENE_ALIAS_MAP.json').read_text(encoding='utf-8'))
def one(rr):
 sm=rr['sample'];f=q/'tables'/(sm+'_all_cell_evidence.tsv.gz')
 c=pd.read_csv(f,sep='\t',encoding='utf-8');assert c.cell_id.is_unique
 rv=review[review['sample'].eq(sm)].set_index('cluster')
 c['cluster_review_class']=c.cluster.map(rv.review_class);c['cluster_review_decision']=c.cluster.map(rv.decision)
 assert c.cluster_review_decision.notna().all()
 c['phase100_expanded_admission']=c.mechanical_expanded_candidate & c.cluster_review_decision.eq('retain_expanded_identity_unresolved')
 c['candidate_ambient_supported']=c.phase100_expanded_admission & c.decontX_lineage_winner.eq('Fibroblast') & c.decontX_detected_Fibroblast.ge(2)
 c['author_tumor_conflict']=c.Ann_Level0.eq('Tumor')
 c['author_CAF_agreement']=c.author_CAF
 c['cluster_identity_conflict']=c.cluster_review_class.isin(['O','C','S','X'])
 c['confirmed_pure_CAF']=False;c['independent_patient_verified']=False
 c=c.rename(columns={'final_expanded_admission':'phase95_final_expanded_admission','manual_cluster_review_pending':'phase95_manual_cluster_review_pending'})
 c['final_expanded_admission']=c.phase100_expanded_admission
 c['manual_cluster_review_pending']=False
 c['review_complete_in_scope']=True
 c.to_csv(p/'tables'/(sm+'_final_evidence.tsv.gz'),sep='\t',encoding='utf-8',index=False)
 admitted=c[c.phase100_expanded_admission].copy()
 summary=dict(sample=sm,cells=len(c),reference_pass=rr['reference_pass'],reference_diagnostics=json.dumps(rr.get('reference_checks',[])),native_status=rr['status'],mechanical_candidates=int(c.mechanical_expanded_candidate.sum()),expanded_candidates=len(admitted),ambient_supported=int(c.candidate_ambient_supported.sum()),author_CAF_agreement=int(admitted.author_CAF.sum()),author_tumor_conflict=int(admitted.author_tumor_conflict.sum()),author_tumor_excluded=int((~admitted.author_tumor_conflict).sum()),putative_group=c.putative_patient_group.iloc[0],patient_verified=False,site='Lung',assay='single_nucleus',completion='reference_not_released' if not rr['reference_pass'] else 'reviewed_candidates_exported' if len(admitted) else 'review_complete_no_retained_candidates')
 expr=[];coverage=[];obj=None
 if len(admitted):
  a=ad.read_h5ad(q/'inputs'/sm/'available_counts.h5ad',backed='r');ix=a.obs_names.get_indexer(admitted.cell_id);assert (ix>=0).all()
  out=a[ix,:].to_memory();a.file.close();assert out.obs_names.tolist()==admitted.cell_id.tolist()
  assert np.array_equal(np.asarray(out.X.sum(1)).ravel(),admitted.n_counts.to_numpy())
  assert np.isfinite(out.X.data).all() and (out.X.data>=0).all() and (out.X.data==np.floor(out.X.data)).all()
  out.obs=admitted.set_index('cell_id').loc[out.obs_names].copy()
  for col in out.obs:
   if out.obs[col].dtype==object:out.obs[col]=out.obs[col].fillna('').astype(str)
  out.uns['admission_scope']='Expanded fibroblast-expression candidates; identity unresolved; no verified independent patient N'
  out.uns['feature_scope']='Author merged available-cell features; original per-library measurement mask unavailable'
  target=p/'objects'/(sm+'_expanded_candidates.h5ad');out.write_h5ad(target,compression='gzip')
  back=ad.read_h5ad(target,backed='r');assert back.shape==out.shape and back.obs_names.equals(out.obs_names);assert (back.X[:]!=out.X).nnz==0;back.file.close()
  obj=dict(sample=sm,path=str(target),cells=out.n_obs,features=out.n_vars,sha256=hashlib.sha256(target.read_bytes()).hexdigest(),raw_counts_roundtrip=True)
  sy=out.var.symbol.astype(str);feature=out.var_names.astype(str)
  for gene in genes:
   g=gene['gene'];match=np.flatnonzero((feature==gene['feature_ids']) | sy.isin(aliases['aliases'][g]) | feature.isin(aliases['aliases'][g]))
   assert len(match)<=1,(sm,g,match)
   coverage.append(dict(sample=sm,gene=g,feature_present_in_author_object=bool(len(match)),reference_row=feature[match[0]] if len(match) else '',original_measurement_mask_available=False))
   for name,mask in [('expanded',np.ones(len(admitted),dtype=bool)),('ambient_supported',admitted.candidate_ambient_supported.to_numpy()),('author_CAF_agreement',admitted.author_CAF.to_numpy()),('author_tumor_excluded',(~admitted.author_tumor_conflict).to_numpy())]:
    den=float(admitted.loc[mask,'n_counts'].sum());n=int(mask.sum());raw=float(out.X[mask][:,match].sum()) if len(match) and n else np.nan
    det=float((out.X[mask][:,match]>0).mean()) if len(match) and n else np.nan
    expr.append(dict(sample=sm,putative_group=summary['putative_group'],stratum=name,gene=g,feature=gene['feature_ids'],cells=n,library_UMIs=den,raw_count=raw,detected_fraction=det,expression=float(np.log2(raw/den*1e4+.1)) if den>0 and np.isfinite(raw) else np.nan,feature_present=bool(len(match)),independent_patient_verified=False,inferential_P=np.nan))
 return summary,admitted,expr,coverage,obj,dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest())
with concurrent.futures.ThreadPoolExecutor(2) as pool:results=list(pool.map(one,run))
pd.DataFrame([z[0] for z in results]).sort_values('sample').to_csv(p/'tables/final_library_disposition.tsv',sep='\t',encoding='utf-8',index=False)
cells=pd.concat([z[1] for z in results],ignore_index=True);assert cells.cell_id.is_unique
cells.to_csv(p/'tables/NCI_expanded_candidate_registry.tsv.gz',sep='\t',encoding='utf-8',index=False)
expr=pd.DataFrame([r for z in results for r in z[2]]);expr.to_csv(p/'tables/NCI_14gene_library_descriptive.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
pd.DataFrame([r for z in results for r in z[3]]).to_csv(p/'tables/NCI_14gene_feature_coverage.tsv',sep='\t',encoding='utf-8',index=False)
pd.DataFrame([z[5] for z in results]).to_csv(p/'resources/EVIDENCE_INPUT_SHA256.tsv',sep='\t',encoding='utf-8',index=False)
(p/'objects/COUNT_COLLECTION.json').write_text(json.dumps([z[4] for z in results if z[4]],indent=2),encoding='utf-8')
groups=[]
for (group,stratum,gene),z in expr.groupby(['putative_group','stratum','gene']):
 den=float(z.library_UMIs.sum());measured=bool(z.loc[z.cells.gt(0),'feature_present'].all()) and bool(z.cells.sum()>0)
 raw=float(z.raw_count.sum()) if measured else np.nan
 groups.append(dict(putative_group=group,stratum=stratum,gene=gene,libraries=';'.join(z['sample']),cells=int(z.cells.sum()),raw_count=raw,library_UMIs=den,expression=float(np.log2(raw/den*1e4+.1)) if den>0 and measured else np.nan,grouping_verified=False,grouping_basis='suffix hypothesis only; descriptive aggregation, not independent biological N',P=np.nan))
pd.DataFrame(groups).to_csv(p/'tables/NCI_14gene_putative_group_descriptive.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
res=dict(status='processing_outputs_complete_QA_pending',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),libraries=19,input_cells=sum(z[0]['cells'] for z in results),expanded_candidates=len(cells),libraries_with_candidates=int(cells['sample'].nunique()),putative_groups_with_candidates=int(cells.putative_patient_group.nunique()),verified_independent_N=None,ambient_supported=int(cells.candidate_ambient_supported.sum()),author_CAF_agreement=int(cells.author_CAF.sum()),author_tumor_conflict=int(cells.author_tumor_conflict.sum()),reference_pass=sum(z[0]['reference_pass'] for z in results),all_counts_roundtrip_pass=True,new_inferential_tests=False)
(p/'BUILD_COMPLETE.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
print(json.dumps(res,indent=2),flush=True)
