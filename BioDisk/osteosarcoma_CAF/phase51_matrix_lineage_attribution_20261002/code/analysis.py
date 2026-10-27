from pathlib import Path
import os,sys,json,hashlib,datetime,gc
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';L=P/'phase45_LR_expression_20261001';O=P/'phase51_matrix_lineage_attribution_20261002'
sys.path.insert(0,str(P/'phase46_spatial_functional_20261001/code'))
from common import test_table,bh
import numpy as np,pandas as pd,anndata as ad
pr=json.loads((O/'resources/PROTOCOL_FROZEN.json').read_text());genes=pr['genes'];receivers=pr['comparators']
partition=['expanded__F'+str(i) for i in range(1,6)]+receivers[:7]
groupnames=partition+['supported__F3','CAF_frozen_all']+receivers[7:]
cw=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False)
symbols=dict(zip(cw.feature,cw.canonical_symbol));quality=dict(zip(cw.feature,~cw.mapping_method.eq('HGNC_Ensembl_approved')))
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=a.obs.copy();a.file.close();allfrozen=set(ob.index)
m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
annotation=(O/'resources/receiver_annotation_frozen.py').read_text()
flag=lambda x:x.astype(str).str.lower().isin(['true','1'])
cache={};sample_rows=[];input_rows=[];mapping=[];identity_rows=[]
for row in m.itertuples(index=False):
 path=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad'
 a=ad.read_h5ad(path,backed='r');pref={}
 for j,en in enumerate(a.var_names.astype(str)):
  g=symbols.get(en,'');key=(quality.get(en,1),en)
  if g in genes and (g not in pref or key<pref[g][0]):pref[g]=(key,j)
 present=np.array([g in pref for g in genes]);idx=[pref[g][1] for g in genes if g in pref]
 x=a.X[:,:].tocsr();raw=np.full((a.n_obs,3),np.nan);raw[:,present]=x[:,idx].toarray()
 cells=pd.read_csv(path.parent/'cells.tsv',sep='\t',keep_default_na=False).set_index('cell_id')
 lib=pd.to_numeric(cells.loc[a.obs_names,'n_counts']).to_numpy();assert np.array_equal(lib,np.asarray(x.sum(1)).ravel())
 meta=a.obs.copy();nonfrozen=~meta.index.isin(allfrozen);exec(annotation)
 here=ob[ob['sample'].astype(str).eq(row.sample)]
 assert here.index.isin(meta.index).all()
 for cl in ['F1','F2','F3','F4','F5']:mask['expanded__'+cl]=meta.index.isin(here.index[here.expanded_class_id.astype(str).eq(cl)])
 vv=here[here.in_supported.astype(bool)]
 mask['supported__F3']=meta.index.isin(vv.index[vv.supported_class_id.astype(str).eq('F3')])
 mask['CAF_frozen_all']=meta.index.isin(here.index)
 osteo=provider.str.contains(r'osteoblast|osteogenic',regex=True).to_numpy()|hist.isin(['osteolike','osteogenic','osteoblast']).to_numpy()
 mask['Annotated_osteogenic_nonfrozen']=osteo&nonfrozen
 mask['CNV_altered_osteogenic_annotated']=osteo&alteration&~reference&nonfrozen
 assert np.all(np.column_stack([mask[g] for g in partition]).sum(1)==1)
 baseline=np.load(L/'objects'/f'{row.sample}__LR_group_pseudobulk.npz',allow_pickle=True)
 for g in receivers[:7]+['supported__F3','expanded__F4']:
  assert int(mask[g].sum())==int(baseline['ncells'][list(baseline['groups']).index(g)])
 stored={}
 for g in groupnames:
  v=mask[g];n=int(v.sum());counts=np.nansum(raw[v],axis=0);den=float(lib[v].sum())
  detected=np.array([int((raw[v,j]>0).sum()) if present[j] else np.nan for j in range(3)])
  both=int(((raw[v,1]>0)&(raw[v,2]>0)).sum()) if present[1:].all() else np.nan
  all3=int((raw[v]>0).all(axis=1).sum()) if present.all() else np.nan
  stored[g]={'counts':counts,'library':den,'cells':n,'present':present.copy(),'detected':detected,'both_LM':both,'all3':all3}
  for j,gene in enumerate(genes):
   sample_rows.append(dict(dataset=row.dataset,sample=row.sample,patient_id=row.donor_record,site=row.site,group=g,gene=gene,cells=n,library=den,assayed=bool(present[j]),counts=counts[j] if present[j] else np.nan,detected_cells=detected[j],fraction_detected=detected[j]/n if n else np.nan,LOX_MMP14_both_cells=both,all3_cells=all3,group_scheme='exhaustive_partition' if g in partition else 'overlapping_view'))
  identity_rows.append(dict(dataset=row.dataset,sample=row.sample,group=g,cells=n,CNV_concordant_alteration=int((v&alteration).sum()),reference_cells=int((v&reference).sum()),annotated_osteogenic=int((v&osteo).sum()),annotation_conflicts=int((v&ambiguous).sum())))
 counts_all=np.nansum(raw,axis=0)
 assert np.array_equal(sum((stored[g]['counts'] for g in partition),np.zeros(3)),counts_all)
 assert sum(stored[g]['cells'] for g in partition)==a.n_obs
 assert np.isclose(sum(stored[g]['library'] for g in partition),lib.sum(),rtol=0,atol=0)
 cache[row.sample]=stored
 for g in genes:mapping.append(dict(sample=row.sample,dataset=row.dataset,gene=g,present=g in pref,feature_id=pref[g][0][1] if g in pref else '',feature_index=pref[g][1] if g in pref else -1))
 input_rows.append(dict(sample=row.sample,dataset=row.dataset,called_cells=a.n_obs,all_library_counts=float(lib.sum()),raw_SHA256=hashlib.sha256(path.read_bytes()).hexdigest(),cells_SHA256=hashlib.sha256((path.parent/'cells.tsv').read_bytes()).hexdigest(),identity_SHA256=hashlib.sha256((P/'phase26_client_reanalysis_20260922/tables'/f'{row.sample}_identity_evidence.tsv.gz').read_bytes()).hexdigest(),partition_verified=True))
 np.savez_compressed(O/'objects'/f'{row.sample}_three_gene_groups.npz',genes=np.array(genes),groups=np.array(groupnames),counts=np.stack([stored[g]['counts'] for g in groupnames]),libraries=np.array([stored[g]['library'] for g in groupnames]),ncells=np.array([stored[g]['cells'] for g in groupnames]),presence=np.stack([stored[g]['present'] for g in groupnames]))
 print('SAMPLE_COMPLETE',row.sample,a.n_obs,flush=True);a.file.close();del a,x,raw;gc.collect()
pd.DataFrame(sample_rows).to_csv(O/'tables/sample_group_expression.tsv.gz',sep='\t',index=False)
pd.DataFrame(identity_rows).to_csv(O/'tables/group_identity_evidence.tsv',sep='\t',index=False)
pd.DataFrame(input_rows).to_csv(O/'tables/input_QA.tsv',sep='\t',index=False)
pd.DataFrame(mapping).to_csv(O/'tables/gene_feature_mapping.tsv',sep='\t',index=False)
def aggregate(samples,group):
 v=[cache[s][group] for s in samples];cells=sum(t['cells'] for t in v);lib=sum(t['library'] for t in v)
 counts=sum((t['counts'] for t in v),np.zeros(3))
 # Presence requires assayed features in all contributing nonempty groups.
 ok=np.logical_and.reduce([t['present'] for t in v if t['cells']>0]) if cells else np.zeros(3,bool)
 log=np.log2(counts/lib*1e4+.1) if lib>0 else np.full(3,np.nan);log[~ok]=np.nan
 return log,counts,cells,lib,ok
rows=[];site_rows=[]
for (ds,donor),t in m.groupby(['dataset','donor_record']):
 for site,ss in t.groupby('site',dropna=False):
  sm=list(ss['sample']);all_counts=sum((sum((cache[s][g]['counts'] for g in partition),np.zeros(3)) for s in sm),np.zeros(3))
  all_present=np.logical_and.reduce([cache[s]['CAF_frozen_all']['present'] for s in sm])
  for group in groupnames:
   log,ct,n,lib,ok=aggregate(sm,group)
   vs=[cache[s][group] for s in sm]
   detected=np.sum(np.stack([v['detected'] for v in vs]),axis=0);both=sum(v['both_LM'] for v in vs);all3=sum(v['all3'] for v in vs)
   for j,gene in enumerate(genes):
    site_rows.append(dict(dataset=ds,patient_id=donor,site=site,group=group,gene=gene,cells=n,library=lib,counts=ct[j] if all_present[j] else np.nan,log_expression=log[j],fraction_detected=detected[j]/n if n else np.nan,gene_capture_share=ct[j]/all_counts[j] if all_present[j] and all_counts[j]>0 else np.nan,LOX_MMP14_fraction=both/n if n else np.nan,all3_fraction=all3/n if n else np.nan,group_scheme='exhaustive_partition' if group in partition else 'overlapping_view'))
 for target in pr['targets']:
  for receiver in receivers:
   deltas=[]
   for _,ss in t.groupby('site',dropna=False):
    sm=list(ss['sample']);deltas.append(aggregate(sm,target)[0]-aggregate(sm,receiver)[0])
   mat=np.vstack(deltas);den=np.isfinite(mat).sum(0);v=np.divide(np.nansum(mat,axis=0),den,out=np.full(3,np.nan),where=den>0)
   for gene,d in zip(genes,v):rows.append(dict(dataset=ds,patient_id=donor,target=target,receiver=receiver,gene=gene,effect=d))
pat=pd.DataFrame(rows);pat.to_csv(O/'tables/patient_matched_effects.tsv',sep='\t',index=False)
st=test_table(pat,['dataset','target','receiver','gene'],'effect',lambda r:'lineage_'+r.dataset+'_'+r.target)
assert st.groupby('family').size().eq(27).all()
st.to_csv(O/'tables/lineage_statistics.tsv',sep='\t',index=False)
site=pd.DataFrame(site_rows);site.to_csv(O/'tables/patient_site_expression_and_shares.tsv',sep='\t',index=False)
patdesc=site.groupby(['dataset','patient_id','group','gene','group_scheme'],dropna=False).agg(log_expression=('log_expression','mean'),fraction_detected=('fraction_detected','mean'),gene_capture_share=('gene_capture_share','mean'),LOX_MMP14_fraction=('LOX_MMP14_fraction','mean'),all3_fraction=('all3_fraction','mean'),captured_cells=('cells','sum')).reset_index()
part=patdesc[patdesc.group_scheme.eq('exhaustive_partition')]
sharesum=part.groupby(['dataset','patient_id','gene']).gene_capture_share.sum(min_count=1)
assert np.allclose(sharesum.dropna(),1,atol=1e-12)
patdesc.to_csv(O/'tables/patient_expression_and_shares.tsv',sep='\t',index=False)
summary=patdesc.groupby(['dataset','group','gene','group_scheme']).agg(N_patients=('patient_id','nunique'),N_expression=('log_expression','count'),mean_log_expression=('log_expression','mean'),median_log_expression=('log_expression','median'),N_detection=('fraction_detected','count'),mean_detection_fraction=('fraction_detected','mean'),N_gene_share=('gene_capture_share','count'),mean_gene_capture_share=('gene_capture_share','mean'),mean_LOX_MMP14_fraction=('LOX_MMP14_fraction','mean'),mean_all3_fraction=('all3_fraction','mean'),captured_cells=('captured_cells','sum')).reset_index()
summary.to_csv(O/'tables/lineage_source_summary.tsv',sep='\t',index=False)
old=pd.read_csv(P/'phase49_audit_revision_20261002/surface/tables/surface_competing_lineage_statistics.tsv',sep='\t',float_precision='round_trip')
old=old[old.gene.eq('MMP14')].copy();old['target']=old.layer+'__'+old.class_id
joined=st[st.gene.eq('MMP14')].merge(old,on=['dataset','target','receiver','gene'],suffixes=('_new','_old'))
assert len(joined)==126
assert np.allclose(joined.mean_effect_new,joined.mean_effect_old,atol=1e-12,rtol=0,equal_nan=True)
assert np.allclose(joined.P_new,joined.P_old,atol=1e-12,rtol=0,equal_nan=True)
joined.to_csv(O/'tables/MMP14_previous_result_QA.tsv',sep='\t',index=False)
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(cache),'called_cells':sum(r['called_cells'] for r in input_rows),'datasets':m.dataset.nunique(),'source_scoped_patient_records':m.donor_record.nunique(),'expanded_frozen_cells':len(ob),'supported_frozen_cells':int(ob.in_supported.astype(bool).sum()),'statistics_rows':len(st),'families':st.family.nunique(),'hypotheses_per_family':27,'partition_count_library_gene_counts_verified':True,'patient_gene_share_sums_verified':True,'phase45_group_counts_reproduced':True,'MMP14_previous_mean_and_P_reproduced':len(joined),'FDR_families_changed_from_surface_analysis':True,'no_cross_source_pooling':True}
assert qa['samples']==107 and qa['called_cells']==876648 and qa['expanded_frozen_cells']==8174
(O/'delivery/ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2));(O/'ANALYSIS_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
