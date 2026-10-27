import sys
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase46_spatial_functional_20261001/code')
from common import *
import anndata as ad,gc
def dump_once(path,value):
 if path.exists():
  old=json.loads(path.read_text());assert {k:v for k,v in old.items() if k!="at"}=={k:v for k,v in value.items() if k!="at"}
 else:dump(path,value)
from scipy import sparse
E=P/'phase47_external_candidate_integration_20261001';B=P/'phase30_two_layer_classification_20260926';old=P/'phase38_GSE237070_expression_validation_20260928';S=P/'phase37_requested_source_search_20260927'
pr=json.loads((O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json').read_text());programs=pr['programs'];sigs=pr['signatures'];genes=pr['external_primary_gene_panel'];surface=['PDGFRB','MCAM','GPC3','THY1','PDGFRA','ITGA11','FAP','LRRC15'];lineage=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text());blind={k:[g for g in v if g not in genes+surface] for k,v in lineage.items()}
excluded=set(sum(programs.values(),[])+sum([v['up']+v['down'] for v in sigs.values()],[])+genes+surface);dump_once(E/'resources/EXTERNAL_PROTOCOL_FROZEN.json',{'at':datetime.datetime.now().astimezone().isoformat(),'functional_signature_SHA256':sha(O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json'),'mapping_excluded_genes':sorted(excluded),'lineage_primary_gene_and_surface_heldout':blind,'broader_signature_lineage_overlap_disclosed':True,'reference_profile_full_library_normalization':True,'statistic':'target_reference-like vs other-reference-like full-library pseudobulk log2(CP10k+0.1), source-patient average within tissue; N<=2 primary exploratory only','no_query_fitting':True,'confirmed_CAF_identity_not_claimed':True,'families':'external tissue/layer/endpoint kind; each F3 or F4 contrast separate; all registered endpoints with NA preserved'})
queries=sorted((S/'objects/GSE237070_reconstructed').glob('*/full_native_called_counts.h5ad'));assert len(queries)==6
first=ad.read_h5ad(queries[0],backed='r');qfeatures=set(first.var_names);first.file.close();a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad');gm=feature_map(a);reverse={j:g for g,j in gm.items()};refs={};mappingrows=[];lib=a.obs.full_library_counts.to_numpy(float);assert np.allclose(np.asarray(a.X.sum(1)).ravel(),lib)
for layer in ['supported','expanded']:
 hv=ad.read_h5ad(B/'objects'/f'{layer}_integrated.h5ad',backed='r');features=[g for g in hv.var_names[hv.var.highly_variable] if g in qfeatures and g in a.var_names and reverse.get(a.var_names.get_loc(g)) not in excluded];hv.file.close();assert len(features)>1
 idx=a.var_names.get_indexer(features);v=a.X[:,idx].astype(float).multiply((1e4/lib)[:,None]).tocsr();v.data=np.log1p(v.data);cl=a.obs[layer+'_class_id'].astype(str).to_numpy();donors=a.obs.donor_record.astype(str).to_numpy();classes=sorted(set(cl)-{''});pro=[]
 for z in classes:
  ds=sorted(set(donors[cl==z]));means=[np.asarray(v[(cl==z)&(donors==d)].mean(0)).ravel() for d in ds];pro.append(np.mean(means,axis=0));mappingrows.append(dict(layer=layer,class_id=z,reference_cells=int((cl==z).sum()),reference_donors=len(ds),mapping_features=len(features)))
 pro=np.asarray(pro);pro-=pro.mean(1,keepdims=True);pro/=np.linalg.norm(pro,axis=1,keepdims=True);refs[layer]=(features,classes,pro);pd.DataFrame({'feature':features}).to_csv(E/'resources'/f'{layer}_mapping_features.tsv',sep='\t',index=False)
 del v
pd.DataFrame(mappingrows).to_csv(E/'tables/reference_mapping_coverage.tsv',sep='\t',index=False);del a;gc.collect();aggregates=[];coverage=[];account=[];coexpr=[]
for path in queries:
 a=ad.read_h5ad(path);sm=str(a.obs['sample'].iloc[0]);patient=str(a.obs.patient.iloc[0]);tissue=str(a.obs.tissue.iloc[0]);x=a.X.tocsr();lib=np.asarray(x.sum(1)).ravel();assert np.allclose(lib,a.obs.n_counts);gm=feature_map(a);obs=a.obs.copy();obs.index.name='cell_id';score=[]
 for label,gs in blind.items():
  have=[g for g in gs if g in gm];v=np.log1p(x[:,[gm[g] for g in have]].toarray()/np.where(lib>0,lib,np.nan)[:,None]*1e4).mean(1) if have else np.full(a.n_obs,np.nan);score.append(v);coverage.append(dict(sample=sm,kind='lineage_blind',endpoint=label,requested=len(gs),present=len(have),missing=';'.join(sorted(set(gs)-set(have)))))
 v=np.column_stack(score);safe=np.where(np.isfinite(v),v,-np.inf);mx=safe.max(1);labels=np.array(list(blind),dtype=object)[safe.argmax(1)];labels[(mx<=0)|(np.sum(safe==mx[:,None],axis=1)!=1)]='Unresolved';obs['lineage_primary_blind']=labels;fib=labels=='Fibroblast'
 for layer,(features,classes,pro) in refs.items():
  ids=a.var_names.get_indexer(features);assert (ids>=0).all();v=x[fib][:,ids].astype(float).multiply((1e4/lib[fib])[:,None]).tocsr();v.data=np.log1p(v.data);norm=np.sqrt(np.maximum(np.asarray(v.multiply(v).sum(1)).ravel()-np.asarray(v.sum(1)).ravel()**2/v.shape[1],0));corrs=np.asarray(v@pro.T)/np.where(norm>0,norm,np.nan)[:,None];safec=np.where(np.isfinite(corrs),corrs,-np.inf);mx=safec.max(1);lab=np.array(classes,dtype=object)[safec.argmax(1)];lab[(norm==0)|(np.sum(safec==mx[:,None],axis=1)!=1)]='Unresolved';obs[layer+'_reference_like']='not_fibroblast_like';obs.loc[fib,layer+'_reference_like']=lab
  oldann=pd.read_csv(old/'tables'/f'{sm}_cell_annotations.tsv.gz',sep='\t').set_index('cell_id');assert oldann.index.equals(obs.index);obs[layer+'_previous_reference_like']=oldann[layer+'_reference_like'];account.append(dict(sample=sm,patient_id=patient,tissue=tissue,layer=layer,fibroblast_like=int(fib.sum()),reference_label_agreement=float(np.mean(obs.loc[fib,layer+'_reference_like']==oldann.loc[fib,layer+'_reference_like'])) if fib.any() else np.nan,mapping_features=len(features),cells=a.n_obs))
  for cl in ['F3','F4']:
   if cl not in classes:continue
   for side,mask in [('target',fib&(obs[layer+'_reference_like'].to_numpy()==cl)),('rest',fib&np.isin(obs[layer+'_reference_like'],[z for z in classes if z!=cl]))]:
    if not mask.any():continue
    den=float(lib[mask].sum());raw=np.asarray(x[mask].sum(0)).ravel();lg=np.log2(raw/den*1e4+.1) if den>0 else np.full(a.n_vars,np.nan)
    endpoints=[('gene',g,lg[gm[g]] if g in gm else np.nan) for g in genes]
    for name,gs in programs.items():
     have=[g for g in gs if g in gm];val=np.mean([lg[gm[g]] for g in have]) if have else np.nan;endpoints.append(('program',name,val));coverage.append(dict(sample=sm,kind='program',endpoint=name,requested=len(gs),present=len(have),missing=';'.join(sorted(set(gs)-set(have)))))
    for name,sg in sigs.items():
     up=[lg[gm[g]] for g in sg['up'] if g in gm];down=[lg[gm[g]] for g in sg['down'] if g in gm];val=np.mean(up)-np.mean(down) if up and down else np.nan;endpoints.append(('signature',name,val));coverage.append(dict(sample=sm,kind='signature',endpoint=name,requested=len(sg['up'])+len(sg['down']),present=len(up)+len(down),missing=';'.join(sorted(set(sg['up']+sg['down'])-set(gm)))))
    for kind,name,val in endpoints:aggregates.append(dict(sample=sm,patient_id=patient,tissue=tissue,layer=layer,class_id=cl,side=side,kind=kind,endpoint=name,value=val,cells=int(mask.sum()),full_library_UMI=den))
   target=fib&(obs[layer+'_reference_like'].to_numpy()==cl)
   for g1,g2 in [('GPC3','THY1'),('PDGFRB','MCAM')]:
    n=int(target.sum());val=int(((x[target,gm[g1]].toarray().ravel()>0)&(x[target,gm[g2]].toarray().ravel()>0)).sum()) if g1 in gm and g2 in gm else np.nan;coexpr.append(dict(sample=sm,patient_id=patient,tissue=tissue,layer=layer,class_id=cl,pair=g1+'_'+g2,cells=n,both_detected=val,fraction_both=val/n if n else np.nan,RNA_not_protein=True))
 obs.to_csv(E/'tables'/f'{sm}_reference_mapping.tsv.gz',sep='\t');print('EXTERNAL_PROCESSED',sm,a.n_obs,int(fib.sum()),flush=True);del a,x,obs;gc.collect()
agg=pd.DataFrame(aggregates);agg.to_csv(E/'tables/external_sample_pseudobulk_values.tsv',sep='\t',index=False);pd.DataFrame(account).to_csv(E/'tables/reference_migration_summary.tsv',sep='\t',index=False);pd.DataFrame(coverage).drop_duplicates().to_csv(E/'tables/external_gene_coverage.tsv',sep='\t',index=False);pd.DataFrame(coexpr).to_csv(E/'tables/external_surface_RNA_codetection.tsv',sep='\t',index=False)
keys=['sample','patient_id','tissue','layer','class_id','kind','endpoint'];w=agg.pivot(index=keys,columns='side',values='value').reset_index();w['effect']=w['target']-w['rest'];pat=w.groupby(['patient_id','tissue','layer','class_id','kind','endpoint']).effect.mean().reset_index();pat.to_csv(E/'tables/external_patient_effects.tsv',sep='\t',index=False);stats=test_table(pat,['tissue','layer','class_id','kind','endpoint'],'effect',lambda r:'external_'+r.tissue+'_'+r.layer+'_'+r.class_id+'_'+r.kind);stats.to_csv(E/'tables/external_patient_statistics.tsv',sep='\t',index=False)
lineage_overlap=[]
for name,sg in sigs.items():
 for k,v in blind.items():lineage_overlap.append(dict(signature=name,lineage=k,remaining_panel_genes=len(v),overlap_with_signature=';'.join(sorted(set(v)&set(sg['up']+sg['down']))),all_signature_heldout_panel_remaining=len(set(v)-set(sg['up']+sg['down']))))
pd.DataFrame(lineage_overlap).to_csv(E/'tables/signature_lineage_circularity.tsv',sep='\t',index=False)
qa={'samples':len(queries),'patients':agg.patient_id.nunique(),'cells':sum(t['cells'] for t in account if t['layer']=='supported'),'primary_patients':agg.loc[agg.tissue.eq('primary'),'patient_id'].nunique(),'statistics':len(stats),'nominal':int(stats.P.lt(.05).sum()),'FDR':int(stats.FDR.lt(.05).sum()),'query_signature_genes_used_in_mapping':0,'identity_confirmed':False,'reference_migration_not_subtype_validation':True,'cross_cohort_patient_independence':'unverified; no pooled biological N','CNV_doublet_unassessed':True};assert qa['patients']==4 and qa['cells']==47059 and qa['primary_patients']==2
dump(E/'delivery/EXTERNAL_QA.json',qa);(E/'EXTERNAL_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
