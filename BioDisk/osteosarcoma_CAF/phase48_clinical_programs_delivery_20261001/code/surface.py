import sys
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase46_spatial_functional_20261001/code')
from common import *
import anndata as ad,gc
Q=P/'phase48_clinical_programs_delivery_20261001';E=P/'phase47_external_candidate_integration_20261001';D=P/'phase44_patient_mechanism_20261001';L=P/'phase45_LR_expression_20261001';genes=['FGFR2','MRC2','SDC1','PTK7','CA12','MMP14','APP'];pairs=[('MRC2','MMP14'),('MRC2','FGFR2'),('FGFR2','CA12')];receivers=['T_NK','Myeloid','Endothelial','Mural','Fibroblast_other','CNV_altered_unresolved','Other_nonfrozen_unresolved'];a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=a.obs.copy();a.file.close();allfrozen=set(ob.index);m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);cw=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);symbols=dict(zip(cw.feature,cw.canonical_symbol));quality=dict(zip(cw.feature,~cw.mapping_method.eq('HGNC_Ensembl_approved')));cache={};co=[];cover=[]
roles=pd.read_csv(D/'resources/intercell/receptor_roles.tsv',sep='\t');accepted=roles[roles.genesymbol.isin(genes)&roles.plasma_membrane_transmembrane.astype(str).str.lower().eq('true')];assert set(accepted.genesymbol)==set(genes);accepted.to_csv(Q/'resources/new_surface_membrane_evidence.tsv',sep='\t',index=False)
dump(Q/'resources/SURFACE_PROTOCOL_FROZEN.json',{'genes':genes,'selection':'supported F3 main raw patient P<.05 upregulation and OmniPath receptor membrane annotation; genome-wide discovery not original8marker gate','pairs':pairs,'statistics':'patient matched same site F3/F4 minus competitive annotated receiver; exact patient signflip; BH across7 genes by layer/source/class/receiver','CNV_altered_not_confirmed_tumor':True,'co_detection':'raw positive counts; no protein positivity claim; absent assayed gene NA','external_features_heldout':'all7 already in frozen F3-up discovery signature excluded from phase47 mapper','not_validated_sorting_gate':True})
# Reuse exactly the receiver annotation expressions from the accepted LR code.
old=(L/'code/stageD_LR_expression.py').read_text().splitlines();start=next(i for i,l in enumerate(old) if l.startswith('  provider='));end=next(i for i,l in enumerate(old) if l.startswith('  # Frozen sender'))
annotation='\n'.join(l[2:] for l in old[start:end]);def_flag=lambda x:x.astype(str).str.lower().isin(['true','1']);flag=def_flag
for row in m.itertuples(index=False):
 f=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad';a=ad.read_h5ad(f,backed='r');pref={}
 for j,en in enumerate(a.var_names.astype(str)):
  g=symbols.get(en,'');key=(quality.get(en,1),en)
  if g in genes and (g not in pref or key<pref[g][0]):pref[g]=(key,j)
 present=np.array([g in pref for g in genes]);idx=[pref[g][1] for g in genes if g in pref];x=a.X[:,:].tocsr();raw=np.full((a.n_obs,len(genes)),np.nan);raw[:,present]=x[:,idx].toarray();cells=pd.read_csv(f.parent/'cells.tsv',sep='\t',keep_default_na=False).set_index('cell_id');lib=pd.to_numeric(cells.loc[a.obs_names,'n_counts']).to_numpy();meta=a.obs.copy();nonfrozen=~meta.index.isin(allfrozen);exec(annotation)
 here=ob[ob['sample'].astype(str).eq(row.sample)]
 for layer in ['supported','expanded']:
  vv=here[here.in_supported.astype(bool)] if layer=='supported' else here
  for cl in ['F3','F4']:mask[layer+'__'+cl]=meta.index.isin(vv.index[vv[layer+'_class_id'].astype(str).eq(cl)])
 stored={};baseline=np.load(L/'objects'/f'{row.sample}__LR_group_pseudobulk.npz',allow_pickle=True)
 for group,v in mask.items():
  assert int(v.sum())==int(baseline['ncells'][list(baseline['groups']).index(group)])
  den=float(lib[v].sum());stored[group]={'counts':np.nansum(raw[v],axis=0),'library':den,'present':present&(v.sum()>0),'cells':int(v.sum())}
  if group in ['supported__F3','expanded__F4']:
   for g1,g2 in pairs:
    n=int(v.sum());available=present[genes.index(g1)]&present[genes.index(g2)];both=int(((raw[v,genes.index(g1)]>0)&(raw[v,genes.index(g2)]>0)).sum()) if available else np.nan;co.append(dict(dataset=row.dataset,sample=row.sample,patient_id=row.donor_record,group=group,pair=g1+'_'+g2,cells=n,both_detected=both,fraction_both=both/n if n else np.nan,assayed=bool(available)))
 cache[row.sample]=stored;cover.append(dict(sample=row.sample,dataset=row.dataset,called_cells=a.n_obs,present_genes=';'.join(g for g,ok in zip(genes,present) if ok),missing_genes=';'.join(g for g,ok in zip(genes,present) if not ok)));a.file.close();del a,x,raw;gc.collect();print('SURFACE_SAMPLE',row.sample,flush=True)
def agg(samples,group):
 v=[cache[s][group] for s in samples];cells=sum(t['cells'] for t in v);lib=sum(t['library'] for t in v);counts=sum((t['counts'] for t in v),np.zeros(len(genes)));pp=np.logical_and.reduce([t['present'] for t in v if t['cells']>0]) if cells else np.zeros(len(genes),bool);log=np.log2(counts/lib*1e4+.1) if lib>0 else np.full(len(genes),np.nan);log[~pp]=np.nan;return log
rows=[]
for layer,cl in [('supported','F3'),('expanded','F4')]:
 for (ds,donor),t in m.groupby(['dataset','donor_record']):
  for receiver in receivers:
   site=[]
   for _,z in t.groupby('site',dropna=False):
    sm=list(z['sample']);left=agg(sm,layer+'__'+cl);right=agg(sm,receiver);site.append(left-right)
   mat=np.vstack(site);den=np.isfinite(mat).sum(0);v=np.divide(np.nansum(mat,axis=0),den,out=np.full(len(genes),np.nan),where=den>0)
   for gene,value in zip(genes,v):rows.append(dict(layer=layer,class_id=cl,dataset=ds,patient_id=donor,receiver=receiver,gene=gene,effect=value))
pat=pd.DataFrame(rows);pat.to_csv(Q/'tables/surface_competing_lineage_patient_effects.tsv',sep='\t',index=False);stats=test_table(pat,['layer','class_id','dataset','receiver','gene'],'effect',lambda r:'surface_'+r.layer+'_'+r.class_id+'_'+r.dataset+'_'+r.receiver);stats.to_csv(Q/'tables/surface_competing_lineage_statistics.tsv',sep='\t',index=False);pd.DataFrame(co).to_csv(Q/'tables/new_surface_pair_RNA_codetection.tsv',sep='\t',index=False);pd.DataFrame(cover).to_csv(Q/'tables/new_surface_gene_coverage.tsv',sep='\t',index=False)
external=[];epairs=[]
for path in sorted((P/'phase37_requested_source_search_20260927/objects/GSE237070_reconstructed').glob('*/full_native_called_counts.h5ad')):
 a=ad.read_h5ad(path);sm=str(a.obs['sample'].iloc[0]);patient=str(a.obs.patient.iloc[0]);tissue=str(a.obs.tissue.iloc[0]);ann=pd.read_csv(E/'tables'/f'{sm}_reference_mapping.tsv.gz',sep='\t').set_index('cell_id').loc[a.obs_names];gm=feature_map(a);x=a.X.tocsr();lib=np.asarray(x.sum(1)).ravel();fib=ann.lineage_primary_blind.eq('Fibroblast').to_numpy()
 for layer,cl in [('supported','F3'),('expanded','F4')]:
  lab=ann[layer+'_reference_like'].astype(str).to_numpy();target=fib&(lab==cl);rest=fib&np.isin(lab,['F1','F2','F3','F4','F5'])&(lab!=cl)
  for gene in genes:
   vals=[]
   for mask in [target,rest]:vals.append(np.log2(float(x[mask,gm[gene]].sum())/lib[mask].sum()*1e4+.1) if gene in gm and mask.any() and lib[mask].sum()>0 else np.nan)
   external.append(dict(layer=layer,class_id=cl,patient_id=patient,tissue=tissue,sample=sm,gene=gene,effect=vals[0]-vals[1],target_cells=int(target.sum()),rest_cells=int(rest.sum())))
  for g1,g2 in pairs:
   n=int(target.sum());available=g1 in gm and g2 in gm;both=int(((x[target,gm[g1]].toarray().ravel()>0)&(x[target,gm[g2]].toarray().ravel()>0)).sum()) if available else np.nan;epairs.append(dict(layer=layer,class_id=cl,patient_id=patient,tissue=tissue,sample=sm,pair=g1+'_'+g2,cells=n,both_detected=both,fraction_both=both/n if n else np.nan))
ef=pd.DataFrame(external);ef=ef.groupby(['layer','class_id','patient_id','tissue','gene']).effect.mean().reset_index();ef.to_csv(Q/'tables/external_surface_patient_effects.tsv',sep='\t',index=False);es=test_table(ef,['layer','class_id','tissue','gene'],'effect',lambda r:'external_surface_'+r.layer+'_'+r.class_id+'_'+r.tissue);es.to_csv(Q/'tables/external_surface_gene_statistics.tsv',sep='\t',index=False);pd.DataFrame(epairs).to_csv(Q/'tables/external_new_surface_pair_RNA_codetection.tsv',sep='\t',index=False)
qa={'samples':len(cache),'called_cells':sum(t['called_cells'] for t in cover),'surface_genes':genes,'genes_statistics':len(stats),'nominal':int(stats.P.lt(.05).sum()),'RNA_pairs':pairs,'cross_lineage_not_assumed_specific':True,'live_protein_gate_unvalidated':True,'source_normalized_cell_groups_match_phase45':True};dump(Q/'delivery/SURFACE_QA.json',qa);(Q/'SURFACE_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
