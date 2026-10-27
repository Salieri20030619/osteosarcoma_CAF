# SYNTAX_CHECK: python -m py_compile stageD_LR_expression.py
# RUN_COMMAND: batch_bench/bin/python -u stageD_LR_expression.py
# EXPECTED_SUCCESS_MARKER: STAGE_D_LR_EXPRESSION_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,hashlib,ast,itertools,collections,warnings
import numpy as np,pandas as pd,anndata as ad
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';D=P/'phase44_patient_mechanism_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase45_LR_expression_20261001'
for sub in ['code','logs','tables','objects','resources','delivery','figures']:(O/sub).mkdir(parents=True,exist_ok=True)
def writej(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
def flag(x):return x.astype(str).str.lower().isin(['true','1'])
source=C/'code/stageC_full_gene_statistics.py';nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in ['mean_valid','tests','fdr_by_family']];signs_cache={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'))
rc=1
try:
 assert (D/'CORE_COMPLETE').read_text().strip()=='SUCCESS' and (D/'qa_resume.exit_code').read_text().strip()=='0'
 assert (D/'INTERCELL_DOWNLOAD_COMPLETE').read_text().strip()=='SUCCESS' and (D/'intercell_download.exit_code').read_text().strip()=='0'
 lig=pd.read_csv(D/'resources/intercell/ligand_roles.tsv',sep='\t',keep_default_na=False);rec=pd.read_csv(D/'resources/intercell/receptor_roles.tsv',sep='\t',keep_default_na=False)
 lig=lig[lig.category.eq('ligand')&flag(lig.transmitter)&(flag(lig.secreted)|flag(lig.plasma_membrane_transmembrane)|flag(lig.plasma_membrane_peripheral))]
 rec=rec[rec.category.eq('receptor')&flag(rec.receiver)&flag(rec.plasma_membrane_transmembrane)]
 edges=pd.read_csv(A/'resources/ligand_receptor.tsv',sep='\t',keep_default_na=False);edges=edges[edges.source.isin(set(lig.uniprot))&edges.target.isin(set(rec.uniprot))&flag(edges.is_directed)].copy()
 cross=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False).set_index('feature');symbols=cross.canonical_symbol.to_dict();quality=cross.mapping_method.eq('HGNC_Ensembl_approved').map({True:0,False:1}).to_dict()
 hgnc=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);hgnc=hgnc[hgnc.status.eq('Approved')];approved=set(hgnc.symbol);alias=collections.defaultdict(set)
 for row in hgnc.itertuples(index=False):
  for col in ['prev_symbol','alias_symbol']:
   for s in getattr(row,col,'').split('|'):
    if s:alias[s].add(row.symbol)
 def canon(s):return s if s in approved else (next(iter(alias[s])) if s in alias and len(alias[s])==1 else s)
 edges['ligand']=edges.source_genesymbol.map(canon);edges['receptor']=edges.target_genesymbol.map(canon)
 simple=edges.ligand.isin(approved)&edges.receptor.isin(approved)&~edges.source.str.startswith('COMPLEX:')&~edges.target.str.startswith('COMPLEX:');unassessed=edges[~simple].copy();unassessed['reason']='complex_or_unresolved_subunit_RNA_readout';save(unassessed,'role_annotated_edges_without_single_gene_readout.tsv')
 edges=edges[simple].drop_duplicates(['ligand','receptor']).copy();edges['pair_id']=edges.ligand+'|'+edges.receptor;save(edges,'role_validated_single_gene_LR_network.tsv')
 genes=sorted(set(edges.ligand)|set(edges.receptor));gidx={x:i for i,x in enumerate(genes)};li=np.asarray([gidx[x] for x in edges.ligand]);ri=np.asarray([gidx[x] for x in edges.receptor]);assert len(edges)>0
 roles=['T_NK','Myeloid','Endothelial','Mural','Fibroblast_other','CNV_altered_unresolved','Other_nonfrozen_unresolved'];senders=[f'{layer}__{cls}' for layer in ['supported','expanded'] for cls in ['F2','F3','F4']];groupnames=[v for s in senders for v in [s,s+'__rest']]+roles;gn={x:i for i,x in enumerate(groupnames)}
 frozen=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=frozen.obs.copy();frozen.file.close();m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);allfrozen=set(ob.index)
 writej(O/'resources/LR_PROTOCOL_FROZEN.json',{'eligible_edges':'all role-annotated directed single-gene ligands with secreted/exposed annotation and transmembrane receptors; no consensus-score threshold; complex/unresolved RNA readout retained separately','receiver_assignment':'existing provider/historical annotations; disagreements unresolved; no marker-based new forced calls; all frozen candidates excluded from receiver groups','CNV_receiver':'concordant alteration, not normal reference and not known normal lineage; labelled CNV_altered_unresolved, not confirmed tumor','sender_contrasts':['F2_vs_other_frozen_candidates','F3_vs_other_frozen_candidates','F4_vs_other_frozen_candidates'],'gene_expression':'full-source library pseudobulk CP10k, using source-specific HGNC-approved deterministic feature','LR_score':'log1p(ligand_sender_CP10k)*log1p(receptor_receiver_CP10k)','effect':'target-sender score minus other-candidate sender score, receiver fixed in same patient/site; sites equal weight within patient','missing':'absent ligand/receptor feature or absent sender/receiver cannot be zero-imputed; all-zero score contrasts P/FDR NA','inference':'independent source-scoped patients; two-sided sign flip exact<=15 otherwiseMC19999+observed','family':'layer|source|all_receiver_all_three_sender_LR_hypotheses; separate from TF/genes/clinical','RNA_compatibility':'expression compatibility only, not observed communication, immunosuppression or treatment-causal effect','C2_required':False,'gene_catalog':len(genes),'single_gene_pairs':len(edges),'complex_unassessed':len(unassessed),'no_new_cell_filters':True})
 groups={};availability=[];receivercounts=[]
 for si,row in enumerate(m.itertuples(index=False),1):
  f=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad';a=ad.read_h5ad(f,backed='r');ids=a.var_names.astype(str).to_numpy();pref={}
  for j,en in enumerate(ids):
   sy=symbols.get(en,'');key=(quality.get(en,1),en)
   if sy in gidx and (sy not in pref or key<pref[sy][0]):pref[sy]=(key,j)
  present=np.asarray([g in pref for g in genes]);ix=np.asarray([pref[g][1] for g in genes if g in pref],dtype=int)
  # Load once per sample, never repeatedly decompress the same candidate matrix.
  x=a.X[:,:].tocsr();selected=x[:,ix];cells=pd.read_csv(f.parent/'cells.tsv',sep='\t',keep_default_na=False).set_index('cell_id');assert set(cells.index)==set(a.obs_names);lib=pd.to_numeric(cells.loc[a.obs_names,'n_counts']).to_numpy();meta=a.obs.copy();nonfrozen=~meta.index.isin(allfrozen)
  provider=meta.get('provider_consensus_celltype_annotation',pd.Series('',index=meta.index)).astype(str).str.lower();hist=meta.get('historical_broad_lineage',pd.Series('',index=meta.index)).astype(str).str.lower()
  known={
  'T_NK':provider.str.contains(r'\bt cell\b|natural killer|innate lymphoid',regex=True).to_numpy()|hist.isin(['t_nk','t','nk','tnk']).to_numpy(),
  'Myeloid':provider.str.contains(r'macrophage|monocyte|myeloid|neutrophil|dendritic|phagocyte',regex=True).to_numpy()|hist.isin(['myeloid','dc']).to_numpy(),
  'Endothelial':provider.str.contains('endothelial').to_numpy()|hist.isin(['endo','endothelial']).to_numpy(),
  'Mural':provider.str.contains(r'pericyte|smooth muscle|perivascular',regex=True).to_numpy()|hist.isin(['mural','pericyte','smooth_muscle']).to_numpy(),
  'Fibroblast_other':provider.str.contains(r'fibroblast|stromal cell',regex=True).to_numpy()|hist.isin(['fiblike','fibroblast','fibroblast_like']).to_numpy()}
  bits=np.column_stack(list(known.values()));ambiguous=bits.sum(1)>1;known_union=bits.any(1);mask={k:v&~ambiguous&nonfrozen for k,v in known.items()}
  evidence=pd.read_csv(P/'phase26_client_reanalysis_20260922/tables'/f'{row.sample}_identity_evidence.tsv.gz',sep='\t',keep_default_na=False).set_index('cell_id');assert evidence.index.is_unique and set(evidence.index)==set(meta.index)
  alteration=evidence.CNV_evidence.reindex(meta.index).astype(str).eq('concordant_alteration').to_numpy();reference=flag(cells.loc[meta.index,'reference']).to_numpy();mask['CNV_altered_unresolved']=alteration&~reference&~known_union&nonfrozen;used=np.column_stack(list(mask.values())).any(1);mask['Other_nonfrozen_unresolved']=nonfrozen&~used
  # Frozen sender classifications are supplied by the accepted layer registry, not re-inferred.
  here=ob[ob['sample'].astype(str).eq(row.sample)]
  for layer in ['supported','expanded']:
   label='supported_class_id' if layer=='supported' else 'expanded_class_id';vv=here[here.in_supported.astype(bool)] if layer=='supported' else here
   for cls in ['F2','F3','F4']:
    key=f'{layer}__{cls}';t=set(vv.index[vv[label].astype(str).eq(cls)]);r=set(vv.index[~vv[label].astype(str).eq(cls)]);mask[key]=meta.index.isin(t);mask[key+'__rest']=meta.index.isin(r)
  counts=np.zeros((len(groupnames),len(genes)),float);libraries=np.zeros(len(groupnames));ncell=np.zeros(len(groupnames),int);pp=np.zeros((len(groupnames),len(genes)),bool)
  for name,v in mask.items():
   k=gn[name];ncell[k]=int(v.sum());libraries[k]=lib[v].sum()
   if ncell[k]:counts[k,present]=np.asarray(selected[v].sum(0)).ravel();pp[k]=present
   if name in roles:receivercounts.append({'sample':row.sample,'dataset':row.dataset,'donor_record':row.donor_record,'receiver':name,'cells':int(v.sum()),'annotation_conflicts_in_full_sample':int(ambiguous.sum())})
  cache=O/'objects'/f'{row.sample}__LR_group_pseudobulk.npz';np.savez_compressed(cache,counts=counts,libraries=libraries,ncells=ncell,presence=pp,genes=np.asarray(genes,dtype=str),groups=np.asarray(groupnames,dtype=str));groups[row.sample]={'counts':counts,'libraries':libraries,'ncells':ncell,'presence':pp}
  availability.append({'sample':row.sample,'dataset':row.dataset,'source_features':a.n_vars,'called_cells':a.n_obs,'LR_gene_catalog':len(genes),'LR_features_measured':int(present.sum()),'ligand_receptor_pairs_measured':int((present[li]&present[ri]).sum())});a.file.close();del x,selected;status('RECEIVER_SENDER_COUNTS',sample=row.sample,samples_done=si,total_samples=len(m))
 save(pd.DataFrame(availability),'LR_source_gene_coverage.tsv');save(pd.DataFrame(receivercounts),'receiver_annotation_cell_counts.tsv')
 def aggregate(samples,group):
  k=gn[group];valid=[groups[s] for s in samples if groups[s]['ncells'][k]>0]
  if not valid:return np.zeros(len(genes)),0.,np.zeros(len(genes),bool),0
  return sum(v['counts'][k] for v in valid),sum(v['libraries'][k] for v in valid),np.logical_and.reduce([v['presence'][k] for v in valid]),sum(v['ncells'][k] for v in valid)
 frames=[];qa=[];qa_sites=[];registry=[];contexts=[]
 for layer in ['supported','expanded']:
  for ds,mm in m.groupby('dataset',sort=False):
   sourceframes=[]
   for cls in ['F2','F3','F4']:
    sender=f'{layer}__{cls}'
    for receiver in roles:
     effects=[];names=[];details=[]
     for donor,dm in mm.groupby('donor_record'):
      sites=[];site_detail=[]
      for site,z in dm.groupby('site'):
       samples=list(z['sample']);tc,tl,tp,tn=aggregate(samples,sender);bc,bl,bp,bn=aggregate(samples,sender+'__rest');rr,rl,rp,rn=aggregate(samples,receiver)
       if tl<=0 or bl<=0 or rl<=0:continue
       valid=tp[li]&bp[li]&rp[ri];target=tc[li]/tl*1e4;baseline=bc[li]/bl*1e4;rv=rr[ri]/rl*1e4;delta=(np.log1p(target)-np.log1p(baseline))*np.log1p(rv);delta[~valid]=np.nan;sites.append(delta);site_detail.append((str(site),target,baseline,rv,valid,tn,bn,rn))
      if sites:
       effects.append(mean_valid(np.vstack(sites)));names.append(str(donor));details.append(site_detail)
     n=len(effects);registry.append({'layer':layer,'dataset':ds,'sender_class':cls,'receiver':receiver,'N_matched_patients':n,'status':'patient_test_available' if n>=2 else ('direction_only' if n else 'not_evaluable_no_sender_receiver_overlap')})
     if not n:continue
     matrix=np.vstack(effects);p,N,method=tests(matrix);frame=pd.DataFrame({'layer':layer,'dataset':ds,'sender_class':cls,'receiver':receiver,'pair_id':edges.pair_id.to_numpy(),'ligand':edges.ligand.to_numpy(),'receptor':edges.receptor.to_numpy(),'N':N,'mean_delta':mean_valid(matrix),'positive_patients':np.sum(matrix>0,0),'negative_patients':np.sum(matrix<0,0),'P':p,'test':method,'family':f'{layer}|{ds}|all_sender_receiver_LR'});sourceframes.append(frame)
     ctx=f'{layer}__{ds}__{cls}__{receiver}';np.savez_compressed(O/'objects'/f'{ctx}__patient_LR_effects.npz',effects=matrix,donors=np.asarray(names,dtype=str),pairs=edges.pair_id.to_numpy(dtype=str));contexts.append({'context_id':ctx,'N_matched_patients':n,'effects_file':f'objects/{ctx}__patient_LR_effects.npz'})
     wanted=np.flatnonzero(edges.pair_id.isin(['TGFB1|TGFBR2','CXCL12|CXCR4','HGF|MET','COL1A1|ITGA1']).to_numpy());extra=np.random.default_rng(20261001).choice(len(edges),min(5,len(edges)),replace=False);chosen=np.unique(np.r_[wanted,extra])
     for j in chosen:
      for donor,value,siteinfo in zip(names,matrix[:,j],details):
       qa.append({'context':ctx,'pair_id':edges.iloc[j].pair_id,'donor':donor,'effect':value,'expected_N':N[j],'expected_P':p[j]})
       for site,target,baseline,rv,valid,tn,bn,rn in siteinfo:qa_sites.append({'context':ctx,'pair_id':edges.iloc[j].pair_id,'donor':donor,'site':site,'ligand_target_CP10k':target[j],'ligand_reference_CP10k':baseline[j],'receptor_CP10k':rv[j],'feature_available':bool(valid[j]),'target_cells':tn,'reference_cells':bn,'receiver_cells':rn})
     status('PATIENT_LR_TEST',layer=layer,dataset=ds,sender=cls,receiver=receiver,patients=n,pairs=len(edges))
   if sourceframes:frames.append(fdr_by_family(pd.concat(sourceframes,ignore_index=True)))
   save(pd.DataFrame(registry),'LR_contrast_availability.tsv')
 allstats=pd.concat(frames,ignore_index=True);save(allstats,'all_patient_LR_compatibility_statistics.tsv.gz');save(allstats[allstats.nominal],'nominal_patient_LR_compatibility_statistics.tsv.gz');save(pd.DataFrame(qa),'independent_LR_patient_P_QA_input.tsv.gz');save(pd.DataFrame(qa_sites),'independent_LR_site_score_reconstruction_input.tsv.gz');writej(O/'resources/LR_EFFECT_CONTEXTS.json',contexts)
 q={'at':datetime.datetime.now().astimezone().isoformat(),'samples_processed':len(m),'single_gene_role_validated_pairs':len(edges),'complex_or_unmapped_pairs_unassessed':len(unassessed),'statistical_rows':len(allstats),'valid_tests':int(allstats.P.notna().sum()),'nominal_rows':int(allstats.P.lt(.05).sum()),'FDR_rows':int(allstats.FDR.lt(.05).sum()),'source_scoped_independent_patients':True,'unknown_cells_not_forced_to_tumor_or_T':True,'RNA_compatibility_not_actual_communication':True,'C2_final_synthesis_pending':True}
 writej(O/'delivery/LR_EXPRESSION_QA.json',q);(O/'COMPLETE').write_text('SUCCESS\n');status('complete',**q);print('STAGE_D_LR_EXPRESSION_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();status('failed');(O/'FAILED').write_text(traceback.format_exc())
finally:(O/'LR_expression.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)