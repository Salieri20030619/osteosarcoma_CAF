from common import *
import anndata as ad,gc
from scipy import sparse
assert (O/'FREEZE_COMPLETE').is_file() and (O/'freeze.exit_code').read_text().strip()=='0'
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad');ob=read(O/'tables/frozen_expanded_cells.tsv.gz').set_index('cell_id').loc[a.obs_names];features=a.var_names.astype(str).to_numpy();presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict();manifest=read(A/'tables/clinical_sample_registry.tsv')
cw=read(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',keep_default_na=False).set_index('feature').loc[features];symbols=cw.canonical_symbol.to_numpy();quality=cw.mapping_method.ne('HGNC_Ensembl_approved').astype(int).to_numpy();sets=json.loads((O/'resources/programs_frozen.json').read_text())['sets'];programs=list(sets)
oldprogram=read(C/'tables/all_patient_program_statistics.tsv.gz');oldprogram=oldprogram[oldprogram.layer.eq('expanded')];registry=[];summaries=[];progs=[];ledgers=[];reuse=[];canonical=[]
for ds in manifest.dataset.unique():
 samples=manifest.loc[manifest.dataset.eq(ds),'sample'].tolist();gi=np.flatnonzero(np.asarray(presence[[sidx[s] for s in samples]].sum(0)).ravel()>0);fe=features[gi];sy=symbols[gi];q=quality[gi];g=ob[ob.dataset.eq(ds)];specs=[{'contrast':cl+'_vs_rest','target':cl,'parent':'','kind':'main'} for cl in ['F1','F2','F3','F4','F5']]
 for parent in ['F1','F2','F3','F4','F5']:
  for state in sorted(ob.loc[ob.class_id.eq(parent),'round1_state'].unique()):
   specs.append({'contrast':state+'_vs_parent_rest','target':state,'parent':parent,'kind':'internal'})
 mapping={}
 for j,(f,s,pr) in enumerate(zip(fe,sy,q)):
  if s and (s not in mapping or (pr,f)<mapping[s][0]):mapping[s]=((pr,f),j)
 canonical.extend({'dataset':ds,'symbol':s,'feature':fe[v[1]],'priority':v[0][0]} for s,v in mapping.items())
 mem={name:[mapping[s][1] for s in members if s in mapping] for name,members in sets.items()};gf=[];pf=[]
 oldfile=C/'tables'/f'expanded__{ds}__full_gene_statistics.tsv.gz';old=read(oldfile) if oldfile.exists() else pd.DataFrame()
 for spec in specs:
  name=spec['contrast'];kind=spec['kind'];parent=spec['parent'];sub=g if kind=='main' else g[g.class_id.eq(parent)];targetmask=sub.class_id.eq(spec['target']) if kind=='main' else sub.round1_state.eq(spec['target']);targetcells=int(targetmask.sum());refcells=int((~targetmask).sum())
  if kind=='internal' and targetcells==0:
   registry.append({'dataset':ds,**spec,'N':0,'target_cells':0,'reference_cells':refcells,'features':len(fe),'status':'state_absent_from_source','effect_file':''});continue
  ctx=f'expanded__{ds}__{name}';oldnpz=C/'objects'/f'{ctx}__patient_effects.npz';oldrows=old[old.contrast.eq(name)] if len(old) else pd.DataFrame();calc='new';donors=[]
  if oldnpz.exists() and len(oldrows):
   z=np.load(oldnpz);assert np.array_equal(z['feature_indices'],gi);d=z['effects'];donors=z['donors'].astype(str).tolist();assert oldrows.feature.tolist()==fe.tolist();pv=oldrows.P.to_numpy();N=oldrows.N.to_numpy();method=oldrows.test.to_numpy();assert np.allclose(mean_valid(d),oldrows.mean_delta,equal_nan=True,atol=1e-10);calc='reused_phase43_exact_same_contrast'
  else:
   rows=[]
   for donor,dg in sub.groupby('donor_record',observed=True):
    siteeffects=[]
    for site,sg in dg.groupby('site',observed=True):
     lm=sg.class_id.eq(spec['target']).to_numpy() if kind=='main' else sg.round1_state.eq(spec['target']).to_numpy();rm=~lm
     if not lm.any() or not rm.any():continue
     il=ob.index.get_indexer(sg.index[lm]);ir=ob.index.get_indexer(sg.index[rm]);ll=sg.loc[sg.index[lm],'full_library_counts'].sum();rr=sg.loc[sg.index[rm],'full_library_counts'].sum()
     if ll<=0 or rr<=0:continue
     lc=np.asarray(a.X[il][:,gi].sum(0)).ravel();rc=np.asarray(a.X[ir][:,gi].sum(0)).ravel();required=sg['sample'].astype(str).unique();valid=np.asarray(presence[[sidx[s] for s in required]][:,gi].sum(0)).ravel()==len(required);effect=np.log2(lc/ll*1e4+.1)-np.log2(rc/rr*1e4+.1);effect[~valid]=np.nan;siteeffects.append(effect)
     ledgers.append({'dataset':ds,'contrast':name,'kind':kind,'parent':parent,'donor_record':donor,'site':site,'target_cells':len(il),'reference_cells':len(ir),'sample_ids':';'.join(required),'valid_features':int(valid.sum())})
    if siteeffects:rows.append(mean_valid(np.vstack(siteeffects)));donors.append(str(donor))
   d=np.vstack(rows) if rows else np.empty((0,len(gi)));pv,N,method=tests(d)
  n=len(donors);family=f'expanded|{ds}|'+('all_F1_F5_genes' if kind=='main' else parent+'_internal_genes')
  with warnings.catch_warnings():warnings.simplefilter('ignore',RuntimeWarning);median=np.nanmedian(d,axis=0) if n else np.full(len(gi),np.nan)
  frame=pd.DataFrame({'layer':'expanded','dataset':ds,'contrast':name,'class_id':spec['target'] if kind=='main' else parent,'kind':kind,'parent':parent,'feature':fe,'symbol':sy,'N':N,'mean_delta':mean_valid(d),'median_delta':median,'positive_patients':np.sum(d>0,axis=0),'negative_patients':np.sum(d<0,axis=0),'P':pv,'test':method,'family':family,'calculation':calc})
  if kind=='main' or n:gf.append(frame)
  np.savez_compressed(O/'objects'/f'{ctx}__patient_effects.npz',effects=d,feature_indices=gi,donors=np.asarray(donors,dtype=str))
  score=np.column_stack([mean_valid(d[:,ix],axis=1) if ix else np.full(n,np.nan) for ix in mem.values()]);nn=np.isfinite(score).sum(0);ppr=np.full(len(programs),np.nan);mm=np.full(len(programs),'not_tested',dtype=object);need=np.ones(len(programs),dtype=bool);op=oldprogram[(oldprogram.dataset==ds)&(oldprogram.contrast==name)].set_index('program');preused=0
  for j,namep in enumerate(programs):
   if namep in op.index:
    prior=op.loc[namep];same=np.isclose(mean_valid(score[:,j]),prior.mean_delta,equal_nan=True,atol=1e-10) and int(nn[j])==int(prior.N)
    if same:ppr[j]=prior.P;mm[j]=prior.test;need[j]=False;preused+=1
  if need.any():ppr[need],nn[need],mm[need]=tests(score[:,need])
  ps=pd.DataFrame({'layer':'expanded','dataset':ds,'contrast':name,'class_id':spec['target'] if kind=='main' else parent,'kind':kind,'parent':parent,'program':programs,'definition_genes':[len(sets[k]) for k in programs],'measured_genes':[len(mem[k]) for k in programs],'missing_genes':[';'.join(s for s in sets[k] if s not in mapping) for k in programs],'N':nn,'mean_delta':mean_valid(score),'positive_patients':np.sum(score>0,axis=0),'negative_patients':np.sum(score<0,axis=0),'P':ppr,'test':mm,'family':family.replace('_genes','_programs')})
  if kind=='main' or n:pf.append(ps)
  np.savez_compressed(O/'objects'/f'{ctx}__program_effects.npz',effects=score,programs=np.asarray(programs),donors=np.asarray(donors,dtype=str))
  registry.append({'dataset':ds,**spec,'N':n,'target_cells':targetcells,'reference_cells':refcells,'features':len(fe),'status':'inference_available' if n>=2 else ('direction_only' if n==1 else 'no_matched_patients'),'effect_file':str(O/'objects'/f'{ctx}__patient_effects.npz')});reuse.append({'dataset':ds,'contrast':name,'gene_calculation':calc,'program_rows_reused':preused,'programs':len(programs)})
  print(json.dumps({'dataset':ds,'contrast':name,'N':n,'reuse':calc}),flush=True)
 full=fdr(pd.concat(gf,ignore_index=True));full.to_csv(O/'tables/genes'/f'expanded__{ds}__all_statistics.tsv.gz',sep='\t',index=False)
 for (contrast,kind,cl),z in full.groupby(['contrast','kind','class_id'],dropna=False):summaries.append({'dataset':ds,'contrast':contrast,'kind':kind,'class_id':cl,'gene_features':len(z),'valid_tests':int(z.P.notna().sum()),'N_max':int(z.N.max()),'nominal_up':int((z.nominal&z.mean_delta.gt(0)).sum()),'nominal_down':int((z.nominal&z.mean_delta.lt(0)).sum()),'FDR_up':int((z.FDR_significant&z.mean_delta.gt(0)).sum()),'FDR_down':int((z.FDR_significant&z.mean_delta.lt(0)).sum())})
 prog=fdr(pd.concat(pf,ignore_index=True));prog.to_csv(O/'tables/programs'/f'expanded__{ds}__all_program_statistics.tsv.gz',sep='\t',index=False);progs.append(prog);pd.DataFrame(registry).to_csv(O/'tables/executed_contrast_registry.tsv',sep='\t',index=False);del full,old,gf,pf;gc.collect()
pd.concat(progs,ignore_index=True).to_csv(O/'tables/all_program_statistics.tsv.gz',sep='\t',index=False);pd.DataFrame(summaries).to_csv(O/'tables/gene_statistics_summary.tsv',sep='\t',index=False);pd.DataFrame(ledgers).to_csv(O/'tables/new_matched_patient_site_roster.tsv',sep='\t',index=False);pd.DataFrame(reuse).to_csv(O/'tables/reuse_ledger.tsv',sep='\t',index=False);pd.DataFrame(canonical).to_csv(O/'resources/source_canonical_feature_map.tsv',sep='\t',index=False)
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cells':a.n_obs,'features':a.n_vars,'main_planned_contrasts':45,'registered_contexts':len(registry),'evaluated_gene_contexts':len(summaries),'source_program_files':len(progs),'source_specific_feature_N':True,'no_new_cell_filters':True,'new_five_class_FDR_families':True};dump(O/'delivery/STATISTICS_ACCEPTANCE.json',qa);(O/'STATISTICS_COMPLETE').write_text('SUCCESS\n');print('STATISTICS_COMPLETE',flush=True)