# SYNTAX_CHECK: python -m py_compile stageC_full_gene_statistics.py
# RUN_COMMAND: batch_bench/bin/python -u stageC_full_gene_statistics.py
# EXPECTED_SUCCESS_MARKER: STAGE_C1_FULL_GENE_STATISTICS_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,traceback,itertools,collections,warnings
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';B=P/'phase42_F3F4_identity_states_20261001';O=P/'phase43_full_transcriptome_discovery_20261001'
for sub in ['code','logs','tables','objects','figures','delivery','resources']:(O/sub).mkdir(parents=True,exist_ok=True)
def writej(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS_C1.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def mean_valid(d,axis=0):
 valid=np.isfinite(d);num=np.where(valid,d,0).sum(axis);den=valid.sum(axis);out=np.full(np.shape(num),np.nan,dtype=float);np.divide(num,den,out=out,where=den>0);return out
signs_cache={}
def tests(d):
 p=np.full(d.shape[1],np.nan);N=np.isfinite(d).sum(0);method=np.full(d.shape[1],'not_tested_insufficient_replication',dtype=object)
 if not d.shape[0]:return p,N,method
 patterns,inv=np.unique(np.isfinite(d).T,axis=0,return_inverse=True)
 for k,pat in enumerate(patterns):
  n=int(pat.sum());ix=np.flatnonzero(inv==k)
  if n<2:continue
  z=d[np.ix_(pat,ix)];nonzero=np.any(z!=0,axis=0);method[ix[~nonzero]]='not_tested_all_zero_effect';ix=ix[nonzero];z=z[:,nonzero]
  if not len(ix):continue
  if n not in signs_cache:
   if n<=15:signs_cache[n]=np.array(list(itertools.product([-1,1],repeat=n)),dtype=float)
   else:signs_cache[n]=np.random.default_rng(20261001).choice([-1.,1.],size=(19999,n))
  sg=signs_cache[n];obs=np.abs(z.mean(0));label='exact_patient_signflip' if n<=15 else 'patient_signflip_MC19999_plus_observed'
  for start in range(0,len(ix),256):
   null=np.abs(sg@z[:,start:start+256]/n);hits=(null>=obs[start:start+256][None,:]-1e-12).sum(0);p[ix[start:start+256]]=hits/len(sg) if n<=15 else (hits+1)/(len(sg)+1)
  method[ix]=label
 return p,N,method
def fdr_by_family(d):
 d=d.copy();d['FDR']=np.nan
 for family,g in d.groupby('family',sort=False):
  q=multipletests(g.P.fillna(1).to_numpy(),method='fdr_bh')[1];d.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 d['nominal']=d.P.lt(.05);d['FDR_significant']=d.FDR.lt(.05);return d
rc=1
try:
 assert (A/'setup.exit_code').read_text().strip()=='0' and (A/'COMPLETE').read_text().strip()=='SUCCESS'
 assert (B/'identity_states.exit_code').read_text().strip()=='0' and (B/'COMPLETE').read_text().strip()=='SUCCESS'
 a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad');ob=a.obs.copy();presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=pd.read_csv(A/'tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample').presence_row.to_dict()
 manifest=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);features=a.var_names.astype(str).to_numpy();rawsym=a.var.symbol.astype(str).to_numpy()
 hgnc=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);hgnc=hgnc[hgnc.status.eq('Approved')]
 enmap={};aliases=collections.defaultdict(set);approved=set(hgnc.symbol)
 for row in hgnc.itertuples(index=False):
  if row.ensembl_gene_id:enmap[row.ensembl_gene_id]=row.symbol
  for field in ['prev_symbol','alias_symbol']:
   for alias in getattr(row,field,'').split('|'):
    if alias:aliases[alias].add(row.symbol)
 symbols=[];mapping=[]
 for en,sy in zip(features,rawsym):
  base=en.split('.')[0]
  if base in enmap:canon=enmap[base];how='HGNC_Ensembl_approved'
  elif sy in approved:canon=sy;how='HGNC_approved_symbol'
  elif sy in aliases and len(aliases[sy])==1:canon=next(iter(aliases[sy]));how='HGNC_unique_previous_or_alias'
  else:canon=sy;how='unresolved_symbol_kept' if sy else 'unmapped'
  symbols.append(canon);mapping.append({'feature':en,'source_consensus_symbol':sy,'canonical_symbol':canon,'mapping_method':how})
 symbols=np.asarray(symbols);save(pd.DataFrame(mapping),'HGNC_full_gene_symbol_crosswalk.tsv.gz')
 states=pd.read_csv(B/'tables/all_F3F4_internal_cell_states.tsv.gz',sep='\t',keep_default_na=False);stproto=json.loads((B/'resources/INTERNAL_STATE_PROTOCOL_FROZEN.json').read_text());assert hashlib.sha256((B/'tables/all_F3F4_internal_cell_states.tsv.gz').read_bytes()).hexdigest()==stproto['mapping_sha256']
 sets={}
 for lib,fname in [('hallmark','hallmark.gmt'),('matrisome','curated_cgp.gmt'),('reactome','reactome.gmt')]:
  for line in (A/'resources'/fname).read_text().splitlines():
   items=line.split('\t');name=items[0]
   if lib=='matrisome' and not name.startswith('NABA_'):continue
   genes=list(dict.fromkeys((s if s in approved else (next(iter(aliases[s])) if s in aliases and len(aliases[s])==1 else s)) for s in items[2:]));sets[f'{lib}|{name}']=genes
 assert any(x.startswith('matrisome|') for x in sets)
 writej(O/'resources/FULL_GENE_STATISTICS_PROTOCOL_FROZEN.json',{'state_mapping_sha256':stproto['mapping_sha256'],'class_contrasts':['F3_vs_rest','F4_vs_rest','F3_vs_F4','F3_vs_F2','F4_vs_F2','F2_vs_rest'],'internal_contrasts':'each frozen internal state versus rest of same parent class','counts':'full raw count features for each source; no shared universe restriction','missing_features':'sample-presence array determines gene-specific valid site/patient N; never treat absent feature as zero expression','normalization':'full-library class pseudobulk CP10k + 0.1; log2 difference; equal matched-site weighting within patient','test':'two-sided patient mean sign flip, exact N<=15, 19999 MC+observed beyond15; all-zero and N<2 NA','gene_family':'layer|source|class gene contrasts or internal gene contrasts per parent','program_family':'layer|source|class program contrasts or internal program contrasts per parent','program_definition':'mean patient pseudobulk log2 difference over canonical gene-set members available for that patient; missing members recorded, no coverage gate','pathway_symbol_collapse':'prefer HGNC approved Ensembl feature; then lexicographic feature ID, never choose by P/effect','program_sets':len(sets),'genes_selected_by':'valid raw P<0.05, both directions; FDR reported not entry gate','new_cell_entry_filters':False,'discovery_label_selection_not_corrected_by_FDR':True,'no_global_patient_pooling':True})
 all_summary=[];all_program=[];coverage=[];registry=[];contexts=[];donors_all=[]
 baseplanned=[('F3_vs_rest','F3','rest'),('F4_vs_rest','F4','rest'),('F3_vs_F4','F3','F4'),('F3_vs_F2','F3','F2'),('F4_vs_F2','F4','F2'),('F2_vs_rest','F2','rest')]
 for layer in ['supported','expanded']:
  label='supported_class_id' if layer=='supported' else 'expanded_class_id';layerob=ob[ob.in_supported.astype(bool)] if layer=='supported' else ob
  statemap=states[states.layer.eq(layer)].set_index('cell_id').state_id.to_dict()
  planned=[{'contrast':n,'target':t,'reference':r,'kind':'class','parent':''} for n,t,r in baseplanned]
  for parent in ['F3','F4']:
   for st in sorted(states[(states.layer==layer)&(states.parent_class==parent)].state_id.unique()):planned.append({'contrast':st+'_vs_parent_rest','target':st,'reference':'parent_rest','kind':'internal','parent':parent})
  for ds in manifest.dataset.unique():
   g=layerob[layerob.dataset.eq(ds)];samples=list(manifest.loc[manifest.dataset.eq(ds),'sample']);gi=np.flatnonzero(np.asarray(presence[[sidx[s] for s in samples]].sum(0)).ravel()>0);fe=features[gi];sy=symbols[gi];pref={}
   # Deterministic mapping of features to one canonical pathway symbol, before any P values.
   for j,(en,symbol) in enumerate(zip(fe,sy)):
    if not symbol:continue
    quality=0 if en.split('.')[0] in enmap and enmap[en.split('.')[0]]==symbol else 1
    key=(quality,en)
    if symbol not in pref or key<pref[symbol][0]:pref[symbol]=(key,j)
   symbol_indices={s:v[1] for s,v in pref.items()};gf=[];pf=[]
   for spec in planned:
    name=spec['contrast'];sub=g[g[label].eq(spec['parent'])] if spec['kind']=='internal' else g
    rows=[];patient_names=[];site_records=[]
    for donor,dg in sub.groupby('donor_record',observed=True):
     deltas=[];used=[]
     for site,z in dg.groupby('site',observed=True):
      if spec['kind']=='internal':left=np.asarray([statemap.get(cid,'')==spec['target'] for cid in z.index]);right=~left
      else:left=z[label].astype(str).eq(spec['target']).to_numpy();right=(~left) if spec['reference']=='rest' else z[label].astype(str).eq(spec['reference']).to_numpy()
      if not left.any() or not right.any():continue
      zzleft=z.index[left];zzright=z.index[right];il=ob.index.get_indexer(zzleft);ir=ob.index.get_indexer(zzright);assert (il>=0).all() and (ir>=0).all()
      lc=np.asarray(a.X[il][:,gi].sum(0)).ravel();rr=np.asarray(a.X[ir][:,gi].sum(0)).ravel();ll=pd.to_numeric(z.loc[zzleft,'full_library_counts']).sum();rl=pd.to_numeric(z.loc[zzright,'full_library_counts']).sum();assert ll>0 and rl>0
      required=list(z.loc[zzleft.union(zzright),'sample'].astype(str).unique());valid=np.asarray(presence[[sidx[s] for s in required]][:,gi].sum(0)).ravel()==len(required)
      delta=np.log2(lc/ll*1e4+.1)-np.log2(rr/rl*1e4+.1);delta[~valid]=np.nan;deltas.append(delta);used.append(str(site))
      site_records.append({'layer':layer,'dataset':ds,'contrast':name,'donor_record':donor,'site':site,'target_cells':len(il),'reference_cells':len(ir),'source_samples':';'.join(required),'valid_features':int(valid.sum())})
     if deltas:rows.append(mean_valid(np.vstack(deltas)));patient_names.append(str(donor))
    n=len(rows);family=f'{layer}|{ds}|'+('class_genes' if spec['kind']=='class' else spec['parent']+'_internal_genes');program_family=family.replace('_genes','_programs')
    registry.append({'layer':layer,'dataset':ds,**spec,'N_matched_patients':n,'gene_features':len(gi),'gene_family':family,'program_family':program_family,'state':'inference_available' if n>=2 else ('direction_only' if n==1 else 'not_evaluable_no_matched_patients')})
    if not n:continue
    status('FULL_GENE_TEST',layer=layer,dataset=ds,contrast=name,patients=n,features=len(gi))
    d=np.vstack(rows);p,N,method=tests(d);mean=mean_valid(d)
    with warnings.catch_warnings():warnings.simplefilter('ignore',RuntimeWarning);median=np.nanmedian(d,axis=0)
    frame=pd.DataFrame({'layer':layer,'dataset':ds,'contrast':name,'kind':spec['kind'],'parent':spec['parent'],'feature':fe,'symbol':sy,'N':N,'median_delta':median,'mean_delta':mean,'positive_patients':np.sum(d>0,axis=0),'negative_patients':np.sum(d<0,axis=0),'P':p,'test':method,'family':family});gf.append(frame)
    safe=f'{layer}__{ds}__{name}';np.savez_compressed(O/'objects'/f'{safe}__patient_effects.npz',effects=d,feature_indices=gi,donors=np.asarray(patient_names,dtype=str));contexts.append({'layer':layer,'dataset':ds,**spec,'N_matched_patients':n,'effect_file':f'objects/{safe}__patient_effects.npz','gene_result_file':f'tables/{layer}__{ds}__full_gene_statistics.tsv.gz','signature_symbol_policy':'deterministic_HGNC_feature_before_P','context_id':safe});donors_all.extend(site_records)
    vals=[];setmeta=[]
    for term,members in sets.items():
     measured=[symbol_indices[s] for s in members if s in symbol_indices]
     if measured:v=mean_valid(d[:,measured],axis=1)
     else:v=np.full(n,np.nan)
     vals.append(v);setmeta.append({'layer':layer,'dataset':ds,'contrast':name,'program':term,'definition_genes':len(members),'measured_genes':len(measured),'missing_genes':';'.join(s for s in members if s not in symbol_indices),'family':program_family})
    pdiff=np.column_stack(vals);pv,pn,pm=tests(pdiff);pr=pd.DataFrame(setmeta);pr['N']=pn;pr['mean_delta']=mean_valid(pdiff);pr['positive_patients']=np.sum(pdiff>0,axis=0);pr['negative_patients']=np.sum(pdiff<0,axis=0);pr['P']=pv;pr['test']=pm;pf.append(pr)
    np.savez_compressed(O/'objects'/f'{safe}__program_patient_effects.npz',effects=pdiff,programs=np.asarray(list(sets),dtype=str),donors=np.asarray(patient_names,dtype=str))
   if gf:
    full=fdr_by_family(pd.concat(gf,ignore_index=True));save(full,f'{layer}__{ds}__full_gene_statistics.tsv.gz')
    for (contrast,kind,parent),z in full.groupby(['contrast','kind','parent'],dropna=False):all_summary.append({'layer':layer,'dataset':ds,'contrast':contrast,'kind':kind,'parent':parent,'gene_features':len(z),'valid_tests':int(z.P.notna().sum()),'N_max':int(z.N.max()),'nominal_up':int((z.nominal&z.median_delta.gt(0)).sum()),'nominal_down':int((z.nominal&z.median_delta.lt(0)).sum()),'FDR_up':int((z.FDR_significant&z.median_delta.gt(0)).sum()),'FDR_down':int((z.FDR_significant&z.median_delta.lt(0)).sum())})
    # All nominal hits are kept; no top-N limitation in discovery.
    save(full[full.nominal],f'{layer}__{ds}__nominal_genes.tsv.gz')
   if pf:all_program.append(fdr_by_family(pd.concat(pf,ignore_index=True)))
   save(pd.DataFrame(registry),'executed_contrast_registry.tsv');save(pd.DataFrame(all_summary),'full_gene_statistics_summary.tsv');writej(O/'resources/EFFECT_CONTEXTS.json',contexts)
 if all_program:save(pd.concat(all_program,ignore_index=True),'all_patient_program_statistics.tsv.gz')
 save(pd.DataFrame(donors_all),'matched_patient_site_roster.tsv');save(pd.DataFrame(registry),'executed_contrast_registry.tsv');save(pd.DataFrame(all_summary),'full_gene_statistics_summary.tsv')
 qa={'at':datetime.datetime.now().astimezone().isoformat(),'registered_contrasts':len(registry),'evaluated_direction_contexts':len(contexts),'inference_contexts':sum(x['N_matched_patients']>=2 for x in contexts),'full_union_features':a.n_vars,'symbol_mappings':len(mapping),'program_definitions':len(sets),'frozen_cells':a.n_obs,'missing_genes_not_zero':True,'gene_specific_N_recorded':True,'all_zero_P_NA':True,'input_stageA_accepted':True,'input_stageB_accepted':True,'all_nominal_genes_retained':True,'cross_source_patients_not_pooled':True,'inference_unit':'independent source-scoped patient','enrichment_not_completed_by_this_marker':True}
 writej(O/'delivery/FULL_GENE_STATISTICS_QA.json',qa);(O/'C1_COMPLETE').write_text('SUCCESS\n');status('complete',**qa);print('STAGE_C1_FULL_GENE_STATISTICS_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();status('failed');(O/'C1_FAILED').write_text(traceback.format_exc())
finally:(O/'full_gene_statistics.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)