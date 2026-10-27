# SYNTAX_CHECK: python -m py_compile stageA_setup.py
# RUN_COMMAND: batch_bench/bin/python -u stageA_setup.py
# EXPECTED_SUCCESS_MARKER: STAGE_A_FULL_TRANSCRIPTOME_SETUP_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,traceback,collections
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=P/'phase26_client_reanalysis_20260922/inputs';T=P/'phase30_two_layer_classification_20260926/tables';O=P/'phase41_full_transcriptome_setup_20261001'
for sub in ['resources','logs','code','tables','objects','delivery']:(O/sub).mkdir(parents=True,exist_ok=True)
def writej(f,x):
 t=f.with_suffix(f.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):
 writej(O/'STATUS.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
rc=1
try:
 manifest=pd.read_csv(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',sep='\t',keep_default_na=False)
 assert len(manifest)==107 and manifest['sample'].is_unique
 exp=pd.read_csv(T/'expanded_classified_cells.tsv.gz',sep='\t',keep_default_na=False).set_index('cell_id');sup=pd.read_csv(T/'supported_classified_cells.tsv.gz',sep='\t',keep_default_na=False).set_index('cell_id')
 assert len(exp)==8174 and len(sup)==1224 and exp.index.is_unique and sup.index.is_unique and set(sup.index).issubset(exp.index)
 exp['expanded_class_id']=exp.class_id;exp['supported_class_id']=sup.class_id.reindex(exp.index).fillna('');exp['supported_cluster']=sup.cluster.reindex(exp.index).fillna('');exp['in_supported']=exp.index.isin(sup.index)
 clinical=pd.read_csv(P/'phase31_customer_completion_20260926/resources/clinical_mapping_frozen.tsv',sep='\t',keep_default_na=False)
 assert clinical['sample'].is_unique and set(manifest['sample'])==set(clinical['sample'])
 save(clinical,'clinical_sample_registry.tsv')
 # Keep each observed feature universe and its original symbol mappings.
 universe={};sampleinfo=[];sample_genes={};feature_symbols=collections.defaultdict(collections.Counter);cross=[];sample_lib=[];selected=[];obs=[];allids=set();candidate_checks=[]
 for no,row in enumerate(manifest.itertuples(index=False),1):
  sm=row.sample;f=I/sm/'full_called_counts.h5ad';assert f.exists(),str(f)
  a=ad.read_h5ad(f,backed='r');ids=a.var_names.astype(str).to_numpy();syms=a.var['symbol'].fillna('').astype(str).to_numpy();assert len(set(ids))==len(ids)
  genehash=hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest();sample_genes[sm]=ids
  if genehash not in universe:universe[genehash]={'sample':sm,'ids':ids,'symbols':syms}
  for en,sym in zip(ids,syms):
   allids.add(en)
   if sym and sym not in ['nan','NA']:feature_symbols[en][sym]+=1
  cells=pd.read_csv(I/sm/'cells.tsv',sep='\t',keep_default_na=False);assert cells.cell_id.is_unique and set(cells.cell_id)==set(a.obs_names)
  g=exp[exp['sample'].eq(sm)].copy();ix=a.obs_names.get_indexer(g.index);assert (ix>=0).all(),sm
  if len(g):
   order=np.argsort(ix);ix=ix[order];g=g.iloc[order];x=a.X[ix,:].tocsr().astype(np.float64)
   assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data)),sm
   lib=np.asarray(x.sum(1)).ravel();recorded=pd.to_numeric(g.full_library_counts).to_numpy()
   err=float(np.max(np.abs(lib-recorded)));assert np.allclose(lib,recorded,rtol=1e-6,atol=.01),(sm,err)
   g['full_library_counts']=lib;selected.append((sm,ids,x));obs.append(g)
   candidate_checks.append({'sample':sm,'candidate_cells':len(g),'max_full_library_count_error':err,'integer_counts':True})
  q=json.loads((I/sm/'input_QA.json').read_text())
  sampleinfo.append({'sample':sm,'dataset':row.dataset,'called_cells':a.n_obs,'source_features':a.n_vars,'unique_symbols':len(set(syms)-{'','nan','NA'}),'duplicate_symbol_rows':int(pd.Series(syms).duplicated().sum()),'candidate_cells':len(g),'candidate_supported':int(g.in_supported.sum()),'gene_universe_sha256':genehash,'input_path':str(f),'input_bytes':f.stat().st_size,'source_sha256_recorded':q.get('source_sha256',''),'counts_sha256_recorded':q.get('counts_sha256',''),'called_cell_library_counts_total':pd.to_numeric(cells.n_counts).sum(),'no_new_cell_filters':True})
  a.file.close();status('INVENTORY',samples_done=no,total_samples=len(manifest),frozen_cells_matched=sum(len(z) for z in obs))
 assert sum(len(z) for z in obs)==len(exp)
 union=sorted(allids);uidx={x:i for i,x in enumerate(union)};symbols=[feature_symbols[x].most_common(1)[0][0] if feature_symbols[x] else '' for x in union]
 var=pd.DataFrame({'symbol':symbols,'symbol_ambiguous':[len(feature_symbols[x])>1 for x in union]},index=pd.Index(union,name='feature_id'))
 blocks=[]
 for sm,ids,x in selected:
  co=x.tocoo();mapped=np.asarray([uidx[v] for v in ids],dtype=np.int64);blocks.append(sparse.csr_matrix((co.data,(co.row,mapped[co.col])),shape=(x.shape[0],len(union))))
 X=sparse.vstack(blocks,format='csr');ob=pd.concat(obs);assert set(ob.index)==set(exp.index)
 ob.index.name='cell_id';a=ad.AnnData(X=X,obs=ob,var=var);a.uns['source_cell_ids_frozen']=True;a.uns['missing_features_not_unexpressed']=True
 a.write_h5ad(O/'objects/all_candidate_full_genes.h5ad',compression='gzip')
 presence=sparse.lil_matrix((len(manifest),len(union)),dtype=np.int8)
 for no,row in enumerate(manifest.itertuples(index=False)):presence[no,[uidx[x] for x in sample_genes[row.sample]]]=1
 sparse.save_npz(O/'objects/sample_feature_presence.npz',presence.tocsr())
 save(pd.DataFrame({'sample':manifest['sample'],'presence_row':np.arange(len(manifest))}),'feature_presence_sample_index.tsv')
 save(pd.DataFrame({'feature':union,'symbol':symbols,'presence_column':np.arange(len(union)),'symbol_ambiguous':var.symbol_ambiguous.to_numpy()}),'feature_presence_gene_index.tsv.gz')
 for h,u in universe.items():
  cross.extend({'universe_sha256':h,'representative_sample':u['sample'],'feature':en,'source_symbol':sy,'consensus_symbol':symbols[uidx[en]],'symbol_ambiguous':bool(var.loc[en,'symbol_ambiguous'])} for en,sy in zip(u['ids'],u['symbols']))
 save(pd.DataFrame(cross),'gene_id_crosswalk.tsv.gz');save(pd.DataFrame(sampleinfo),'source_sample_gene_coverage.tsv');save(pd.DataFrame(candidate_checks),'frozen_candidate_count_alignment_QA.tsv')
 coverage=[]
 for ds,m in manifest.groupby('dataset'):
  sets=[set(sample_genes[x]) for x in m['sample']];inter=set.intersection(*sets);uni=set.union(*sets);raw=pd.DataFrame(sampleinfo);raw=raw[raw.dataset.eq(ds)]
  coverage.append(dict(dataset=ds,samples=len(m),called_cells=int(raw.called_cells.sum()),min_source_features=int(raw.source_features.min()),max_source_features=int(raw.source_features.max()),within_source_all_sample_common=len(inter),within_source_union=len(uni),candidate_expanded=int(raw.candidate_cells.sum()),candidate_supported=int(raw.candidate_supported.sum())))
 save(pd.DataFrame(coverage),'source_gene_coverage.tsv')
 alias=[]
 for sm,ids in sample_genes.items():
  row=manifest.set_index('sample').loc[sm];u=universe[next(z['gene_universe_sha256'] for z in sampleinfo if z['sample']==sm)];s=set(u['symbols'])
  for name in ['CTGF','CCN2','ANGPT2','CA9','LDHA']:
   alias.append({'sample':sm,'dataset':row.dataset,'symbol':name,'present_in_full_source':name in s,'in_old_shared_symbol_universe':False})
 save(pd.DataFrame(alias),'previously_missing_mechanism_gene_coverage.tsv')
 donors=clinical.groupby(['dataset','donor_record','patient_id'],dropna=False).agg(samples=('sample',lambda x:';'.join(x)),sample_records=('sample','size')).reset_index()
 donors['canonical_patient_identifier']=donors.patient_id.astype(str).str.strip()
 dsby=donors.groupby('canonical_patient_identifier').dataset.agg(lambda x:';'.join(sorted(set(x)))).to_dict()
 donors['datasets_with_same_reported_id']=donors.canonical_patient_identifier.map(dsby)
 donors['cross_source_status']=donors.datasets_with_same_reported_id.map(lambda x:'potential_overlap_same_reported_id_do_not_pool' if ';' in x else 'source_scoped_identity_no_global_deduplication_claim')
 save(donors,'patient_overlap_registry.tsv')
 overlap=P/'phase35_external_validation_prep_20260927/tables/OsteoCAR_discovery_overlap_review.tsv'
 if overlap.exists():save(pd.read_csv(overlap,sep='\t',keep_default_na=False),'prior_external_source_overlap_review.tsv')
 contrasts=[];roster=[]
 planned=[('F3_vs_rest','F3','rest'),('F4_vs_rest','F4','rest'),('F3_vs_F4','F3','F4'),('F3_vs_F2','F3','F2'),('F4_vs_F2','F4','F2'),('F2_vs_rest','F2','rest')]
 for layer in ['supported','expanded']:
  label='supported_class_id' if layer=='supported' else 'expanded_class_id';v=ob[ob.in_supported] if layer=='supported' else ob
  for ds in manifest.dataset.unique():
   g=v[v.dataset.eq(ds)]
   for name,tgt,ref in planned:
    used=[]
    for donor,dg in g.groupby('donor_record'):
     sites=[]
     for site,sg in dg.groupby('site'):
      left=sg[label].eq(tgt);right=~left if ref=='rest' else sg[label].eq(ref)
      if left.any() and right.any():sites.append(str(site))
     if sites:
      used.append(donor);roster.append({'layer':layer,'dataset':ds,'contrast':name,'donor_record':donor,'matched_sites':';'.join(sites)})
    n=len(used);contrasts.append({'layer':layer,'dataset':ds,'contrast':name,'target':tgt,'reference':ref,'N_independent_within_source_patients':n,'status':'evaluable_exploratory' if 2<=n<10 else ('evaluable' if n>=10 else 'insufficient_independent_replication'),'gene_family':f'{layer}|{ds}|all_full_gene_planned_contrasts','effect':'log2(CP10k+0.1) target-minus-reference, sites averaged inside patient','test':'exact_two_sided_patient_sign_flip_N_le_15_else_MC19999_plus_observed','sample_roster':';'.join(used)})
 save(pd.DataFrame(contrasts),'contrast_registry.tsv');save(pd.DataFrame(roster),'contrast_patient_roster.tsv')
 protocol={'at':datetime.datetime.now().astimezone().isoformat(),'phase':'A','plan':'CAF_NEXT_STAGE_FULL_PLAN_20260930_v1.txt','user_order':'execute sequentially; downloads background direct, no proxy','source_scoped_patients':True,'cross_source_pooling':False,'layers_are_not_independent_replicates':True,'frozen_expanded':8174,'frozen_supported':1224,'union_features':len(union),'candidate_genes':'full raw source features, missing features recorded by sample presence; no old global intersection','normalization':'full original per-cell/library counts; class pseudobulk aggregation; within same patient/site; sites equal weight inside patient','count_integrity':'raw nonnegative integer counts; full selected row sums match frozen full_library_counts','new_hypothesis_families':'per layer and source, all six planned gene contrasts; program, ORA-library/direction, TF, communication, clinical and spatial families separate','selection':'target raw P<0.05; all P/FDR retained; FDR<0.05 corrected; N<10 exploratory; n<2 or all-zero NA','class_discovery':'same data labels, tests discovery not independent confirmation','software_filters':'no new cell entry thresholds, no UMI/MT/confidence filters','alias_mapping':'source/consensus symbol crosswalk preserved; HGNC resource normalization applied in a later registered step','downloads':'independent background, trust_env False, no proxy forwarding'}
 writej(O/'resources/PROTOCOL_FROZEN.json',protocol)
 writej(O/'resources/STATISTICAL_FAMILIES.json',{'gene':'layer|dataset|six planned full-gene contrasts','program':'layer|dataset|all planned program contrasts','ORA':'layer|dataset|contrast|library|direction','TF':'layer|dataset|all registered TF contrasts','spatial':'frozen before spatial P values, primary program retained or multiple endpoints family','external':'separate from discovery','NA':'retained; declared family rows may use 1 only for BH accounting, reported NA'})
 qa={'samples':len(manifest),'called_cells':int(pd.DataFrame(sampleinfo).called_cells.sum()),'datasets':len(coverage),'expanded_cells':a.n_obs,'supported_cells':int(ob.in_supported.sum()),'features_union':a.n_vars,'unique_feature_universes':len(universe),'all_frozen_ids_matched':True,'frozen_row_counts_match':True,'raw_integer_counts':True,'new_cell_filters':False,'source_scoped_patient_overlap_registry_present':True,'global_patient_deduplication_complete':False,'contrasts_registered':len(contrasts),'evaluable_contrasts':int((pd.DataFrame(contrasts).N_independent_within_source_patients>=2).sum())}
 writej(O/'delivery/SETUP_QA.json',qa)
 writej(O/'delivery/EVIDENCE_MANIFEST.json',[{'file':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for folder in ['tables','resources','objects'] for f in sorted((O/folder).iterdir()) if f.is_file() and not f.name.endswith('.part')])
 (O/'COMPLETE').write_text('SUCCESS\n');status('complete',**qa);print('STAGE_A_FULL_TRANSCRIPTOME_SETUP_COMPLETE',flush=True);rc=0
except Exception:
 traceback.print_exc();status('failed');(O/'FAILED').write_text(traceback.format_exc())
finally:(O/'setup.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)