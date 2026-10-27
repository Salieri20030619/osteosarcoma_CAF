# SYNTAX_CHECK: python -m py_compile stageD_robustness_tf.py
# RUN_COMMAND: batch_bench/bin/python -u stageD_robustness_tf.py
# EXPECTED_SUCCESS_MARKER: STAGE_D_PATIENT_ROBUSTNESS_TF_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,ast,hashlib,itertools,collections,warnings
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase44_patient_mechanism_20261001'
for sub in ['code','logs','tables','objects','resources','delivery','figures']:(O/sub).mkdir(parents=True,exist_ok=True)
def writej(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS_ROBUSTNESS_TF.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
# Reuse the independently R-verified mathematical functions, not the stage-C main program.
source=C/'code/stageC_full_gene_statistics.py';tree=ast.parse(source.read_text());functions=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in ['mean_valid','tests','fdr_by_family']];signs_cache={};exec(compile(ast.Module(body=functions,type_ignores=[]),str(source),'exec'))
rc=1
try:
 assert (C/'C1_COMPLETE').read_text().strip()=='SUCCESS' and (C/'INDEPENDENT_R_QA_COMPLETE').read_text().strip()=='SUCCESS' and (C/'full_gene_statistics.exit_code').read_text().strip()=='0'
 signatures=json.loads((C/'resources/CODE_SHA256.json').read_text());assert hashlib.sha256(source.read_bytes()).hexdigest()==signatures[source.name]
 contexts=json.loads((C/'resources/EFFECT_CONTEXTS.json').read_text());cross=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);features=cross.feature.to_numpy();symbols=cross.canonical_symbol.to_numpy();qpriority=cross.mapping_method.eq('HGNC_Ensembl_approved').map({True:0,False:1}).to_numpy()
 programs=pd.read_csv(C/'tables/all_patient_program_statistics.tsv.gz',sep='\t');stability=[];program_stability=[];vectors=[];cached=None;nominal=None
 protocol={'input_dependency':'accepted C1 full-gene and patient-program statistics; C2 ongoing, not used as accepted discovery','gene_selection':'all valid raw P<0.05 across registered contrasts; no FDR/positive-direction gate','patient_leaveout':'omit one actually observed patient and recalculate mean effect; direction stability and influence descriptive, no new P, not independent validation','source_reproduction':'source-specific direction and P retained; unknown source overlap not pooled or called independent validation','TF_measure':'signed mean target response over available canonical target genes; inferred regulon response proxy, not direct TF protein activity or causal treatment effect','TF_network':'CollecTRI via downloaded OmniPath with commercial resource query; consensus directed signed edges; conflicting directions unavailable','TF_targets':'all defined regulators, no measured-target coverage threshold; missing members reported','TF_test':'same independent-patient sign-flip test as R-verified C1; all-zero or N<2 P/FDR NA','TF_families':'layer|source|class TF contrasts or parent-specific internal TF contrasts, separate from genes/programs/clinical','cell_and_state_membership':'unchanged'}
 writej(O/'resources/ROBUSTNESS_TF_PROTOCOL_FROZEN.json',protocol)
 def leaveout(d):
  valid=np.isfinite(d);N=valid.sum(0);sm=np.where(valid,d,0).sum(0);base=mean_valid(d);den=N[None,:]-valid
  loo=np.full(d.shape,np.nan);np.divide(sm[None,:]-np.where(valid,d,0),den,out=loo,where=valid&(den>0))
  same=np.sum((np.sign(loo)==np.sign(base)[None,:])&valid,0);fraction=np.full(len(N),np.nan);np.divide(same,N,out=fraction,where=N>=2)
  absden=np.where(valid,np.abs(d),0).sum(0);largest=np.zeros(len(N));np.divide(np.where(valid,np.abs(d),0).max(0),absden,out=largest,where=absden>0)
  with warnings.catch_warnings():warnings.simplefilter('ignore',RuntimeWarning);lo=np.nanmin(loo,axis=0);hi=np.nanmax(loo,axis=0)
  return N,base,fraction,lo,hi,largest
 for ci,ctx in enumerate(contexts,1):
  ef=np.load(C/ctx['effect_file']);gi=ef['feature_indices'];d=ef['effects'];lookup={features[j]:i for i,j in enumerate(gi)}
  fn=ctx['gene_result_file'].replace('__full_gene_statistics.tsv.gz','__nominal_genes.tsv.gz')
  if fn!=cached:nominal=pd.read_csv(C/fn,sep='\t');cached=fn
  g=nominal[nominal.contrast.eq(ctx['contrast'])].copy()
  if len(g):
   col=[lookup[x] for x in g.feature];N,base,fraction,lo,hi,largest=leaveout(d[:,col]);assert np.all(N==g.N.to_numpy())
   g['patient_leaveout_same_direction_fraction']=fraction;g['patient_leaveout_mean_min']=lo;g['patient_leaveout_mean_max']=hi;g['largest_absolute_patient_effect_fraction']=largest;g['source_replication_claim']='within-source discovery; cross-source independence unresolved unless separately documented';stability.append(g)
  pg=programs[(programs.layer==ctx['layer'])&(programs.dataset==ctx['dataset'])&(programs.contrast==ctx['contrast'])&programs.P.lt(.05)].copy()
  pf=C/'objects'/f"{ctx['context_id']}__program_patient_effects.npz"
  if len(pg):
   pe=np.load(pf);ix={s:i for i,s in enumerate(pe['programs'])};N,base,fraction,lo,hi,largest=leaveout(pe['effects'][:,[ix[s] for s in pg.program]])
   pg['patient_leaveout_same_direction_fraction']=fraction;pg['patient_leaveout_mean_min']=lo;pg['patient_leaveout_mean_max']=hi;pg['largest_absolute_patient_effect_fraction']=largest;program_stability.append(pg)
  for row in ctx.items():pass
  status('PATIENT_STABILITY',context=ci,total=len(contexts),nominal_genes=len(g),nominal_programs=len(pg))
 save(pd.concat(stability,ignore_index=True),'all_nominal_gene_patient_stability.tsv.gz');save(pd.concat(program_stability,ignore_index=True),'all_nominal_program_patient_stability.tsv.gz')
 # Cross-source direction summaries are descriptive. Every participating source retains its own N/P/FDR.
 pairs=pd.concat(stability,ignore_index=True)
 rep=pairs.groupby(['layer','contrast','feature','symbol'],dropna=False).agg(nominal_sources=('dataset','nunique'),sources=('dataset',lambda x:';'.join(sorted(set(x)))),minimum_P=('P','min'),minimum_FDR=('FDR','min'),positive_sources=('mean_delta',lambda x:int((x>0).sum())),negative_sources=('mean_delta',lambda x:int((x<0).sum()))).reset_index();rep['independent_external_confirmation']=False;save(rep,'nominal_gene_cross_source_descriptive_summary.tsv.gz')
 network=pd.read_csv(A/'resources/collectri.tsv',sep='\t',keep_default_na=False)
 def flag(v):return v.astype(str).str.lower().isin(['true','1'])
 signed=network[flag(network.consensus_direction)&(flag(network.consensus_stimulation)^flag(network.consensus_inhibition))].copy();signed['weight']=np.where(flag(signed.consensus_stimulation),1.,-1.)
 # Canonicalize network gene names using approved/unique aliases without splitting TF complexes.
 hgnc=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);hgnc=hgnc[hgnc.status.eq('Approved')];approved=set(hgnc.symbol);alias=collections.defaultdict(set)
 for row in hgnc.itertuples(index=False):
  for col in ['prev_symbol','alias_symbol']:
   for x in getattr(row,col,'').split('|'):
    if x:alias[x].add(row.symbol)
 def canon(s):return s if s in approved else (next(iter(alias[s])) if s in alias and len(alias[s])==1 else s)
 signed['regulator']=signed.source_genesymbol.map(canon);signed['target_symbol']=signed.target_genesymbol.map(canon);signed=signed.drop_duplicates(['regulator','target_symbol','weight'])
 conflicted=signed.groupby(['regulator','target_symbol']).weight.nunique();conflicted=set(conflicted[conflicted>1].index);signed=signed[[((r,t) not in conflicted) for r,t in zip(signed.regulator,signed.target_symbol)]]
 save(signed,'signed_TF_target_network_used.tsv.gz');regulons={r:(z.target_symbol.to_numpy(),z.weight.to_numpy()) for r,z in signed.groupby('regulator')};tf_frames=[];coverage=[];qa_rows=[];summary=[]
 for ci,ctx in enumerate(contexts,1):
  ef=np.load(C/ctx['effect_file']);gi=ef['feature_indices'];d=ef['effects'];pref={}
  for k,j in enumerate(gi):
   sy=symbols[j]
   if not sy:continue
   candidate=(qpriority[j],features[j])
   if sy not in pref or candidate<pref[sy][0]:pref[sy]=(candidate,k)
  cols=[];rr=[]
  for regulator,(targets,weights) in regulons.items():
   hit=[k for k,t in enumerate(targets) if t in pref];ii=[pref[targets[k]][1] for k in hit];ww=weights[hit]
   if len(ii):
    z=d[:,ii];valid=np.isfinite(z);num=np.where(valid,z,0)@ww;den=valid@np.abs(ww);v=np.full(len(d),np.nan);np.divide(num,den,out=v,where=den>0)
   else:v=np.full(len(d),np.nan)
   cols.append(v);rr.append({'layer':ctx['layer'],'dataset':ctx['dataset'],'contrast':ctx['contrast'],'kind':ctx['kind'],'parent':ctx['parent'],'regulator':regulator,'defined_targets':len(targets),'measured_targets':len(ii),'missing_targets':';'.join(t for t in targets if t not in pref),'family':f"{ctx['layer']}|{ctx['dataset']}|"+('class_TF' if ctx['kind']=='class' else ctx['parent']+'_internal_TF'),'score_type':'signed_target_response_proxy'})
  matrix=np.column_stack(cols);p,N,method=tests(matrix);frame=pd.DataFrame(rr);frame['N']=N;frame['mean_delta']=mean_valid(matrix);frame['positive_patients']=np.sum(matrix>0,0);frame['negative_patients']=np.sum(matrix<0,0);frame['P']=p;frame['test']=method;tf_frames.append(frame)
  np.savez_compressed(O/'objects'/f"{ctx['context_id']}__TF_patient_scores.npz",scores=matrix,regulators=np.asarray(list(regulons),dtype=str),donors=ef['donors'])
  # Independent-QA export: predefined regulators plus deterministic random targets, with original patient vectors.
  keep=np.flatnonzero(frame.regulator.isin(['SMAD3','SMAD2','TEAD1','JUN','FOS','STAT3','HIF1A','RUNX2','NFKB1','RELA']).to_numpy());extra=np.random.default_rng(20261001).choice(len(frame),size=min(5,len(frame)),replace=False);keep=np.unique(np.r_[keep,extra])
  for j in keep:
   row=frame.iloc[j]
   for donor,value in zip(ef['donors'],matrix[:,j]):qa_rows.append({'context':ctx['context_id'],'regulator':row.regulator,'donor':donor,'score':value,'expected_N':row.N,'expected_mean':row.mean_delta,'expected_P':row.P})
  summary.append({'context_id':ctx['context_id'],'regulators':len(frame),'N_matched_patients':ctx['N_matched_patients'],'valid_TF_tests':int(frame.P.notna().sum()),'nominal_TF_tests':int(frame.P.lt(.05).sum())});status('TF_SIGNED_PROGRAM',context=ci,total=len(contexts),**summary[-1])
 tf=fdr_by_family(pd.concat(tf_frames,ignore_index=True));save(tf,'all_TF_signed_patient_statistics.tsv.gz');save(tf[tf.nominal],'nominal_TF_signed_patient_statistics.tsv');save(pd.DataFrame(summary),'TF_context_summary.tsv');save(pd.DataFrame(qa_rows),'independent_TF_QA_patient_scores.tsv.gz')
 qa={'at':datetime.datetime.now().astimezone().isoformat(),'stability_gene_rows':len(pairs),'stability_program_rows':sum(len(x) for x in program_stability),'TF_regulators':len(regulons),'signed_edges_used':len(signed),'network_rows_input':len(network),'TF_statistics_rows':len(tf),'TF_nominal_rows':int(tf.P.lt(.05).sum()),'TF_FDR_rows':int(tf.FDR.lt(.05).sum()),'independence_claim':'source-scoped patient tests; leaveout is stability, not independent validation','C2_dependency':False,'new_cell_or_target_coverage_filters':False,'causal_TF_activity_not_established':True}
 writej(O/'delivery/ROBUSTNESS_TF_QA.json',qa);(O/'ROBUSTNESS_TF_COMPLETE').write_text('SUCCESS\n');status('complete',**qa);print('STAGE_D_PATIENT_ROBUSTNESS_TF_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();status('failed');(O/'ROBUSTNESS_TF_FAILED').write_text(traceback.format_exc())
finally:(O/'robustness_tf.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)