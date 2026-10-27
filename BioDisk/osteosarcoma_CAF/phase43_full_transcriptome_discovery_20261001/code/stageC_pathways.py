# SYNTAX_CHECK: python -m py_compile stageC_pathways.py
# RUN_COMMAND: batch_bench/bin/python -u stageC_pathways.py
# EXPECTED_SUCCESS_MARKER: STAGE_C2_OPEN_PATHWAY_DISCOVERY_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,hashlib,collections,gc
import numpy as np,pandas as pd,gseapy as gp
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';O=P/'phase43_full_transcriptome_discovery_20261001'
for sub in ['tables/pathways','objects/gsea','resources','logs','delivery']:(O/sub).mkdir(parents=True,exist_ok=True)
def writej(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS_C2.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def bh(p):
 v=np.asarray(p,float);q=multipletests(np.where(np.isfinite(v),v,1),method='fdr_bh')[1];q[~np.isfinite(v)]=np.nan;return q
rc=1
try:
 assert (O/'full_gene_statistics.exit_code').read_text().strip()=='0' and (O/'C1_COMPLETE').read_text().strip()=='SUCCESS'
 maps=pd.read_csv(O/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False).set_index('feature');quality=maps.mapping_method.eq('HGNC_Ensembl_approved').map({True:0,False:1}).to_dict()
 hgnc=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);hgnc=hgnc[hgnc.status.eq('Approved')];approved=set(hgnc.symbol);aliases=collections.defaultdict(set)
 for row in hgnc.itertuples(index=False):
  for field in ['prev_symbol','alias_symbol']:
   for alias in getattr(row,field,'').split('|'):
    if alias:aliases[alias].add(row.symbol)
 def canon(s):return s if s in approved else (next(iter(aliases[s])) if s in aliases and len(aliases[s])==1 else s)
 libraries={}
 for lib,fname in [('hallmark','hallmark.gmt'),('reactome','reactome.gmt'),('go_bp','go_bp.gmt'),('matrisome','curated_cgp.gmt')]:
  sets={}
  for line in (A/'resources'/fname).read_text().splitlines():
   items=line.split('\t')
   if lib=='matrisome' and not items[0].startswith('NABA_'):continue
   sets[items[0]]=list(dict.fromkeys(canon(s) for s in items[2:]))
  libraries[lib]=sets
 writej(O/'resources/OPEN_PATHWAY_PROTOCOL_FROZEN.json',{'libraries':{lib:len(s) for lib,s in libraries.items()},'version':'MSigDB 2025.1.Hs, archived downloaded resources SHA256','ORA':'one-sided hypergeometric, raw gene P<0.05, up and down separate, full testable canonical-symbol background for each contrast','ORA_FDR_family':'layer|source|contrast|library|direction, all predefined library terms, unavailable terms retain NA','GSEA':'complete valid gene ranking by mean patient log2-pseudobulk effect; canonical symbol collapsed deterministically before observing P, full recorded library sets','GSEA_permutations':1000,'GSEA_seed':20261001,'GSEA_min_size':1,'GSEA_max_size':'background minus one, mathematical nonempty gene-set/complement','GSEA_P_type':'gene-set permutation/model P, not an independent biological treatment-effect P','GSEA_zero_P':'retain library-reported zero and finite-permutation warning; not an exact zero probability','GSEA_FDR':'retain native library q; separately report BH of model P in each library/contrast','small_N':'exploratory; patient-level program tests are separate primary biological statistics','gene_or_pathway_FDR_not_entry_gate':True,'no_upload_of_project_genes':True,'GSEA_dictionaries_from_local_downloads':True})
 contexts=json.loads((O/'resources/EFFECT_CONTEXTS.json').read_text());summary=[];coverage=[];failed=[];cached=None;gene=None
 for ci,ctx in enumerate(contexts,1):
  if ctx['N_matched_patients']<2:summary.append({**ctx,'state':'direction_only_no_enrichment_inference'});continue
  file=ctx['gene_result_file']
  if file!=cached:
   gene=pd.read_csv(O/file,sep='\t',keep_default_na=True);cached=file
  d=gene[gene.contrast.eq(ctx['contrast'])].copy();d=d[d.P.notna()&d.N.ge(2)&d.symbol.notna()&d.symbol.ne('')]
  d['symbol_priority']=d.feature.map(quality).fillna(1);d=d.sort_values(['symbol','symbol_priority','feature']).drop_duplicates('symbol',keep='first');background=set(d.symbol);M=len(background)
  if M<2:summary.append({**ctx,'state':'insufficient_testable_gene_background'});continue
  rank=d[['symbol','mean_delta']].sort_values(['mean_delta','symbol'],ascending=[False,True]);query={'up':set(d.loc[d.P.lt(.05)&d.median_delta.gt(0),'symbol']),'down':set(d.loc[d.P.lt(.05)&d.median_delta.lt(0),'symbol'])}
  for lib,sets in libraries.items():
   status('OPEN_PATHWAY_TEST',context=ci,total_contexts=len(contexts),layer=ctx['layer'],dataset=ctx['dataset'],contrast=ctx['contrast'],library=lib,background=M)
   orows=[];gsets={}
   for term,members in sets.items():
    gs=set(members)&background;K=len(gs)
    coverage.append({'context_id':ctx['context_id'],'library':lib,'term':term,'defined_genes':len(members),'tested_background_members':K,'background_genes':M,'missing_or_untestable_members':';'.join(s for s in members if s not in background)})
    if 0<K<M:gsets[term]=sorted(gs)
    for direction,q in query.items():
     n=len(q);over=q&gs;k=len(over);p=float(hypergeom.sf(k-1,M,K,n)) if n and K else np.nan
     orows.append({'layer':ctx['layer'],'dataset':ctx['dataset'],'contrast':ctx['contrast'],'kind':ctx['kind'],'library':lib,'direction':direction,'term':term,'background':M,'query_genes':n,'set_background_genes':K,'overlap_genes':k,'fold_enrichment':(k/n)/(K/M) if n and K else np.nan,'P':p,'overlap_symbols':';'.join(sorted(over)),'N_matched_patients':ctx['N_matched_patients'],'family':f"{ctx['layer']}|{ctx['dataset']}|{ctx['contrast']}|{lib}|{direction}"})
   ora=pd.DataFrame(orows);ora['FDR']=np.nan
   for direction,z in ora.groupby('direction'):ora.loc[z.index,'FDR']=bh(z.P)
   ora.to_csv(O/'tables/pathways'/f"{ctx['context_id']}__{lib}__ORA.tsv.gz",sep='\t',index=False)
   entry={**ctx,'library':lib,'background':M,'nominal_up_genes':len(query['up']),'nominal_down_genes':len(query['down']),'ORA_nominal_terms':int(ora.P.lt(.05).sum()),'ORA_FDR_terms':int(ora.FDR.lt(.05).sum()),'ranked_sets_evaluable':len(gsets)}
   if gsets and np.any(rank.mean_delta!=0):
    try:
     path=O/'objects/gsea'/f"{ctx['context_id']}__{lib}";path.mkdir(exist_ok=True)
     res=gp.prerank(rnk=rank,gene_sets=gsets,outdir=str(path),min_size=1,max_size=M-1,permutation_num=1000,threads=4,seed=20261001,no_plot=True,verbose=False)
     r=res.res2d.copy();r['layer']=ctx['layer'];r['dataset']=ctx['dataset'];r['contrast']=ctx['contrast'];r['kind']=ctx['kind'];r['library']=lib;r['N_matched_patients']=ctx['N_matched_patients'];r['P_type']='gene_set_permutation_model';r['GSEA_model_P']=pd.to_numeric(r['NOM p-val']);r['model_BH_FDR']=bh(r.GSEA_model_P);r['finite_permutation_zero_P']=r.GSEA_model_P.eq(0)
     r.to_csv(O/'tables/pathways'/f"{ctx['context_id']}__{lib}__GSEA.tsv.gz",sep='\t',index=False)
     entry.update(GSEA_rows=len(r),GSEA_nominal_terms=int(r.GSEA_model_P.lt(.05).sum()),GSEA_state='complete')
    except Exception as e:entry.update(GSEA_state='failed',error=repr(e));failed.append(entry.copy());traceback.print_exc()
   else:entry.update(GSEA_state='not_evaluable_no_rank_or_membership')
   summary.append(entry);pd.DataFrame(summary).to_csv(O/'tables/pathways/pathway_execution_summary.tsv',sep='\t',index=False);writej(O/'delivery/PATHWAY_PROGRESS.json',{'contexts_total':len(contexts),'library_outcomes':len(summary),'failed':failed,'at':datetime.datetime.now().astimezone().isoformat()});gc.collect()
 pd.DataFrame(coverage).to_csv(O/'tables/pathways/pathway_gene_coverage.tsv.gz',sep='\t',index=False)
 qa={'contexts':len(contexts),'library_outcomes':len(summary),'failed_GSEA':len(failed),'resource_download_proxy':False,'ORA_full_background':True,'all_ORA_terms_retained':True,'ranked_library_terms_math_evaluable_only':True,'patient_program_statistics_separate':True,'GSEA_model_P_not_biological_N_P':True,'libraries':{k:len(v) for k,v in libraries.items()}}
 writej(O/'delivery/OPEN_PATHWAY_QA.json',qa)
 if failed:raise RuntimeError(f'{len(failed)} GSEA outcomes failed; all failures retained')
 (O/'C2_COMPLETE').write_text('SUCCESS\n');(O/'COMPLETE').write_text('SUCCESS\n');status('complete',**qa);print('STAGE_C2_OPEN_PATHWAY_DISCOVERY_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();status('failed');(O/'C2_FAILED').write_text(traceback.format_exc())
finally:(O/'pathways.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)