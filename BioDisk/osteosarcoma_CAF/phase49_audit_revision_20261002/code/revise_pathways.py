from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,hashlib,collections,shutil,traceback,concurrent.futures
import numpy as np,pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';R=P/'phase49_audit_revision_20261002'
def put(path,x):path.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
def bh(v):
 v=np.asarray(v,float);ok=np.isfinite(v);q=multipletests(np.where(ok,v,1),method='fdr_bh')[1];q[~ok]=np.nan;return q
def worker(job):
 import gseapy as gp
 cid,lib=job['context_id'],job['library'];rank=pd.read_csv(R/job['rank_file'],sep='	',float_precision='round_trip');sets=json.loads((R/'resources/libraries_canonical.json').read_text())[lib];background=set(rank.symbol);M=len(background);gsets={k:sorted(set(v)&background) for k,v in sets.items() if 0<len(set(v)&background)<M};dest=R/'tables/pathways'/f'{cid}__{lib}__GSEA.tsv.gz'
 if not gsets:return {**job,'state':'not_evaluable_no_sets'}
 d=R/'objects/gsea'/f'{cid}__{lib}';d.mkdir(exist_ok=True)
 res=gp.prerank(rnk=rank,gene_sets=gsets,outdir=str(d),min_size=1,max_size=M-1,permutation_num=1000,threads=4,seed=20261001,no_plot=True,verbose=False)
 out=res.res2d.copy()
 for k in ['layer','dataset','contrast','kind','N_matched_patients']:out[k]=job[k]
 out['library']=lib;out['P_type']='gene_set_permutation_model';out['GSEA_model_P']=pd.to_numeric(out['NOM p-val']);out['model_BH_FDR']=bh(out.GSEA_model_P);out['finite_permutation_zero_P']=out.GSEA_model_P.eq(0);out.to_csv(dest,sep='	',index=False);res.ranking.to_csv(d/'actual_ranking.tsv',sep='	',header=['rank_score'],index_label='symbol')
 return {**job,'state':'recomputed','terms':len(out),'GSEA_nominal':int(out.GSEA_model_P.lt(.05).sum())}
def main():
 cross=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='	',keep_default_na=False);priority=dict(zip(cross.feature,np.where(cross.mapping_method.eq('HGNC_Ensembl_approved'),0,1)))
 h=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='	',dtype=str,keep_default_na=False);h=h[h.status.eq('Approved')];approved=set(h.symbol);aliases=collections.defaultdict(set)
 for r in h.itertuples():
  for sy in (r.alias_symbol+'|'+r.prev_symbol).split('|'):
   if sy:aliases[sy].add(r.symbol)
 def canon(s):return s if s in approved else next(iter(aliases[s])) if s in aliases and len(aliases[s])==1 else s
 libraries={}
 for lib,fname in [('hallmark','hallmark.gmt'),('reactome','reactome.gmt'),('go_bp','go_bp.gmt'),('matrisome','curated_cgp.gmt')]:
  sets={}
  for line in (A/'resources'/fname).read_text().splitlines():
   v=line.split('	')
   if lib=='matrisome' and not v[0].startswith('NABA_'):continue
   sets[v[0]]=list(dict.fromkeys(canon(s) for s in v[2:]))
  libraries[lib]=sets
 put(R/'resources/libraries_canonical.json',libraries)
 protocol={'at':datetime.datetime.now().astimezone().isoformat(),'mapping':'per source full feature universe; HGNC Ensembl priority then lexicographic feature; frozen before P/effect/evaluability','background':'the frozen representative features with valid patient P and N>=2; no fallback to other feature after P NA','ORA_direction':'sign(mean_delta), target of patient signflip; median_delta retained as descriptive','families':'unchanged layer/source/contrast/library/direction','GSEA':'same1000gene_set_permutations_seed20261001_threads4; rerun changed ranking/background only','workers':2,'max_worker_threads':8,'original_outputs_preserved':True,'discovery_and_spatial_primary_statistics_unchanged':True}
 put(R/'resources/PATHWAY_REVISION_PROTOCOL_FROZEN.json',protocol)
 maps=[];summaries=[];jobs=[];registry=[];switches=[];coverage=[];ora_compare=[];ctxs=json.loads((C/'resources/EFFECT_CONTEXTS.json').read_text());cached=None
 for ctx in ctxs:
  cid=ctx['context_id'];file=ctx['gene_result_file']
  if file!=cached:
   full=pd.read_csv(C/file,sep='	',low_memory=False,float_precision='round_trip');features=full[['feature','symbol']].drop_duplicates();features=features[features.symbol.notna()&features.symbol.ne('')].copy();features['priority']=features.feature.map(priority).fillna(1);canonical=features.sort_values(['symbol','priority','feature']).drop_duplicates('symbol');canonical['layer']=ctx['layer'];canonical['dataset']=ctx['dataset'];maps.append(canonical);cached=file
  z=full[full.contrast.eq(ctx['contrast'])].copy();nom=z.P.lt(.05)
  summaries.append(dict(layer=ctx['layer'],dataset=ctx['dataset'],contrast=ctx['contrast'],kind=ctx['kind'],parent=ctx['parent'],gene_features=len(z),valid_tests=int(z.P.notna().sum()),N_max=int(z.N.max()),nominal_up=int((nom&z.mean_delta.gt(0)).sum()),nominal_down=int((nom&z.mean_delta.lt(0)).sum()),nominal_zero_mean=int((nom&z.mean_delta.eq(0)).sum()),median_zero_nominal=int((nom&z.median_delta.eq(0)).sum()),FDR_up=int((z.FDR.lt(.05)&z.mean_delta.gt(0)).sum()),FDR_down=int((z.FDR.lt(.05)&z.mean_delta.lt(0)).sum())))
  selected=z[z.feature.isin(canonical.feature)].copy();valid=selected[selected.P.notna()&selected.N.ge(2)];old=z[z.P.notna()&z.N.ge(2)&z.symbol.notna()&z.symbol.ne('')].copy();old['priority']=old.feature.map(priority).fillna(1);old=old.sort_values(['symbol','priority','feature']).drop_duplicates('symbol')
  keys=['symbol','mean_delta'];rank=valid[keys].sort_values(['mean_delta','symbol'],ascending=[False,True]);rankold=old[keys].sort_values(['mean_delta','symbol'],ascending=[False,True]);changed=not rank.reset_index(drop=True).equals(rankold.reset_index(drop=True))
  mapdict=canonical.set_index('symbol').feature.to_dict()
  for rr in old.itertuples():
   if mapdict[rr.symbol]!=rr.feature:switches.append(dict(context_id=cid,symbol=rr.symbol,old_feature=rr.feature,frozen_feature=mapdict[rr.symbol]))
  registry.append({**ctx,'mapping_rank_changed':changed,'old_background':len(old),'new_background':len(valid),'nominal_old_median_up':int((old.P.lt(.05)&old.median_delta.gt(0)).sum()),'nominal_mean_up':int((valid.P.lt(.05)&valid.mean_delta.gt(0)).sum()),'nominal_mean_down':int((valid.P.lt(.05)&valid.mean_delta.lt(0)).sum())})
  if ctx['N_matched_patients']<2 or len(valid)<2:continue
  rankfile=f'resources/{cid}__rank.tsv';rank.to_csv(R/rankfile,sep='	',index=False);background=set(valid.symbol);M=len(background);query={direction:set(valid.loc[valid.P.lt(.05)&(valid.mean_delta.gt(0) if direction=='up' else valid.mean_delta.lt(0)),'symbol']) for direction in ['up','down']}
  assert len(query['up']|query['down'])==int(valid.P.lt(.05).sum())
  for lib,sets in libraries.items():
   rows=[]
   for term,genes in sets.items():
    gs=set(genes)&background;K=len(gs);coverage.append(dict(context_id=cid,library=lib,term=term,defined_genes=len(genes),tested_background_members=K,background_genes=M,missing_or_untestable_members=';'.join(g for g in genes if g not in background)))
    for direction,q in query.items():
     k=len(q&gs);n=len(q);pv=float(hypergeom.sf(k-1,M,K,n)) if n and K else np.nan
     rows.append(dict(layer=ctx['layer'],dataset=ctx['dataset'],contrast=ctx['contrast'],kind=ctx['kind'],library=lib,direction=direction,term=term,background=M,query_genes=n,set_background_genes=K,overlap_genes=k,fold_enrichment=(k/n)/(K/M) if n and K else np.nan,P=pv,overlap_symbols=';'.join(sorted(q&gs)),N_matched_patients=ctx['N_matched_patients'],family=f"{ctx['layer']}|{ctx['dataset']}|{ctx['contrast']}|{lib}|{direction}"))
   ora=pd.DataFrame(rows);ora['FDR']=ora.groupby('family').P.transform(bh);ora.to_csv(R/'tables/pathways'/f'{cid}__{lib}__ORA.tsv.gz',sep='	',index=False)
   prev=pd.read_csv(C/'tables/pathways'/f'{cid}__{lib}__ORA.tsv.gz',sep='	');cmp=ora[['direction','term','P','FDR']].merge(prev[['direction','term','P','FDR']],on=['direction','term'],suffixes=('_new','_old'))
   ora_compare.append(dict(context_id=cid,library=lib,rows=len(cmp),nominal_old=int(cmp.P_old.lt(.05).sum()),nominal_new=int(cmp.P_new.lt(.05).sum()),nominal_status_changed=int((cmp.P_old.lt(.05)!=cmp.P_new.lt(.05)).sum()),FDR_old=int(cmp.FDR_old.lt(.05).sum()),FDR_new=int(cmp.FDR_new.lt(.05).sum()),FDR_status_changed=int((cmp.FDR_old.lt(.05)!=cmp.FDR_new.lt(.05)).sum())))
   if cid=='expanded__SCPCP000017__F4_vs_rest':cmp.to_csv(R/'tables'/f'F4_{lib}_ORA_before_after.tsv',sep='	',index=False)
   oldg=C/'tables/pathways'/f'{cid}__{lib}__GSEA.tsv.gz';job={**ctx,'library':lib,'rank_file':rankfile}
   if changed:jobs.append(job)
   elif oldg.exists():shutil.copy2(oldg,R/'tables/pathways'/oldg.name);jobs.append({**job,'reuse':True})
 pd.concat(maps,ignore_index=True).drop_duplicates().to_csv(R/'resources/source_canonical_feature_map.tsv',sep='	',index=False);pd.DataFrame(summaries).to_csv(R/'tables/full_gene_statistics_summary.tsv',sep='	',index=False);pd.DataFrame(registry).to_csv(R/'tables/pathway_input_revision_registry.tsv',sep='	',index=False);pd.DataFrame(switches).to_csv(R/'tables/representative_feature_changes.tsv',sep='	',index=False);pd.DataFrame(ora_compare).to_csv(R/'tables/ORA_before_after_summary.tsv',sep='	',index=False);pd.DataFrame(coverage).to_csv(R/'tables/pathways/pathway_gene_coverage.tsv.gz',sep='	',index=False)
 outcomes=[{**j,'state':'reused_identical_inputs'} for j in jobs if j.get('reuse')];run=[j for j in jobs if not j.get('reuse')];put(R/'resources/GSEA_REVISION_JOBS.json',run);put(R/'delivery/PATHWAY_PROGRESS.json',{'state':'GSEA_running','jobs':len(run),'reused':len(outcomes),'done':0});print('ORA_COMPLETE_GSEA_START',len(run),'REUSED',len(outcomes),flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
  future={pool.submit(worker,j):j for j in run}
  for f in concurrent.futures.as_completed(future):
   result=f.result();outcomes.append(result);put(R/'delivery/PATHWAY_PROGRESS.json',{'state':'GSEA_running','jobs':len(run),'done':sum(x['state']=='recomputed' for x in outcomes),'last':result['context_id']+'|'+result['library']});print('GSEA_DONE',result['context_id'],result['library'],flush=True)
 pd.DataFrame(outcomes).to_csv(R/'tables/pathway_execution_summary.tsv',sep='	',index=False);put(R/'delivery/PATHWAY_REVISION_QA.json',{'state':'complete','contexts':len(ctxs),'ORA_library_contexts':len(ora_compare),'GSEA_recomputed':len(run),'GSEA_reused':len(jobs)-len(run),'representative_changes':len(switches),'all_nominal_mean_direction_retained':True});(R/'PATHWAY_REVISION_COMPLETE').write_text('SUCCESS');print('PATHWAY_REVISION_COMPLETE',flush=True)
if __name__=='__main__':
 rc=1
 try:main();rc=0
 except Exception:traceback.print_exc();(R/'logs/pathway_revision_failure.txt').write_text(traceback.format_exc())
 (R/'pathway_revision.exit_code').write_text(str(rc));raise SystemExit(rc)
