from common import *
import os,select,subprocess,concurrent.futures,shutil
from scipy.stats import hypergeom
if not (O/'statistics.exit_code').exists():
 try:fd=os.pidfd_open(int((O/'statistics.pid').read_text()));select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
assert (O/'statistics.exit_code').read_text().strip()=='0' and (O/'STATISTICS_COMPLETE').is_file()
libs=json.loads((O/'resources/libraries_canonical.json').read_text());mp=read(O/'resources/source_canonical_feature_map.tsv');registry=read(O/'tables/executed_contrast_registry.tsv');registry=registry[registry.kind.eq('main')];jobs=[];summary=[];orasum=[];coverage=[]
(O/'resources/gsea_jobs').mkdir(exist_ok=True)
for ds,rs in registry.groupby('dataset',sort=False):
 full=read(O/'tables/genes'/f'expanded__{ds}__all_statistics.tsv.gz');mapping=mp[mp.dataset.eq(ds)]
 for row in rs.itertuples(index=False):
  ctx=f'expanded__{ds}__{row.contrast}'
  if row.N<2:summary.append({'context_id':ctx,'dataset':ds,'class_id':row.target,'N':row.N,'library':'all','state':'not_tested_N_lt2'});continue
  z=full[full.contrast.eq(row.contrast)];picked=mapping[['symbol','feature']].merge(z.drop(columns='symbol'),on='feature',how='left',validate='one_to_one');assert len(picked)==len(mapping);valid=picked[picked.P.notna()&picked.N.ge(2)].copy();rank=valid[['symbol','mean_delta']].sort_values(['mean_delta','symbol'],ascending=[False,True]);rankfile=O/'resources'/f'{ctx}__rank.tsv';rank.to_csv(rankfile,sep='\t',index=False);bg=set(rank.symbol);M=len(bg);qgenes={'up':set(valid.loc[valid.P.lt(.05)&valid.mean_delta.gt(0),'symbol']),'down':set(valid.loc[valid.P.lt(.05)&valid.mean_delta.lt(0),'symbol'])}
  oldrank=R/'resources'/f'{ctx}__rank.tsv';same=False
  if oldrank.exists():
   old=read(oldrank);same=old.symbol.tolist()==rank.symbol.tolist() and np.allclose(old.mean_delta,rank.mean_delta,rtol=0,atol=1e-12)
  for lib,sets in libs.items():
   oldgsea=R/'tables/pathways'/f'{ctx}__{lib}__GSEA.tsv.gz';out=O/'tables/gsea'/f'{ctx}__{lib}__GSEA.tsv.gz';base={'context_id':ctx,'dataset':ds,'class_id':row.target,'contrast':row.contrast,'N':row.N,'library':lib,'ranking_genes':M}
   if same and oldgsea.exists():shutil.copy2(oldgsea,out);summary.append({**base,'state':'verified_same_full_rank_reused_phase49','old_rank_max_abs_error':float(np.max(np.abs(old.mean_delta.to_numpy()-rank.mean_delta.to_numpy()))),'output':str(out),'SHA256':sha(out)})
   elif M>1 and np.any(rank.mean_delta!=0):
    job={**base,'job_id':ctx+'__'+lib,'rank':str(rankfile),'output':str(out)};jf=O/'resources/gsea_jobs'/f'{job["job_id"]}.json';dump(jf,job);jobs.append((job,jf))
   else:summary.append({**base,'state':'not_testable_empty_or_zero_rank'})
   orows=[]
   for term,defined in sets.items():
    members=set(defined)&bg;K=len(members);coverage.append({'context_id':ctx,'library':lib,'term':term,'definition_genes':len(defined),'tested_background_members':K,'background_genes':M})
    for direction,query in qgenes.items():
     nn=len(query);over=members&query;k=len(over);pv=float(hypergeom.sf(k-1,M,K,nn)) if M>1 and nn and 0<K<M else np.nan
     orows.append({'layer':'expanded','dataset':ds,'contrast':row.contrast,'class_id':row.target,'N_patients':row.N,'library':lib,'direction':direction,'term':term,'background':M,'query_genes':nn,'set_background_genes':K,'overlap_genes':k,'fold_enrichment':(k/nn)/(K/M) if M and nn and K else np.nan,'P':pv,'overlap_symbols':';'.join(sorted(over)),'family':ctx+'|'+lib+'|'+direction,'P_type':'overrepresentation_model'})
   ora=fdr(pd.DataFrame(orows));ora.to_csv(O/'tables/ora'/f'{ctx}__{lib}__ORA.tsv.gz',sep='\t',index=False);orasum.append({**base,'query_up':len(qgenes['up']),'query_down':len(qgenes['down']),'nominal_terms':int(ora.nominal.sum()),'FDR_terms':int(ora.FDR_significant.sum())})
pd.DataFrame(coverage).to_csv(O/'tables/pathway_coverage.tsv.gz',sep='\t',index=False);pd.DataFrame(orasum).to_csv(O/'tables/ORA_summary.tsv',sep='\t',index=False);dump(O/'resources/GSEA_EXECUTION_PLAN.json',{'new_jobs':len(jobs),'reused_results':sum(x['state'].startswith('verified') for x in summary),'max_concurrent_jobs':2,'threads_per_job':4,'ranking':'all valid tested canonical genes; no P threshold prefilter','permutations':1000,'seed':20261001})
def run(item):
 job,jf=item;log=O/'logs'/f'{job["job_id"]}__GSEA.log'
 with log.open('w') as f:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-u',str(O/'code/gsea_worker.py'),str(jf)],stdout=f,stderr=subprocess.STDOUT)
 (O/'logs'/f'{job["job_id"]}__GSEA.exit_code').write_text(str(r.returncode));return {**job,'exit_code':r.returncode,'state':'computed' if r.returncode==0 else 'failed','log':str(log)}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for f in concurrent.futures.as_completed([pool.submit(run,item) for item in jobs]):
  outcome=f.result();summary.append(outcome);pd.DataFrame(summary).to_csv(O/'tables/GSEA_execution_summary.tsv',sep='\t',index=False);print(json.dumps({'job':outcome['job_id'],'state':outcome['state']}),flush=True)
assert not any(x.get('exit_code',0)!=0 for x in summary)
pd.DataFrame(summary).to_csv(O/'tables/GSEA_execution_summary.tsv',sep='\t',index=False);dump(O/'delivery/GSEA_ACCEPTANCE.json',{'new_jobs':len(jobs),'reused':sum(x['state'].startswith('verified') for x in summary),'outcomes':len(summary),'gene_set_permutations':1000,'full_ranking_before_any_nominal_selection':True,'patient_statistics_separate':True});(O/'GSEA_COMPLETE').write_text('SUCCESS\n');print('GSEA_COMPLETE',flush=True)