from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,pandas as pd,numpy as np,sys,gseapy as gp,subprocess,concurrent.futures
from statsmodels.stats.multitest import multipletests
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase61_F3_internal_portraits_20261005');jobs=json.loads((O/'resources/GSEA_JOBS.json').read_text())
if len(sys.argv)>1:
 job=jobs[int(sys.argv[1])];libs=json.loads((O/'resources/libraries_canonical.json').read_text());rank=pd.read_csv(job['rank'],sep='\t');bg=set(rank.symbol);M=len(bg);sets={k:sorted(set(v)&bg) for k,v in libs[job['library']].items() if 0<len(set(v)&bg)<M};res=gp.prerank(rnk=rank[['symbol','mean_delta']],gene_sets=sets,outdir=None,min_size=1,max_size=M-1,permutation_num=1000,threads=4,seed=20261005,no_plot=True,verbose=False);x=res.res2d.copy();x['GSEA_model_P']=pd.to_numeric(x['NOM p-val']);x['model_BH_FDR']=multipletests(x.GSEA_model_P.fillna(1),method='fdr_bh')[1]
 for k in ['context','matrix','contrast','scope','library','target_cells','reference_cells','target_donors','reference_donors']:x[k]=job[k]
 x['P_type']='gene_set_permutation_model_not_biological_patient_test';x['finite_permutation_zero_P']=x.GSEA_model_P.eq(0);x.to_csv(job['output'],sep='\t',index=False);print('JOB_COMPLETE',job['context'],job['library'],len(x),flush=True)
else:
 assert (O/'PREPARE_COMPLETE').exists()
 def run(i):
  log=O/'logs'/f'gsea_{i:02d}.log'
  with log.open('w') as f:r=subprocess.run([sys.executable,str(Path(__file__)),str(i)],stdout=f,stderr=subprocess.STDOUT)
  (O/'logs'/f'gsea_{i:02d}.exit_code').write_text(str(r.returncode));return i,r.returncode
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for i,code in pool.map(run,range(len(jobs))):
   print('GSEA',i,code,flush=True)
   if code:raise RuntimeError(f'GSEA job {i} failed')
 frames=[pd.read_csv(j['output'],sep='\t') for j in jobs];pd.concat(frames,ignore_index=True).to_csv(O/'tables/all_GSEA_results.tsv.gz',sep='\t',index=False);(O/'GSEA_COMPLETE').write_text('SUCCESS\n');print('GSEA_COMPLETE',flush=True)
