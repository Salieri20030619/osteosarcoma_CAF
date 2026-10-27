from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,gzip,shutil,subprocess,datetime,fcntl,ast
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.io import mmwrite,mmread
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');OLD=P/'phase16_scpca_supplement_20260920';SUP=P/'phase22_no_UMI_sensitivity_20260921';O=SUP/'CNV611';SM='SCPCS000611'
for d in ['code','resources','inputs','cnv_inputs','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
lock=(O/'runner.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(s,**kw):(O/'STATUS.json').write_text(json.dumps({'status':s,'updated':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw},indent=2));print(s,kw,flush=True)
state('PREPARING')
meta=pd.read_csv(OLD/'manifest.tsv',sep='\t').set_index('scpca_sample_id').loc[SM]
q=pd.read_csv(OLD/'tables'/f'{SM}_all_cell_qc.tsv.gz',sep='\t').set_index('cell_id');thr=pd.read_csv(OLD/'tables/qc_summary.tsv',sep='\t').set_index('sample').loc[SM]
keep=q.official_processed_member&q.official_singlet&(q.n_genes>=thr.min_genes)&(q.pct_mt<=thr.max_mt)
q=q[keep].copy();assert len(q)==3787 and (~q.analysis_keep).sum()==9
src=ad.read_h5ad(meta.filtered,backed='r');x=src[q.barcode.tolist()].X.tocsr();genes=pd.DataFrame({'gene_id':src.var_names,'ensembl_id':src.var_names,'symbol':src.var.gene_symbol.to_numpy()});src.file.close()
assert np.allclose(np.asarray(x.sum(1)).ravel(),q.n_counts)
ob=pd.read_csv(OLD/'review'/f'{SM}_cells.tsv.gz',sep='\t').set_index('cell_id').reindex(q.index)
new=~q.analysis_keep;added=pd.read_csv(SUP/'tables'/f'{SM}_added_expression_evidence.tsv.gz',sep='\t').set_index('cell_id')
for col in ['sample','patient_id','site','n_counts','n_genes','barcode','scDblFinder_score']:
 ob[col]=q[col]
ob.loc[new,'cluster']='extra_no_UMI';ob.loc[new,'lineage_phase3']='Unresolved';ob.loc[new,'phase2_label']='Unresolved';ob.loc[new,'CTA_context']=False;ob.loc[new,'cnv_role']='target'
for c in added:
 if c.startswith(('det_','pct_')):ob.loc[added.index,c]=added[c]
normal_ids=added.index[added.normal_reference_molecular_compatible]
assert len(normal_ids)==2 and (ob.loc[normal_ids,'scDblFinder_score']<=q.loc[q.analysis_keep,'scDblFinder_score'].quantile(.9)).all()
assert (ob.loc[normal_ids,'det_CTA']<2).all()
ob.loc[normal_ids,'lineage_phase3']='Endo';ob.loc[normal_ids,'phase2_label']='Endo';ob.loc[normal_ids,'cnv_role']='ref_endo'
refs=ob.index[ob.cnv_role.str.startswith('ref_')].to_numpy();holds=ob.index[ob.cnv_role.str.startswith('holdout_')]
assert len(refs)==41 and len(holds)==10
rng=np.random.default_rng(270231);move=rng.choice(refs,10,replace=False)
for ix in move:ob.loc[ix,'cnv_role']=ob.loc[ix,'cnv_role'].replace('ref_','holdout_')
assert ob.cnv_role.str.startswith('ref_').sum()==31 and ob.cnv_role.str.startswith('holdout_').sum()==20
ob['short']=SM;ob['ambient_group']=ob.lineage_phase3;ob['tool_id']=[SM+'_noUMI_'+str(i).zfill(6) for i in range(len(ob))];ob['sample']=SM;ob['cell_id']=ob.index
for name in ['infercnv_gene_order_GRCh38_98.tsv','coordinate_audit.tsv','coordinate_provenance.json']:shutil.copy2(P/'geo_cnv_20260918/resources'/name,O/'resources'/name)
go=pd.read_csv(O/'resources/infercnv_gene_order_GRCh38_98.tsv',sep='\t',header=None,names=['id','chr','start','end'],dtype={'chr':str}).set_index('id');expressed=np.asarray((x>0).sum(0)).ravel()>=3;ch=go.reindex(genes.gene_id)['chr'].fillna('').to_numpy();cov={str(i):int((expressed&(ch==str(i))).sum()) for i in range(1,23)};assert min(cov.values())>=30
inp=O/'inputs'/SM;inp.mkdir(exist_ok=True);ob.to_csv(inp/'cells.tsv',sep='\t',index=False);genes.to_csv(inp/'genes.tsv',sep='\t',index=False)
with gzip.open(inp/'counts.mtx.gz','wb',compresslevel=1) as f:mmwrite(f,x.T,field='integer')
(inp/'COMPLETE').write_text('COUNTS_VERIFIED\n');cnv=O/'cnv_inputs'/SM
if not cnv.exists():cnv.symlink_to(inp,target_is_directory=True)
ref=ob.cnv_role.str.startswith('ref_');pd.DataFrame({'tool_id':ob.tool_id,'annotation':np.where(ref,ob.cnv_role,'C'+ob.cluster.astype(str)+'_'+ob.lineage_phase3)}).to_csv(inp/'annotations.tsv',sep='\t',index=False,header=False)
(inp/'INPUT_COMPLETE').write_text('CNV_REFERENCE_INPUT_READY\n')
for name in ['run_cnv.R','run_cnv_job.sh','run_decontx.R']:
 code=(OLD/'code'/name).read_text().replace(str(OLD),str(O));(O/'code'/name).write_text(code)
review={'QC_cells':len(ob),'added_cells':9,'added_normal_candidates':normal_ids.tolist(),'train_reference':31,'heldout_calibration':20,'minimum_autosome_gene_coverage':min(cov.values()),'normal_reference_UMI750_rule_retained':True,'new_unresolved_cells':7,'original_results_overwritten':False}
(O/'INPUT_REVIEW.json').write_text(json.dumps(review,indent=2));state('CNV_AND_AMBIENT_RUNNING',**review)
def run(method):
 cmd=['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/run_decontx.R'),SM] if method=='decontx' else ['bash',str(O/'code/run_cnv_job.sh'),method,SM]
 with (O/'logs'/f'{method}.log').open('w') as f:rc=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT).returncode
 (O/'logs'/f'{method}.exit_code').write_text(str(rc)+'\n');assert rc==0,(method,rc)
 assert (O/method/SM/'COMPLETE').exists();print('METHOD_COMPLETE',method,flush=True)
 return method
with ThreadPoolExecutor(max_workers=3) as pool:
 for fut in [pool.submit(run,m) for m in ['infercnv','copykat','decontx']]:fut.result()
# Reuse frozen molecular gate functions, evaluate this sample only.
tree=ast.parse((OLD/'code/finalize_followup.py').read_text());functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['evidence','conflict']]
panels=json.loads((P/'phase10_batch_gates_20260919/resources/panels.json').read_text());ns={'np':np,'pd':pd,'panels':panels};exec(compile(ast.Module(body=functions,type_ignores=[]),'<frozen_gate_functions>','exec'),ns)
sy=genes.symbol.astype(str).to_numpy();raw,lib=ns['evidence'](x,sy)
with gzip.open(O/'decontx'/SM/'corrected_counts.mtx.gz','rb') as f:cor=sparse.csr_matrix(mmread(f).T)
corr,clib=ns['evidence'](cor,sy);raw.index=ob.index;corr.index=ob.index
calls={};cal=[]
for method,low in [('infercnv','low_signal'),('copykat','diploid')]:
 v=pd.read_csv(O/method/SM/'cell_results.tsv',sep='\t').set_index('cell_id').loc[ob.index];assert v.tool_id.tolist()==ob.tool_id.tolist();calls[method]=v.call
 hold=ob.cnv_role.str.startswith('holdout_');valid=v.call.isin(['low_signal','elevated_signal','ambiguous_signal'] if method=='infercnv' else ['diploid','aneuploid']);n=int((hold&valid).sum());fraction=float(v.loc[hold&valid,'call'].eq(low).mean());cal.append({'method':method,'holdout_evaluable':n,'low_fraction':fraction,'pass':n>=20 and fraction>=.9})
pd.DataFrame(cal).to_csv(O/'delivery/reference_acceptance.tsv',sep='\t',index=False)
ic=calls['infercnv'];cc=calls['copykat'];dual=ic.eq('low_signal')&cc.eq('diploid');high=ic.eq('elevated_signal')|cc.eq('aneuploid');low=ic.eq('low_signal')|cc.eq('diploid');refpass=all(r['pass'] for r in cal)
frac=pd.Series(np.nan,index=ob.index)
for cluster,g in ob.groupby('cluster'):
 valid=ic.loc[g.index].isin(['low_signal','elevated_signal'])&cc.loc[g.index].isin(['diploid','aneuploid'])
 if valid.any():frac.loc[g.index]=high.loc[g.index][valid].mean()
strict=(ob.lineage_phase3=='FibLike')&dual&refpass&~ob.CTA_context.astype(bool)&(raw.det_CTA<2)&(corr.det_CTA<2)&(ob.n_genes>=500)&(raw.det_Fib>=2)&(corr.det_Fib>=2)&(frac.fillna(1)<.2)&(clib>0)&~ns['conflict'](raw)&~ns['conflict'](corr)
old_ev=pd.read_csv(P/'phase17_all_cohort_expanded_20260921/tables/all_cell_layer_evidence.tsv.gz',sep='\t');broad=old_ev.set_index('cell_id').expanded.reindex(ob.index).fillna(False).astype(bool)
broad.loc[added.index]=added.raw_molecular_gate
broad=broad&~(high&~low);supported=broad&low&~high&refpass
customer=json.loads((OLD/'resources/customer_panels.json').read_text())['myCAF_customer'];rx=x.tocsc();cx=cor.tocsc();rv=np.column_stack([np.asarray(rx[:,sy==g].sum(1)).ravel() for g in customer]);cv=np.column_stack([np.asarray(cx[:,sy==g].sum(1)).ravel() for g in customer]);my=(rv>=1).sum(1)>=2;mycor=my&((cv>=1).sum(1)>=2)
out=pd.DataFrame({'cell_id':ob.index,'baseline_QC':q.analysis_keep.to_numpy(),'infer_call':ic.to_numpy(),'copykat_call':cc.to_numpy(),'reference_QC_pass':refpass,'strict':strict.to_numpy(),'supported':supported.to_numpy(),'broad':broad.to_numpy(),'raw_myCAF':my,'raw_and_corrected_myCAF':mycor});out.to_csv(O/'delivery/identity_evidence.tsv.gz',sep='\t',index=False)
result={**review,'all_methods_complete':True,'reference_QC_pass':refpass,'strict_candidates':int(strict.sum()),'supported_candidates':int(supported.sum()),'broad_candidates':int(broad.sum()),'strict_raw_myCAF':int((strict&my).sum()),'supported_raw_myCAF':int((supported&my).sum()),'strict_raw_corrected_myCAF':int((strict&mycor).sum()),'both_CNV_low_cells':int(dual.sum()),'CNV_high_cells':int(high.sum()),'interpretation':'single-sample no-UMI QC sensitivity; original global clusters and primary statistics not refit'}
(O/'delivery/QA.json').write_text(json.dumps(result,indent=2));(O/'ANALYSIS_COMPLETE').write_text('NO_UMI_SAMPLE_CNV_RECOVERY_COMPLETE\n');state('COMPLETE',**result)
