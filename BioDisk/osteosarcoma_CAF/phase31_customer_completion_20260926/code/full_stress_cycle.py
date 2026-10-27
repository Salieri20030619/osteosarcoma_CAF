from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,gc,hashlib
import pandas as pd,numpy as np,anndata as ad,scanpy as sc
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase31_customer_completion_20260926';B=P/'phase26_client_reanalysis_20260922';M=P/'phase30_two_layer_classification_20260926'
D=O/'stress_cycle';D.mkdir(exist_ok=True)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
rc=1
try:
 ev=read(M/'tables/expanded_classified_cells.tsv.gz').set_index('cell_id')
 st=pd.read_csv(O/'resources/stress_core_author.csv');cc=read(O/'resources/seurat_cc_updated2019.tsv')
 assert len(st)==512 and (st.logFC>0).sum()==507 and (st.logFC<0).sum()==5
 manifest=read(B/'tables/analysis_source_manifest.tsv');scores=[];coverage=[];qc=[]
 for i,row in enumerate(manifest.to_dict('records')):
  sm=row['sample'];a=ad.read_h5ad(B/'inputs'/sm/'full_called_counts.h5ad',backed='r')
  ob=a.obs.copy()
  for field in ['n_counts','n_genes','pct_mt','subsets_mito_percent']:
   if field in ob:
    x=pd.to_numeric(ob[field],errors='coerce');q=x.quantile([0,.25,.5,.75,1])
    qc.append(dict(sample=sm,field=field,cells=len(ob),missing=int(x.isna().sum()),minimum=q.iloc[0],q25=q.iloc[1],median=q.iloc[2],q75=q.iloc[3],maximum=q.iloc[4],new_filter_applied=False))
  ids=ev.index[ev['sample'].eq(sm)]
  if not len(ids):a.file.close();continue
  # Full measured gene space and all source cells define native score control bins.
  a=a.to_memory();x=a.X.astype(np.float32);den=np.asarray(x.sum(1)).ravel()
  a.X=x.multiply((1e4/np.maximum(den,1))[:,None]).tocsr();a.X.data=np.log1p(a.X.data)
  up=[z for z in st.loc[st.logFC>0,'ensembl_gene_id'] if z in a.var_names];down=[z for z in st.loc[st.logFC<0,'ensembl_gene_id'] if z in a.var_names]
  sy=a.var.symbol.astype(str)
  s=a.var_names[sy.isin(cc.loc[cc.phase=='S','gene'])].tolist();g=a.var_names[sy.isin(cc.loc[cc.phase=='G2M','gene'])].tolist()
  rowcov=dict(sample=sm,stress_up_requested=507,stress_down_requested=5,stress_up_measured=len(up),stress_down_measured=len(down),S_requested=int((cc.phase=='S').sum()),G2M_requested=int((cc.phase=='G2M').sum()),S_measured=len(s),G2M_measured=len(g),coverage_cutoff_applied=False)
  try:
   if not up or not down or not s or not g:raise ValueError('A required gene set has no measured genes; no replacement score invented')
   sc.tl.score_genes(a,up,score_name='stress507_up',random_state=20260926,use_raw=False)
   sc.tl.score_genes(a,down,score_name='stress5_down',random_state=20260926,use_raw=False)
   a.obs['stress512_signed']=a.obs.stress507_up-a.obs.stress5_down
   sc.tl.score_genes_cell_cycle(a,s_genes=s,g2m_genes=g,random_state=20260926,use_raw=False)
   cols=['stress507_up','stress5_down','stress512_signed','S_score','G2M_score','phase']
   out=ev.loc[ids,['sample','dataset','donor_record','class_id','in_supported','dual_reference']].join(a.obs.loc[ids,cols])
   assert len(out)==len(ids) and out.index.is_unique
   scores.append(out.reset_index(names='cell_id'));rowcov['status']='complete'
  except Exception as ex:rowcov.update(status='score_failed',error=str(ex))
  coverage.append(rowcov);writej(O/'SCORES_STATUS.json',{'status':'RUNNING','at':now(),'finished_source_index':i+1,'scored_samples':len(scores)})
  del a,x;gc.collect()
 out=pd.concat(scores,ignore_index=True);out.to_csv(D/'full512_and_cell_cycle_scores.tsv.gz',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(D/'score_coverage.tsv',sep='\t',index=False);pd.DataFrame(qc).to_csv(D/'all107_QC_distributions.tsv',sep='\t',index=False)
 for layer,sel in [('expanded',np.ones(len(out),bool)),('supported',out.in_supported.to_numpy())]:
  sub=out.loc[sel].copy()
  if layer=='supported':
   supported_labels=read(M/'tables/supported_classified_cells.tsv.gz').set_index('cell_id');sub['class_id']=sub.cell_id.map(supported_labels.class_id)
  sub.groupby(['class_id','dataset','donor_record'],observed=True)[['stress512_signed','S_score','G2M_score']].mean().to_csv(D/f'{layer}_patient_class_scores.tsv',sep='\t')
  pd.crosstab(sub.class_id,sub.phase).to_csv(D/f'{layer}_class_cell_cycle_composition.tsv',sep='\t')
 qa={'at':now(),'requested_expanded_cells':8174,'scored_expanded_cells':len(out),'scored_supported_cells':int(out.in_supported.sum()),'source_records_QC':len(manifest),'stress_requested':512,'Seurat_S':int((cc.phase=='S').sum()),'Seurat_G2M':int((cc.phase=='G2M').sum()),'software':'scanpy score_genes / score_genes_cell_cycle native defaults except deterministic seed and explicit count normalization','no_primary_cells_deleted':True,'score_failed_samples':[r for r in coverage if r['status']!='complete']}
 writej(D/'QA.json',qa);writej(O/'SCORES_STATUS.json',dict(status='SCORES_FINISHED_REVIEW_REQUIRED',**qa));(D/'COMPLETE').write_text(now()+chr(10));rc=0 if len(out)==8174 else 1
except Exception as ex:writej(O/'SCORES_STATUS.json',{'status':'FAILED','error':str(ex),'at':now()});traceback.print_exc()
finally:(D/'exit_code').write_text(str(rc)+chr(10))
raise SystemExit(rc)
