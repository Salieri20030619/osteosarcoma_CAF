from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,gzip,traceback,datetime,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase31_customer_completion_20260926';B=P/'phase26_client_reanalysis_20260922';M=P/'phase30_two_layer_classification_20260926';D=O/'ambient_comparison';D.mkdir(exist_ok=True)
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
def writej(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
def now():return datetime.datetime.now().astimezone().isoformat()
try:
 ev=read(M/'tables/expanded_classified_cells.tsv.gz').set_index('cell_id')
 su=read(M/'tables/supported_classified_cells.tsv.gz').set_index('cell_id')
 panels=json.loads((O/'resources/customer_program_panels.json').read_text())
 panels={k:v for k,v in panels.items() if not k.endswith('descriptive')}
 panels['Mural_review']=['RGS5','CSPG4','NOTCH3','DES','MCAM','PDGFRB']
 panels['Antigen_presentation_review']=['CD74','HLA-DRA','HLA-DRB1','CIITA','CD80','CD86']
 markers=sorted(set(sum(panels.values(),[])+['FAP','THY1','PDGFRA','PDGFRB','MCAM','JAG1','CD36','PTPRC','PECAM1','EPCAM']))
 results=[];cellrows=[];coverage=[];outcomes=[];doublets=[]
 for job in json.loads((O/'METHOD_JOBS.json').read_text()):
  sm=job['sample'];kind=job['method'];d=O/kind/sm
  ok=(d/'COMPLETE').exists() and (d/'exit_code').exists() and (d/'exit_code').read_text().strip()=='0'
  if not ok:outcomes.append(dict(sample=sm,method=kind,status='upstream_failed_or_unverified'));continue
  if kind=='doublet':
   r=read(d/'doublet_results.tsv').set_index('cell_id');ids=ev.index.intersection(r.index)
   r=r.loc[ids].join(ev[['sample','dataset','donor_record','in_supported','doublet_evidence']]);doublets.append(r.reset_index(names='cell_id'));outcomes.append(dict(sample=sm,method=kind,status='compared',candidate_cells=len(ids)));continue
  genes=read(B/'inputs'/sm/'genes.tsv');cells=read(B/'inputs'/sm/'cells.tsv');cont=read(d/'cell_contamination.tsv')
  assert cont.cell_id.tolist()==cells.cell_id.tolist()
  ids=ev.index[ev['sample'].eq(sm)]
  with gzip.open(d/'corrected_counts.mtx.gz','rb') as fh:
   with threadpool_limits(limits=4):y=mmread(fh).tocsr().T.tocsr()
  assert y.shape==(len(cells),len(genes)) and np.isfinite(y.data).all() and (y.data>=0).all()
  a=ad.read_h5ad(B/'inputs'/sm/'full_called_counts.h5ad',backed='r')
  assert list(a.obs_names)==cells.cell_id.tolist() and list(a.var_names)==genes.ensembl.tolist()
  # Full-sample ambient markers, without excluding any cells.
  original_den=cont.original_UMI.to_numpy(float);corrected_den=cont.corrected_UMI.to_numpy(float)
  assert np.allclose(np.asarray(y.sum(1)).ravel(),corrected_den,rtol=1e-5,atol=1e-4)
  measured={g:np.flatnonzero(genes.symbol.eq(g).to_numpy()) for g in markers}
  for gene,ix in measured.items():
   if not len(ix):coverage.append(dict(sample=sm,gene=gene,measured=False));continue
   coverage.append(dict(sample=sm,gene=gene,measured=True))
   raw=np.asarray(a[:,ix].X.sum(1)).ravel();corr=np.asarray(y[:,ix].sum(1)).ravel()
   results.append(dict(sample=sm,gene=gene,population='all_source_cells',cells=len(cells),original_total_counts=raw.sum(),corrected_total_counts=corr.sum(),original_positive_fraction=(raw>0).mean(),corrected_positive_fraction=(corr>0).mean(),original_mean_logCP10k=np.log1p(raw/np.maximum(original_den,1)*1e4).mean(),corrected_mean_logCP10k=np.log1p(corr/np.maximum(corrected_den,1)*1e4).mean()))
  if len(ids):
   index=pd.Index(cells.cell_id).get_indexer(ids);assert (index>=0).all()
   raw=a[index,:].X.tocsr();corr=y[index,:];rd=original_den[index];cd=corrected_den[index]
   out=ev.loc[ids,['sample','dataset','donor_record','class_id','in_supported']].copy().rename(columns={'class_id':'expanded_class_id'})
   out['supported_class_id']=su.class_id.reindex(ids);out['contamination']=cont.contamination.to_numpy()[index]
   for gene,ix in measured.items():
    if not len(ix):continue
    out['raw_'+gene]=np.asarray(raw[:,ix].sum(1)).ravel();out['corrected_'+gene]=np.asarray(corr[:,ix].sum(1)).ravel()
   for k,geneset in panels.items():
    ix=np.flatnonzero(genes.symbol.isin(geneset).to_numpy())
    if not len(ix):out['raw_program_'+k]=np.nan;out['corrected_program_'+k]=np.nan;continue
    r=raw[:,ix].multiply((1e4/np.maximum(rd,1))[:,None]).tocsr();r.data=np.log1p(r.data)
    c=corr[:,ix].multiply((1e4/np.maximum(cd,1))[:,None]).tocsr();c.data=np.log1p(c.data)
    out['raw_program_'+k]=np.asarray(r.mean(1)).ravel();out['corrected_program_'+k]=np.asarray(c.mean(1)).ravel()
   cellrows.append(out.reset_index(names='cell_id'))
  a.file.close();del a,y;gc.collect()
  outcomes.append(dict(sample=sm,method=kind,status='compared',candidate_cells=len(ids)))
  writej(O/'AMBIENT_STATUS.json',{'status':'RUNNING','outcomes':len(outcomes),'at':now()})
 pd.DataFrame(results).to_csv(D/'all_source_marker_before_after.tsv.gz',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(D/'measured_marker_coverage.tsv',sep='\t',index=False)
 pd.DataFrame(outcomes).to_csv(D/'method_output_integrity.tsv',sep='\t',index=False)
 out=pd.concat(cellrows,ignore_index=True) if cellrows else pd.DataFrame()
 if len(out):
  assert out.cell_id.is_unique
  out.to_csv(D/'candidate_cell_marker_program_before_after.tsv.gz',sep='\t',index=False)
  cols=[x for x in out if '_program_' in x]+['contamination']
  for layer,mask,cls in [('expanded',np.ones(len(out),bool),'expanded_class_id'),('supported',out.in_supported.to_numpy(bool),'supported_class_id')]:
   out.loc[mask].groupby(['dataset','donor_record',cls],dropna=False)[cols].mean().to_csv(D/f'{layer}_patient_class_before_after.tsv',sep='\t')
 if doublets:pd.concat(doublets,ignore_index=True).to_csv(D/'updated_doublet_evidence_candidates.tsv.gz',sep='\t',index=False)
 qa={'at':now(),'scheduled_methods':len(outcomes),'verified_decontx_samples':sum(r['method']=='decontx' and r['status']=='compared' for r in outcomes),'verified_doublet_samples':sum(r['method']=='doublet' and r['status']=='compared' for r in outcomes),'candidate_cells_compared':len(out),'expected_expanded_cells':8174,'all_failures_retained':[r for r in outcomes if r['status']!='compared'],'no_contamination_cell_exclusion_threshold':True,'all_primary_memberships_unchanged':True}
 writej(D/'QA.json',qa);(D/'COMPLETE').write_text(now());(D/'exit_code').write_text('0\n')
 writej(O/'AMBIENT_STATUS.json',dict(status='COMPARISON_COMPUTED_REVIEW_REQUIRED',**qa))
except Exception as ex:
 traceback.print_exc();(D/'exit_code').write_text('1\n');writej(O/'AMBIENT_STATUS.json',{'status':'FAILED','at':now(),'error':str(ex)});raise
