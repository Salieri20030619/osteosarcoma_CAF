from pathlib import Path
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']: os.environ[key]='2'
os.environ['MPLBACKEND']='Agg'
import gzip,json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
from scipy import sparse
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
threadpool_limits(2)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'geo_identity_20260919';Q=P/'geo_preprocess_20260918';C=P/'geo_cnv_20260918'
assert (O/'AMBIENT_COMPLETE').exists()
assert (O/'logs/pipeline.exit_code').read_text().strip()=='0'
M=pd.read_csv(O/'tables/sample_provenance.tsv',sep='\t').fillna('')
panels=json.loads((Q/'review/lineage_rules.json').read_text())['panels']
panels['OC_specific']=['OCSTAMP','DCSTAMP','ATP6V0D2','SIGLEC15','SUCNR1','CALCR']
(O/'gene_panels.json').write_text(json.dumps(panels,indent=2))
cnv=pd.read_csv(C/'result_review_20260919/tables/cell_review.tsv.gz',sep='\t',index_col='cell_id',dtype={'cluster':str})
assert cnv.index.is_unique and len(cnv)==97238
final=[]; summaries=[]; references=[]; clusters=[]; markerrows=[]; coverage=[]; integrity=[]; ambient=[]; transitions=[]
genes_show=['DCN','LUM','PDGFRA','DPT','COL3A1','COL1A1','FAP','POSTN','RUNX2','SP7','ALPL','IBSP','RGS5','NOTCH3','MYH11','PECAM1','PTPRC','LYZ','EPCAM','SFTPC','MAGEA3']
def evidence(mat,symbols):
 totals=np.asarray(mat.sum(1)).ravel()
 result={}
 for name,genes in panels.items():
  n=np.zeros(mat.shape[0],dtype=int);npos=n.copy();sumc=np.zeros(mat.shape[0])
  for g in genes:
   v=np.asarray(mat[:,np.flatnonzero(symbols==g)].sum(1)).ravel()
   n+=(v>=1);npos+=(v>0);sumc+=v
  result['det_'+name]=n;result['pos_'+name]=npos
  result['pct_'+name]=sumc/np.maximum(totals,1)*100
 return pd.DataFrame(result),totals
def competing(e,det='det_',osteo=False):
 mask=np.zeros(len(e),dtype=bool)
 for name,n,pct in [('Mural',3,.2),('OC',2,.1),('T_NK',3,.2),('Myeloid',4,.2),('Endo',3,.2),('Epi',3,.2),('Myogenic',2,.2)]:
  pc='OC_specific' if name=='OC' else name
  mask|=(e[det+name]>=n)&(e['pct_'+pc]>pct)
 if osteo:mask|=(e[det+'Osteo']>=3)&(e.pct_Osteo>.2)
 return np.asarray(mask)
def mixed(e,det='det_'):
 # Fraction guard prevents trace transcripts being mistaken for a mixed cell.
 immune=sum(e['pct_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC'])
 mes=e.pct_Fib+e.pct_Osteo+e.pct_Mural
 return np.asarray(((immune>.2)&(mes>.2))|((e.pct_Endo>.2)&(mes>.2))|((e.pct_Epi>.2)&(mes>.2)))
for row in M.to_dict('records'):
 sm=row['sample'];a=ad.read_h5ad(Q/'accepted'/sm/'analysis_ready.h5ad')
 o=pd.read_csv(C/'review'/(sm+'_cells.tsv.gz'),sep='\t',dtype={'cluster':str}).set_index('cell_id',drop=False)
 assert a.obs_names.tolist()==o.index.tolist()
 o.index.name=None
 for col in ['infer_call','copy_call','joint','infer_score','copy_score','infer_cutoff','altered_chromosomes']:
  o[col]=cnv[col].reindex(o.index)
 for col in ['infer_call','copy_call','joint']:o[col]=o[col].fillna('not_assessed')
 for col in ['technology','source_type','patient_independence_status','library_mapping_status','clinical_pairing','cross_cohort_patient_pooling','treatment_raw','disease_stage_raw']:
  o[col]=row[col]
 o['dataset']=row['dataset'];o['site']=row['site'];o['patient_id']=row['patient_id']
 x=a.layers['counts'].tocsr();symbols=a.var.symbol.astype(str).to_numpy()
 raw,tot=evidence(x,symbols);raw.index=o.index
 assert np.allclose(tot,o.n_counts)
 for name in ['Fib','Osteo','Mural','CTA']:
  assert np.array_equal(raw['det_'+name],o['det_'+name+'_review']), (sm,name)
 for name,gs in panels.items():
  coverage.append({'sample':sm,'panel':name,'requested':len(gs),'covered':sum(g in set(symbols) for g in gs),'missing':';'.join(g for g in gs if g not in set(symbols))})
 for col in raw:o['raw_'+col]=raw[col]
 good={}
 for method,norm in [('infer','low_signal'),('copy','diploid')]:
  ok=sm in set(cnv['sample'])
  used=sorted(o.loc[o.cnv_role.str.startswith('ref_'),'cnv_role'].unique())
  ok=ok and bool(used)
  for role in used+[r.replace('ref_','holdout_') for r in used]:
   cc=o.loc[o.cnv_role==role,method+'_call'];ev=cc[~cc.isin(['not_evaluable','not_assessed'])]
   fraction=float((ev==norm).mean()) if len(ev) else np.nan
   passed=len(ev)>=10 and fraction>=.9
   references.append({'sample':sm,'method':method,'role':role,'total':len(cc),'evaluable':len(ev),'normal':int((ev==norm).sum()),'normal_fraction_evaluable':fraction,'missing_fraction':float(cc.isin(['not_evaluable','not_assessed']).mean()) if len(cc) else np.nan,'pass':passed})
   ok=ok and passed
  good[method]=bool(ok)
 o['reference_QC_pass']=good['infer'] and good['copy']
 o['cluster_both_high_fraction']=np.nan;o['cluster_any_high_fraction']=np.nan
 for cl,g in o.groupby('cluster'):
  binary=g.infer_call.isin(['low_signal','elevated_signal'])&g.copy_call.isin(['diploid','aneuploid'])
  ev=g.infer_call.isin(['low_signal','elevated_signal','ambiguous_signal'])&g.copy_call.isin(['diploid','aneuploid'])
  both=float(g.loc[ev,'joint'].eq('both_high').mean()) if ev.any() else np.nan
  anyhi=float((g.loc[binary,'infer_call'].eq('elevated_signal')|g.loc[binary,'copy_call'].eq('aneuploid')).mean()) if binary.any() else np.nan
  o.loc[g.index,'cluster_both_high_fraction']=both;o.loc[g.index,'cluster_any_high_fraction']=anyhi
  clusters.append({'sample':sm,'cluster':cl,'lineage':g.lineage_phase3.iloc[0],'cells':len(g),'joint_evaluable':int(ev.sum()),'binary_evaluable':int(binary.sum()),'both_high_fraction':both,'any_high_fraction_binary':anyhi})
 fib=o.lineage_phase3.eq('FibLike'); low=o.joint.eq('both_low')
 clean=(~mixed(raw))&(~o.CTA_context)&(raw.det_CTA<2)
 common=fib&(raw.det_Fib>=2)&low&clean&~competing(raw)&(o.cluster_both_high_fraction.fillna(1)<.2)&(o.n_genes>=500)&o.reference_QC_pass
 o['legacy_rule_candidate']=common
 o['raw_osteo_conflict']=(raw.det_Osteo>=3)&(raw.pct_Osteo>.2)
 o['raw_strict_candidate']=common&~o.raw_osteo_conflict&(o.cluster_any_high_fraction.fillna(1)<.2)
 o['ambient_checked']=False;o['decontx_contamination']=np.nan
 o['corrected_lineage_supported']=False;o['corrected_positive_sensitivity']=False
 corrected=None
 d=O/'decontx'/sm
 if (d/'COMPLETE').exists():
  assert (O/'logs'/(sm+'.exit_code')).read_text().strip()=='0'
  inputs=O/'inputs'/sm;ids=pd.read_csv(inputs/'cells.tsv',sep='\t');vg=pd.read_csv(inputs/'genes.tsv',sep='\t')
  assert ids.cell_id.tolist()==o.index.tolist() and vg.gene_id.tolist()==a.var_names.tolist()
  with gzip.open(d/'corrected_counts.mtx.gz','rb') as f:corrected=sparse.csr_matrix(mmread(f).T)
  assert corrected.shape==x.shape and np.isfinite(corrected.data).all() and (corrected.data>=0).all()
  diff=corrected-x
  assert not len(diff.data) or diff.data.max()<1e-5,sm
  ct=pd.read_csv(d/'cell_contamination.tsv',sep='\t'); assert ct.cell_id.tolist()==o.index.tolist()
  cor,ctot=evidence(corrected,symbols);cor.index=o.index
  assert np.allclose(ctot,ct.corrected_UMI) and np.allclose(tot,ct.original_UMI)
  o['ambient_checked']=True;o['decontx_contamination']=ct.contamination.to_numpy()
  for col in cor:o['corrected_'+col]=cor[col]
  o['corrected_lineage_supported']=(cor.det_Fib>=2)&~competing(cor,osteo=True)&~mixed(cor)&(cor.det_CTA<2)&(ctot>0)
  o['corrected_positive_sensitivity']=(cor.pos_Fib>=2)&~competing(cor,det='pos_',osteo=True)&~mixed(cor,det='pos_')&(cor.pos_CTA<2)&(ctot>0)
  ambient.append({'sample':sm,'cells':len(o),'fib_cells':int(fib.sum()),'median_contamination':float(ct.contamination.median()),'fib_median_contamination':float(o.loc[fib,'decontx_contamination'].median()),'raw_strict_median_contamination':float(o.loc[o.raw_strict_candidate,'decontx_contamination'].median()) if o.raw_strict_candidate.any() else np.nan})
 o['caf_supported_candidate']=o.raw_strict_candidate&o.corrected_lineage_supported
 o['caf_supported_positive_sensitivity']=o.raw_strict_candidate&o.corrected_positive_sensitivity
 o['technology_validation_open']=str(row['technology']).startswith('10x_FLEX')
 o['CAF_confirmed_independent_truth']=False
 o['identity_decision']='outside_existing_FibLike_scope'
 o.loc[fib,'identity_decision']='identity_unresolved'
 o.loc[fib&~low,'identity_decision']='CNV_unresolved_or_elevated'
 o.loc[fib&o.joint.eq('not_assessed'),'identity_decision']='CNV_unassessed'
 o.loc[fib&low,'identity_decision']='dual_low_additional_evidence_insufficient'
 o.loc[o.raw_strict_candidate,'identity_decision']='raw_supported_ambient_sensitive'
 o.loc[o.caf_supported_candidate,'identity_decision']='supported_CAF_candidate'
 o['exclusion_reasons']=''
 masks=[
  (~fib,'outside_FibLike_scope'),(~low,'not_dual_low'),(~o.reference_QC_pass,'reference_QC_not_passed'),
  (raw.det_Fib<2,'fewer_than_2_Fib_genes'),(o.n_genes<500,'fewer_than_500_genes'),
  (mixed(raw),'mixed_lineage_fraction'),(competing(raw),'competing_lineage_program'),
  (o.raw_osteo_conflict,'osteo_program'),(o.CTA_context|(raw.det_CTA>=2),'CTA_context_or_expression'),
  (o.cluster_both_high_fraction.fillna(1)>=.2,'cluster_concordant_high_risk'),
  (o.cluster_any_high_fraction.fillna(1)>=.2,'cluster_any_method_high_risk'),
  (~o.ambient_checked,'ambient_not_assessed'),
  (o.ambient_checked&~o.corrected_lineage_supported,'corrected_lineage_not_supported')]
 for mask,reason in masks:o.loc[mask,'exclusion_reasons']+=reason+';'
 # Keep each excluded cell, even when multiple exclusion gates apply.
 assert (~o.caf_supported_candidate|o.raw_strict_candidate).all()
 assert not (o.caf_supported_candidate&~low).any()
 assert o.loc[o.caf_supported_candidate,'exclusion_reasons'].eq('').all()
 o['decision_scope']='operational_tumor_associated_fibroblast_candidate; no new patient pooling'
 o.to_csv(O/'tables'/(sm+'_identity.tsv.gz'),sep='\t',index=False)
 final.append(o.reset_index(drop=True))
 summary={'sample':sm,'dataset':row['dataset'],'site':row['site'],'technology':row['technology'],'cells':len(o),'FibLike':int(fib.sum()),'dual_low_FibLike':int((fib&low).sum()),'legacy_rule_candidate':int(o.legacy_rule_candidate.sum()),'raw_strict_candidate':int(o.raw_strict_candidate.sum()),'supported_candidate':int(o.caf_supported_candidate.sum()),'positive_count_sensitivity':int(o.caf_supported_positive_sensitivity.sum()),'reference_QC_pass':good['infer'] and good['copy'],'ambient_checked':corrected is not None}
 summaries.append(summary)
 for decision,g in o[fib].groupby('identity_decision'):
  transitions.append({'sample':sm,'decision':decision,'cells':len(g)})
 for version,mat in [('raw',x),('decontx',corrected)]:
  if mat is None:continue
  total=np.asarray(mat.sum(1)).ravel()
  for group,mask in [('all_FibLike',fib.to_numpy()),('raw_strict',o.raw_strict_candidate.to_numpy()),('supported',o.caf_supported_candidate.to_numpy()),('dual_low_rejected',(fib&low&~o.caf_supported_candidate).to_numpy())]:
   if not mask.any():continue
   for gene in genes_show:
    val=np.asarray(mat[:,np.flatnonzero(symbols==gene)].sum(1)).ravel()[mask]
    log=np.log1p(val/np.maximum(total[mask],1)*1e4)
    markerrows.append({'sample':sm,'version':version,'group':group,'gene':gene,'cells':int(mask.sum()),'count_ge1_pct':float((val>=1).mean()*100),'count_gt0_pct':float((val>0).mean()*100),'mean_count':float(val.mean()),'mean_log1p':float(log.mean())})
 if o.caf_supported_candidate.any():
  mask=o.caf_supported_candidate.to_numpy();b=a[mask].copy()
  b.obs=o.loc[mask].drop(columns='cell_id').copy()
  for col in b.obs.select_dtypes(include='object'):b.obs[col]=b.obs[col].fillna('').astype(str)
  b.layers['decontx_counts']=corrected[mask].copy()
  b.uns['identity_status']='operational CAF candidates; independent biological truth not established'
  b.uns['identity_rules']='geo_identity_20260919/decision_rules.json'
  b.uns['annotation_status']='joint_lineage_CNV_and_ambient_review_complete; operational_candidates'
  b.uns['CAF_status']='SUPPORTED_OPERATIONAL_CANDIDATE_NOT_INDEPENDENTLY_CONFIRMED'
  b.uns['patient_pooling']='not permitted by this identity review; cross-cohort overlap unresolved'
  b.write_h5ad(O/'objects'/(sm+'_supported_CAF_candidates.h5ad'),compression='gzip')
  check=ad.read_h5ad(O/'objects'/(sm+'_supported_CAF_candidates.h5ad'))
  assert check.obs_names.equals(o.index[mask])
  assert (check.layers['counts']!=x[mask]).nnz==0
  assert (check.layers['decontx_counts']!=corrected[mask]).nnz==0
  del b,check
 integrity.append({'sample':sm,'accepted_cells':len(o),'cell_order_match':True,'raw_counts_conserved':True,'raw_panel_counts_reproduced':True,'corrected_order_values_verified':corrected is not None,'supported_object_verified':bool(o.caf_supported_candidate.any())})
 print('IDENTITY',json.dumps(summary),flush=True)
 del a,x,corrected,o
allcells=pd.concat(final,ignore_index=True)
assert len(allcells)==101487 and allcells.cell_id.is_unique
assert allcells.lineage_phase3.eq('FibLike').sum()==8870
assert (allcells.lineage_phase3.eq('FibLike')&allcells.joint.eq('both_low')).sum()==3028
allcells.to_csv(O/'tables/all_cell_identity_evidence.tsv.gz',sep='\t',index=False)
allcells.loc[allcells.caf_supported_candidate,['cell_id','sample','patient_id','dataset','site','technology','identity_decision','technology_validation_open','patient_independence_status','library_mapping_status']].to_csv(O/'tables/supported_CAF_candidates.tsv',sep='\t',index=False)
for name,data in [('sample_identity_summary',summaries),('reference_QC',references),('cluster_CNV_risk',clusters),('marker_before_after',markerrows),('gene_panel_coverage',coverage),('integrity_checks',integrity),('ambient_summary',ambient),('identity_transitions',transitions)]:
 pd.DataFrame(data).to_csv(O/'tables'/(name+'.tsv'),sep='\t',index=False)
s=pd.DataFrame(summaries);s=s[s.FibLike>0]
fig,ax=plt.subplots(figsize=(12,6));yy=np.arange(len(s));left=np.zeros(len(s))
for col,lab,color in [('supported_candidate','Supported candidate','#287D65'),('raw_strict_candidate','Ambient-sensitive raw candidate','#DFA83D'),('dual_low_FibLike','Other dual-low','#718DAC'),('FibLike','CNV elevated / unresolved','#C8CDD2')]:
 if col=='supported_candidate':v=s[col].to_numpy()
 elif col=='raw_strict_candidate':v=(s[col]-s.supported_candidate).to_numpy()
 elif col=='dual_low_FibLike':v=(s[col]-s.raw_strict_candidate).to_numpy()
 else:v=(s[col]-s.dual_low_FibLike).to_numpy()
 assert (v>=0).all()
 ax.barh(yy,v,left=left,label=lab,color=color);left+=v
ax.set_yticks(yy,[v.split('__')[1]+(' [lung]' if site=='lung' else '') for v,site in zip(s['sample'],s.site)])
ax.invert_yaxis();ax.set_xlabel('Cells (descriptive counts)');ax.set_title('CAF identity review | original FibLike candidates')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2,frameon=False);fig.tight_layout();fig.savefig(O/'figures/identity_funnel.png',dpi=160,bbox_inches='tight');plt.close(fig)
mk=pd.DataFrame(markerrows);ms=mk[(mk.group=='supported')&(mk.version=='decontx')]
if len(ms):
 order=[v for v in M['sample'] if v in set(ms['sample'])]
 fig,ax=plt.subplots(figsize=(13,max(4,len(order)*.45+2)))
 for i,sm in enumerate(order):
  mm=ms[ms['sample']==sm].set_index('gene').reindex(genes_show)
  pts=ax.scatter(np.arange(len(genes_show)),np.repeat(i,len(genes_show)),s=mm.count_ge1_pct*.9,c=mm.mean_log1p,cmap='Reds',vmin=0,vmax=3)
 ax.set_xticks(range(len(genes_show)),genes_show,rotation=90);ax.set_yticks(range(len(order)),[v.split('__')[1] for v in order]);ax.invert_yaxis()
 ax.set_title('Supported candidates | DecontX expression\nDot area: fraction with corrected count >=1; color: mean log1p(CP10k)')
 fig.colorbar(pts,ax=ax,label='Mean log1p(CP10k)');fig.tight_layout();fig.savefig(O/'figures/supported_markers.png',dpi=160,bbox_inches='tight');plt.close(fig)
state={'status':'IDENTITY_TABLES_COMPLETE_VISUAL_QA_PENDING','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':16,'cells':len(allcells),'FibLike':8870,'dual_low_FibLike':3028,'legacy_rule_candidates':int(allcells.legacy_rule_candidate.sum()),'raw_strict_candidates':int(allcells.raw_strict_candidate.sum()),'supported_CAF_candidates':int(allcells.caf_supported_candidate.sum()),'positive_count_sensitivity':int(allcells.caf_supported_positive_sensitivity.sum()),'confirmed_CAF_independent_truth':None,'new_pooled_independent_patient_N':None,'inference_tests_performed':False,'original_labels_changed':False}
(O/'STATUS.json').write_text(json.dumps(state,indent=2));(O/'TABLES_COMPLETE').write_text('IDENTITY_TABLES_COMPLETE_VISUAL_QA_PENDING\n')
print(json.dumps(state),flush=True)
