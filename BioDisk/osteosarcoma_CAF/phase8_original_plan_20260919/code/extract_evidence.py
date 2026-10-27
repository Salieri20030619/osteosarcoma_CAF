from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,itertools,datetime
import numpy as np,pandas as pd,anndata as ad
import scipy.sparse as sp
from threadpoolctl import threadpool_limits
threadpool_limits(2)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase8_original_plan_20260919';D=P/'phase7_review_states_20260917';Q=P/'geo_preprocess_20260918'
panels=json.loads((O/'resources/programs_frozen.json').read_text())
cfg=json.loads((O/'resources/surface_candidates_frozen.json').read_text())
genes=sorted(set(sum(panels.values(),[])+cfg['surface_candidates']+cfg['controls']))
(O/'resources/extracted_gene_order.json').write_text(json.dumps(genes))
h=ad.read_h5ad(D/'integration_review/strict_CAF_harmony_review.h5ad',backed='r')
hlabels=h.obs['harmony_theta2_r0.5'].astype(str).copy();assert len(hlabels)==2047;h.file.close()
jobs=[]
for row in pd.read_csv(D/'review/input_summary.tsv',sep='\t').to_dict('records'):
 jobs.append({'sample':row['sample'],'source':'prior','path':P/'phase5_nonpdx_qc_20260916/results'/row['sample']/'annotated.h5ad','identity':D/'identity'/(row['sample']+'_identity.tsv.gz')})
for row in pd.read_csv(Q/'manifest.tsv',sep='\t').to_dict('records'):
 jobs.append({'sample':row['sample'],'source':'new_GEO','path':Q/'accepted'/row['sample']/'analysis_ready.h5ad','identity':P/'geo_identity_20260919/tables'/(row['sample']+'_identity.tsv.gz')})
rows=[];coverage=[];pb=[];pr=[];gate=[];integrity=[];summ=[]
candidates=cfg['surface_candidates'];gates=[(g,) for g in candidates]+list(itertools.combinations(candidates,2))
for job in jobs:
 sm=job['sample'];a=ad.read_h5ad(job['path'])
 o=pd.read_csv(job['identity'],sep='\t',dtype={'cluster':str}).set_index('cell_id',drop=False)
 assert o.index.tolist()==a.obs_names.tolist()
 symbols=a.var.symbol.astype(str).to_numpy();ix=[];jx=[]
 for j,g in enumerate(genes):
  v=np.flatnonzero(symbols==g);ix.extend(v);jx.extend([j]*len(v))
 mapping=sp.csr_matrix((np.ones(len(ix)),(ix,jx)),shape=(a.n_vars,len(genes)))
 x=a.layers['counts'].tocsr();raw=np.asarray((x@mapping).toarray(),dtype=np.float64);tot=np.asarray(x.sum(1)).ravel()
 assert np.isfinite(raw).all() and (raw>=0).all()
 log=np.log1p(raw/np.maximum(tot[:,None],1)*1e4)
 present=np.array([g in set(symbols) for g in genes])
 labelcol='lineage_reviewed' if job['source']=='prior' else 'lineage_phase3'
 line=o[labelcol].astype(str)
 if job['source']=='prior':
  strict=o.caf_strict.to_numpy();sens=np.zeros(len(o),bool)
  cls=o.final_identity.astype(str).copy()
  cls.loc[cls=='Other']=line[cls=='Other']
 else:
  strict=o.caf_supported_candidate.to_numpy();sens=o.count_threshold_sensitive_candidate.to_numpy()
  cls=line.copy()
  hi=o.joint.eq('both_high')&line.isin(['FibLike','OsteoLike','MesAmb','Mixed','Myogenic_pending'])
  cls.loc[hi]='CNV_high_mesenchymal'
  cls.loc[line.isin(['FibLike','MesAmb','Mixed','Unresolved'])&~hi]='Identity_pending'
 cls.loc[strict]='CAF_strict';cls.loc[sens]='CAF_detection_sensitive'
 hs=hlabels.reindex(o.index).fillna('not_discovery')
 group=cls.copy();group.loc[strict&(hs=='0')]='CAF_H0';group.loc[strict&(hs=='2')]='CAF_H2'
 patient=str(o.patient_id.iloc[0]);dataset=str(o.dataset.iloc[0]);site=str(o.site.iloc[0])
 meta=pd.DataFrame({'cell_id':o.index,'sample':sm,'dataset':dataset,'patient_id':patient,'site':site,'lineage':line.to_numpy(),'cell_class':cls.to_numpy(),'analysis_group':group.to_numpy(),'strict_CAF':strict,'sensitivity_CAF':sens,'discovery_state':hs.to_numpy(),'n_counts':tot,'n_genes':o.n_genes.to_numpy()})
 meta.to_csv(O/'tables'/(sm+'_cell_metadata.tsv.gz'),sep='\t',index=False)
 np.savez_compressed(O/'objects'/(sm+'_selected_gene_counts.npz'),raw=raw,totals=tot,present=present,genes=np.asarray(genes))
 for g,pres in zip(genes,present):coverage.append({'sample':sm,'dataset':dataset,'gene':g,'present':bool(pres)})
 for label in sorted(set(group)):
  mask=(group==label).to_numpy();n=int(mask.sum());lib=tot[mask].sum();counts=raw[mask].sum(0)
  for j,g in enumerate(genes):
   rows.append({'sample':sm,'dataset':dataset,'patient_id':patient,'site':site,'source':job['source'],'group':label,'cells':n,'gene':g,'present':bool(present[j]),'pct_detected':float((raw[mask,j]>0).mean()*100) if present[j] else np.nan,'mean_log1p_CP10K':float(log[mask,j].mean()) if present[j] else np.nan,'pseudobulk_count':float(counts[j]) if present[j] else np.nan,'pseudobulk_CPM':float(counts[j]/max(lib,1)*1e6) if present[j] else np.nan})
  for name,gs in panels.items():
   jj=[genes.index(g) for g in gs if g in genes and present[genes.index(g)]]
   if not jj:continue
   pr.append({'sample':sm,'dataset':dataset,'patient_id':patient,'site':site,'version':'raw','group':label,'cells':n,'program':name,'requested':len(gs),'covered':len(jj),'mean_log1p_CP10K':float(log[mask][:,jj].mean()),'pct_detect_ge2':float(((raw[mask][:,jj]>0).sum(1)>=2).mean()*100),'pct_detect_ge3':float(((raw[mask][:,jj]>0).sum(1)>=3).mean()*100)})
 # Candidate-only original/corrected records for external expression validation and robustness.
 ids=np.flatnonzero(strict)
 if len(ids):
  if job['source']=='new_GEO':
   b=ad.read_h5ad(P/'geo_identity_20260919/objects'/(sm+'_supported_CAF_candidates.h5ad'))
   corrected=b.layers['decontx_counts'].tocsr()
  else:
   cp=D/'remediation_20260917/decontx_comparison'/(sm+'_strict_corrected.h5ad')
   b=ad.read_h5ad(cp) if cp.exists() else None
   corrected=b.X.tocsr() if b is not None else None
  if b is not None:
   order=pd.Index(b.obs_names).get_indexer(o.index[ids]);assert (order>=0).all()
   bs=b.var.symbol.astype(str).to_numpy();ii=[];jj=[]
   for j,g in enumerate(genes):
    v=np.flatnonzero(bs==g);ii.extend(v);jj.extend([j]*len(v))
   mp=sp.csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(b.n_vars,len(genes)))
   cr=(corrected[order]@mp).toarray();ct=np.asarray(corrected[order].sum(1)).ravel()
   clog=np.log1p(cr/np.maximum(ct[:,None],1)*1e4)
   np.savez_compressed(O/'objects'/(sm+'_strict_gene_counts.npz'),cell_ids=o.index[ids].to_numpy(dtype=str),genes=np.asarray(genes),present=present,raw=raw[ids],corrected=cr,totals=tot[ids],corrected_totals=ct,discovery_states=hs.iloc[ids].to_numpy(dtype=str))
   for version,mat,lmat in [('raw',raw[ids],log[ids]),('corrected',cr,clog)]:
    for lab in ['all_strict']+(['H0','H2','other_CAF'] if job['source']=='prior' else []):
     mask=np.ones(len(ids),bool) if lab=='all_strict' else ((hs.iloc[ids]==lab[-1]).to_numpy() if lab in ['H0','H2'] else ~hs.iloc[ids].isin(['0','2']).to_numpy())
     if not mask.any():continue
     for name,gs in panels.items():
      js=[genes.index(g) for g in gs if g in genes and present[genes.index(g)]]
      if not js:continue
      pb.append({'sample':sm,'dataset':dataset,'patient_id':patient,'site':site,'group':lab,'version':version,'cells':int(mask.sum()),'program':name,'covered':len(js),'requested':len(gs),'mean_log1p_CP10K':float(lmat[mask][:,js].mean()),'pct_detect_ge2_gt0':float(((mat[mask][:,js]>0).sum(1)>=2).mean()*100),'pct_detect_ge2_count1':float(((mat[mask][:,js]>=1).sum(1)>=2).mean()*100)})
   del b,corrected,cr,clog
 # RNA co-detection screen only: no inferred protein gate performance.
 for threshold in [1,2]:
  dump=np.ones(len(o),bool)
  for g in ['PTPRC','PECAM1','EPCAM']:dump&=(raw[:,genes.index(g)]<threshold)
  for gs in gates:
   if not all(present[genes.index(g)] for g in gs):continue
   mask=dump.copy()
   for g in gs:mask&=(raw[:,genes.index(g)]>=threshold)
   n=int(mask.sum())
   for target in ['CAF_H0','CAF_H2','CAF_all']:
    truth=strict if target=='CAF_all' else (group==target).to_numpy()
    if not truth.any() and target!='CAF_all':continue
    mal=cls.isin(['Malignant_OS_candidate','Malignant_OS','CNV_high_mesenchymal']).to_numpy()
    unknown=cls.isin(['Identity_pending','CAF_detection_sensitive','Unresolved','MesAmb','Mixed']).to_numpy()
    gate.append({'sample':sm,'dataset':dataset,'patient_id':patient,'site':site,'source':job['source'],'target':target,'positive_genes':'+'.join(gs),'detection_count_threshold':threshold,'target_cells':int(truth.sum()),'gate_cells':n,'target_in_gate':int((mask&truth).sum()),'all_strict_CAF_in_gate':int((mask&strict).sum()),'CNV_high_mesenchymal_in_gate':int((mask&mal).sum()),'unresolved_in_gate':int((mask&unknown).sum()),'other_cells_in_gate':int((mask&~strict&~mal&~unknown).sum()),'target_coverage':float((mask&truth).sum()/truth.sum()) if truth.any() else np.nan,'known_CAF_fraction_gate':float((mask&strict).sum()/n) if n else np.nan,'CNV_high_fraction_gate':float((mask&mal).sum()/n) if n else np.nan})
 summ.append({'sample':sm,'source':job['source'],'dataset':dataset,'patient_id':patient,'site':site,'cells':len(o),'strict_CAF':int(strict.sum()),'sensitivity_CAF':int(sens.sum()),'H0':int((group=='CAF_H0').sum()),'H2':int((group=='CAF_H2').sum())})
 integrity.append({'sample':sm,'cells':len(o),'identity_order_match':True,'count_totals_match':bool(np.allclose(tot,o.n_counts))})
 assert np.allclose(tot,o.n_counts)
 print('EXTRACTED',sm,len(o),int(strict.sum()),flush=True)
 del a,x,raw,log,o,meta
pd.DataFrame(rows).to_csv(O/'tables/sample_group_gene_expression.tsv.gz',sep='\t',index=False)
pd.DataFrame(pr).to_csv(O/'tables/full_cell_programs.tsv',sep='\t',index=False)
pd.DataFrame(pb).to_csv(O/'tables/strict_programs_raw_corrected.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'tables/gene_coverage.tsv',sep='\t',index=False)
pd.DataFrame(gate).to_csv(O/'tables/RNA_marker_combination_screen.tsv.gz',sep='\t',index=False)
pd.DataFrame(summ).to_csv(O/'tables/analysis_sample_inventory.tsv',sep='\t',index=False)
pd.DataFrame(integrity).to_csv(O/'tables/extraction_integrity.tsv',sep='\t',index=False)
assert sum(r['H0'] for r in summ)==511 and sum(r['H2'] for r in summ)==448
(O/'EXTRACTION_COMPLETE').write_text('ORIGINAL_PLAN_DATA_EXTRACTION_COMPLETE\n')
print('EXTRACTION_COMPLETE',len(jobs),flush=True)
