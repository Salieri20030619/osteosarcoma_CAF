from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='4'
import json, shutil, hashlib, datetime
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase16_scpca_supplement_20260920'
assert (O/'QC_COMPLETE').exists() and (O/'logs/preprocess.exit_code').read_text().strip()=='0'
for d in ['cnv_inputs','decontx','infercnv','copykat','review','delivery']: (O/d).mkdir(exist_ok=True)
old=P/'geo_cnv_20260918'
for name in ['infercnv_gene_order_GRCh38_98.tsv','coordinate_provenance.json']:shutil.copy2(old/'resources'/name,O/'resources'/name)
code=(old/'code/run_cnv.R').read_text().replace(str(old),str(O)).replace("file.path(root,'inputs',sample_id)","file.path(root,'cnv_inputs',sample_id)")
(O/'code/run_cnv.R').write_text(code)
code=(P/'phase10_batch_gates_20260919/code/run_decontx.R').read_text().replace(str(P/'phase10_batch_gates_20260919'),str(O)).replace("z<-as.integer(factor(cells$cluster))","z<-as.integer(factor(cells$ambient_group))")
code=code.replace("background='none; per-sample broad reviewed lineage z'","background='none; same-sample frozen broad lineage groups; sensitivity correction'")
(O/'code/run_decontx.R').write_text(code)
go=pd.read_csv(O/'resources/infercnv_gene_order_GRCh38_98.tsv',sep='\t',header=None,names=['id','chr','start','end'],dtype={'chr':str}).set_index('id')
manifest=pd.read_csv(O/'manifest.tsv',sep='\t',keep_default_na=False)
rows=[];audits=[]
for r in manifest.to_dict('records'):
 sm=r['scpca_sample_id'];a=ad.read_h5ad(O/'objects'/f'{sm}_qc.h5ad');o=a.obs.copy();o['cell_id']=a.obs_names
 assert not o.index.duplicated().any() and o.official_processed_member.all() and o.official_singlet.all() and o.analysis_keep.all()
 assert np.allclose(np.asarray(a.layers['counts'].sum(1)).ravel(),o.n_counts)
 # Preserve initial scope. Molecular evidence is additionally required per cell; no target-only CNV selection.
 o['lineage_phase3']=o.lineage_candidate.astype(str);o['phase2_label']=o.lineage_phase3;o['short']=sm
 o['ambient_group']=o.lineage_phase3
 if o.ambient_group.nunique()<2:o['ambient_group']='cluster_'+o.cluster.astype(str)
 ambient=bool(len(o)>=30 and o.ambient_group.nunique()>=2)
 o['CTA_context']=False
 for cl,idx in o.groupby('cluster',observed=True).groups.items():o.loc[idx,'CTA_context']=bool((o.loc[idx,'det_CTA']>=2).mean()>=.1)
 mes=o.pct_Fib+o.pct_Osteo+o.pct_Mural
 dbl=o.scDblFinder_score
 base=(o.n_genes>=500)&(o.n_counts>=750)&dbl.notna()&(dbl<=dbl.quantile(.9))&(mes<=.2)&((o.pct_Epi+o.pct_Myogenic)<=.2)&(o.pct_OC_specific<=.1)&~o.CTA_context
 immune=np.zeros(len(o),bool)
 for typ in ['T_NK','B_Plasma','Myeloid','DC']:
  immune|=((o.lineage_phase3==typ)&(o['det_'+typ]>=3)&(o['pct_'+typ]>=.5)&base&(o.pct_Endo<=.2)).to_numpy()
 ipct=o.pct_T_NK+o.pct_B_Plasma+o.pct_Myeloid+o.pct_DC
 endo=((o.lineage_phase3=='Endo')&(o.det_Endo>=3)&(o.pct_Endo>=.5)&base&(ipct<=.2)).to_numpy()
 o['cnv_role']='target';rng=np.random.default_rng(270231)
 for typ,mask,cap in [('immune',immune,400),('endo',endo,160)]:
  ids=np.flatnonzero(mask);rng.shuffle(ids);ids=ids[:cap]
  if len(ids)<20:continue
  n=max(10,int(len(ids)*.2));o.iloc[ids[:n],o.columns.get_loc('cnv_role')]='holdout_'+typ;o.iloc[ids[n:],o.columns.get_loc('cnv_role')]='ref_'+typ
 ref=o.cnv_role.str.startswith('ref_');hold=o.cnv_role.str.startswith('holdout_')
 if hold.sum()<20 and ref.sum()+hold.sum()>=50:
  ids=np.flatnonzero(ref);rng.shuffle(ids)
  for j in ids[:20-int(hold.sum())]:o.iloc[j,o.columns.get_loc('cnv_role')]=o.iloc[j].cnv_role.replace('ref_','holdout_')
 ref=o.cnv_role.str.startswith('ref_');hold=o.cnv_role.str.startswith('holdout_')
 assert not (ref&hold).any() and (mes[ref|hold]<=.2).all()
 assert o.loc[ref|hold,'lineage_phase3'].isin(['T_NK','B_Plasma','Myeloid','DC','Endo']).all()
 x=a.layers['counts'];expressed=np.asarray((x>0).sum(0)).ravel()>=3;ids=a.var_names
 mapped=go.reindex(ids);cov={str(ch):int((expressed&(mapped['chr'].fillna('').to_numpy()==str(ch))).sum()) for ch in range(1,23)}
 matched=int(mapped['chr'].notna().sum());audits.append({'sample':sm,'input_genes':a.n_vars,'GRCh38_98_id_matched':matched,'expressed_gene_coverage_by_autosome':json.dumps(cov),'minimum_autosomal_coverage':min(cov.values()),'mapping':'stable Ensembl ID; original ScPCA source metadata retained'})
 eligible=bool(ref.sum()>=30 and hold.sum()>=20 and min(cov.values())>=30)
 o['tool_id']=[sm+'_c'+str(i).zfill(6) for i in range(len(o))]
 o.to_csv(O/'review'/f'{sm}_cells.tsv.gz',sep='\t',index=False)
 # Add annotation columns without changing counts or cell order.
 o.to_csv(O/'inputs'/sm/'cells.tsv',sep='\t',index=False)
 if eligible:
  dest=O/'cnv_inputs'/sm;dest.mkdir(exist_ok=True)
  for name in ['counts.mtx.gz','genes.tsv','cells.tsv']:
   target=dest/name
   if not target.exists():target.symlink_to(O/'inputs'/sm/name)
  pd.DataFrame({'tool_id':o.tool_id,'annotation':np.where(ref,o.cnv_role,'C'+o.cluster.astype(str)+'_'+o.lineage_phase3)}).to_csv(dest/'annotations.tsv',sep='\t',index=False,header=False)
  (dest/'INPUT_COMPLETE').write_text('REFERENCE_HOLDOUT_AND_COUNTS_VALIDATED\n')
 rows.append({'sample':sm,'cells':len(o),'patient':r['participant_id'],'site':r['tissue_location'],'ambient_eligible':ambient,'cnv_eligible':eligible,'reference_cells':int(ref.sum()),'heldout_cells':int(hold.sum()),'min_autosomal_genes':min(cov.values()),'reason':'PASS' if eligible else 'insufficient_reference_holdout_or_gene_coverage','strict_candidate_scope':int(o.lineage_phase3.eq('FibLike').sum())})
 print('FOLLOWUP_INPUT',json.dumps(rows[-1]),flush=True)
pd.DataFrame(rows).to_csv(O/'review/input_summary.tsv',sep='\t',index=False)
pd.DataFrame(audits).to_csv(O/'resources/coordinate_audit.tsv',sep='\t',index=False)
protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'ambient':'per-sample DecontX, broad lineage grouping; no-background sensitivity; no single-group forced fit','CNV':'existing infercnv and CopyKAT code, stable Ensembl IDs, raw full-gene counts, independent heldout reference calibration','reference':'same-sample immune/endothelial only; minimum training30, heldout20; each autosome >=30 genes expressed in3 cells','CAF':'frozen FibLike scope plus exact phase10 raw/corrected molecular rules, dual low, reference performance and cluster risk gates; no unassessed-CNV rescue','negative_control_acceptance':'heldout evaluable>=20 for each method, >=90% low/diploid among evaluable; holdout call table retained','patient_statistics':'descriptive paired pseudobulk; no P values; multiple lung specimens pooled within patient and site; no cross-study independence assumed','completion':'every eligible job exits0 and complete marker/output integrity; unavailable data explicitly recorded; final report and file checksums'}
(O/'resources/followup_protocol.json').write_text(json.dumps(protocol,indent=2))
(O/'FOLLOWUP_INPUTS_READY').write_text('FOLLOWUP_INPUTS_READY\n')
