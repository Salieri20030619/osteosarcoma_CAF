from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad,json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase22_no_UMI_sensitivity_20260921'
for d in ['tables','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'change':'remove project QC n_counts/UMI lower gate only; freeze original per-sample gene and MT thresholds','all52_samples_audited':True,'retained':'official processed membership and singlet calls for ScPCA; existing biological identity gates not changed','legacy_extra_doublets':'previous UMI-excluded cells did not receive original downstream doublet calls; report as pending, not newly accepted singlets','scope':'one-factor QC and downstream eligibility sensitivity, with measured added-cell gene programs for ScPCA; no relabeling of unassessed CNV or automatic assignment to C7/C9','UMI_reference_selection':'normal-reference quality rules are retained; cancelling QC UMI does not cancel independent normal-reference criteria','baseline':'frozen299483 QC cells and Phase17/20 labels/statistics preserved'}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
legacy=pd.read_csv(O/'legacy_UMI_audit.tsv',sep='\t');base=pd.read_csv(P/'phase17_all_cohort_expanded_20260921/tables/sample_layer_counts.tsv',sep='\t').set_index('sample')
rows=[];addedmeta=[];programs=[]
for r in legacy.to_dict('records'):
 sm=r['sample'];q=pd.read_csv(r['QC_table'],sep='\t',low_memory=False);thr=json.loads((Path(r['QC_table']).parent/'qc_thresholds.json').read_text())
 mask=(q.n_genes>=thr['min_genes'])&(q.pct_mt<=thr['max_pct_mt'])&(q.n_counts<thr['min_counts'])
 a=q[mask].copy();a['dataset']=base.loc[sm,'dataset'];a['supplement_status']='QC_rescue_pending_doublet_and_CNV';addedmeta.append(a[['cell_id','sample','dataset','n_genes','n_counts','pct_mt','supplement_status']])
 rows.append({'sample':sm,'dataset':base.loc[sm,'dataset'],'baseline_QC':int(base.loc[sm,'QC_cells']),'UMI_only_additions':int(mask.sum()),'added_known_singlets':0,'added_doublet_unassessed':int(mask.sum()),'no_UMI_QC_ceiling':int(base.loc[sm,'QC_cells']+mask.sum()),'additional_genes_ge300':int((a.n_genes>=300).sum()),'additional_raw_molecular_gate':np.nan,'normal_reference_compatible_additions':np.nan,'interpretation':'QC upper bound; additions not accepted as final singlets'})
root=P/'phase16_scpca_supplement_20260920';qsum=pd.read_csv(root/'tables/qc_summary.tsv',sep='\t');manifest=pd.read_csv(root/'manifest.tsv',sep='\t').set_index('scpca_sample_id');panels=json.loads((P/'phase10_batch_gates_20260919/resources/panels.json').read_text());customer=json.loads((root/'resources/customer_panels.json').read_text())
for r in qsum.to_dict('records'):
 sm=r['sample'];q=pd.read_csv(root/'tables'/f'{sm}_all_cell_qc.tsv.gz',sep='\t',low_memory=False)
 keep=q.official_processed_member&q.official_singlet&(q.n_genes>=r['min_genes'])&(q.pct_mt<=r['max_mt']);add=keep&~q.analysis_keep;a=q[add].copy()
 rr={'sample':sm,'dataset':r['project'],'baseline_QC':int(q.analysis_keep.sum()),'UMI_only_additions':int(add.sum()),'added_known_singlets':int(add.sum()),'added_doublet_unassessed':0,'no_UMI_QC_ceiling':int(keep.sum()),'additional_genes_ge300':int((a.n_genes>=300).sum()),'additional_raw_molecular_gate':0,'normal_reference_compatible_additions':0,'interpretation':'ScPCA singlets pass remaining frozen QC; new CNV not assessed'}
 if len(a):
  src=ad.read_h5ad(manifest.loc[sm,'filtered'],backed='r');ids=a.barcode.astype(str).tolist();assert set(ids)<=set(src.obs_names)
  x=src[ids,:].X.tocsr();sy=src.var.gene_symbol.astype(str).to_numpy();src.file.close();den=np.asarray(x.sum(1)).ravel();assert np.allclose(den,a.n_counts)
  z=pd.DataFrame(index=a.cell_id.to_numpy());xc=x.tocsc()
  for name,gs in {**panels,**customer}.items():
   vs=[np.asarray(xc[:,sy==g].sum(1)).ravel() for g in gs if (sy==g).any()];v=np.column_stack(vs) if vs else np.zeros((len(a),0))
   z['det_'+name]=(v>=1).sum(1);z['pct_'+name]=v.sum(1)/np.maximum(den,1)*100
   if name in customer:programs.append({'sample':sm,'program':name,'added_cells':len(a),'raw_ge2_cells':int(((v>=1).sum(1)>=2).sum()),'coverage':v.shape[1]/len(gs),'identity_status':'all added singlets, not confirmed CAF'})
  competing=np.zeros(len(a),bool)
  for name,ng in [('T_NK',3),('B_Plasma',3),('Myeloid',3),('DC',3),('Endo',3),('Epi',3),('Myogenic',3),('OC_specific',2)]:competing|=((z['det_'+name]>=ng)&(z['pct_'+name]>1)).to_numpy()
  molecular=(a.n_genes.to_numpy()>=300)&((z.det_Fib>=2)|((z.det_Fib>=1)&(z.det_ECM>=3))).to_numpy()&~competing
  rr['additional_raw_molecular_gate']=int(molecular.sum());z['raw_molecular_gate']=molecular
  mes=z.pct_Fib+z.pct_Osteo+z.pct_Mural;immune=z.pct_T_NK+z.pct_B_Plasma+z.pct_Myeloid+z.pct_DC
  normalbase=(a.n_genes.to_numpy()>=500)&(den>=750)&(mes<=.2)&((z.pct_Epi+z.pct_Myogenic)<=.2)&(z.pct_OC_specific<=.1)&(z.det_CTA<2)
  normal=np.zeros(len(a),bool)
  for name in ['T_NK','B_Plasma','Myeloid','DC']:normal|=((z['det_'+name]>=3)&(z['pct_'+name]>=.5)&normalbase&(z.pct_Endo<=.2)).to_numpy()
  normal|=((z.det_Endo>=3)&(z.pct_Endo>=.5)&normalbase&(immune<=.2)).to_numpy()
  z['normal_reference_molecular_compatible']=normal;rr['normal_reference_compatible_additions']=int(normal.sum())
  z['sample']=sm;z.to_csv(O/'tables'/f'{sm}_added_expression_evidence.tsv.gz',sep='\t',index_label='cell_id')
  a['dataset']=r['project'];a['supplement_status']='official_singlet_QC_rescued_CNV_unassessed';addedmeta.append(a[['cell_id','sample','dataset','n_genes','n_counts','pct_mt','supplement_status']])
 rows.append(rr);print('SAMPLE',rr,flush=True)
S=pd.DataFrame(rows);S.to_csv(O/'tables/no_UMI_sample_sensitivity.tsv',sep='\t',index=False);pd.concat(addedmeta,ignore_index=True).to_csv(O/'tables/UMI_only_added_records.tsv.gz',sep='\t',index=False);pd.DataFrame(programs).to_csv(O/'tables/ScPCA_added_program_support.tsv',sep='\t',index=False)
assert len(S)==52 and S.baseline_QC.sum()==299483 and S.UMI_only_additions.sum()==6566
q={'samples':52,'baseline_QC':299483,'UMI_only_additions':int(S.UMI_only_additions.sum()),'additional_official_singlets':int(S.added_known_singlets.sum()),'legacy_added_doublet_unassessed':int(S.added_doublet_unassessed.sum()),'QC_ceiling_if_all_legacy_additions_pass':int(S.no_UMI_QC_ceiling.sum()),'baseline_plus_known_ScPCA_singlets':int(299483+S.added_known_singlets.sum()),'additional_genes_ge300_all':int(S.additional_genes_ge300.sum()),'ScPCA_additional_raw_molecular_candidates':int(S.additional_raw_molecular_gate.sum()),'ScPCA_normal_reference_compatible_additions':int(S.normal_reference_compatible_additions.sum()),'new_confirmed_CAF':None,'new_CNV_jobs_run':0,'original_results_overwritten':False}
(O/'delivery/QA.json').write_text(json.dumps(q,indent=2));(O/'ANALYSIS_COMPLETE').write_text('NO_UMI_QC_AND_ELIGIBILITY_SENSITIVITY_COMPLETE\n');print(json.dumps(q),flush=True)
