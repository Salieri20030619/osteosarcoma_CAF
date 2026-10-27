"""SYNTAX_CHECK: python -m py_compile code/finalize_lineages_and_export.py
RUN_COMMAND: python code/finalize_lineages_and_export.py
EXPECTED_SUCCESS_MARKER: PHASE6_INPUTS_READY
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,gzip,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits
threadpool_limits(limits=4)
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p5=r/'phase5_nonpdx_qc_20260916';p6=r/'phase6_nonpdx_cnv_20260916'
# Explicit evidence review. Two named markers must remain in each cluster's top-50 list.
OVERRIDES=[
('GSE299023__001_BX_COM','2','OC_like','SIGLEC15,DCSTAMP','coherent osteoclast program'),
('GSE299023__001_BX_COM','9','DC','CD1C,FCER1A','dendritic lineage'),
('GSE299023__001_BX_COM','15','Endo','DLL4,ESM1','coherent arterial endothelial program'),
('GSE299023__001_BX_COM','22','DC','CLEC4C,LILRA4','plasmacytoid dendritic program'),
('GSE299023__003_BX_COM','3','OC_like','DCSTAMP,SIGLEC15','coherent osteoclast program'),
('GSE299023__003_BX_COM','14','OC_like','SIGLEC15,SLC9B2','osteoclast-associated program; CNV still evaluated'),
('GSE299023__003_BX_COM','16','T_NK','FOXP3,CTLA4','regulatory T-cell program'),
('GSE299023__010_BX_COM','7','MesAmb','MMP13,PTGDS','matrix/osteogenic overlap; CTSK alone does not prove osteoclast'),
('GSE299023__OS.C.21.BX.OCT','0','Myeloid','IL1B,CCL3','inflammatory myeloid program'),
('GSE299023__OS.C.21.BX.OCT','1','Myeloid','LYZ,S100A8','myeloid lineage'),
('GSE299023__OS.C.32.BX','4','Myeloid','LYZ,S100A9','myeloid lineage'),
('GSE152048__BC10','2','Chondrogenic_pending','COL2A1,CNMD','chondrogenic program; malignancy pending'),
('GSE152048__BC11','6','OsteoLike','PHOSPHO1,IFITM5','osteogenic program; malignancy pending'),
('GSE152048__BC11','11','FibLike','PI16,COL15A1','fibroblast-associated program; CNV required'),
('GSE152048__BC11','13','Epi','KRT6A,RAB25','epithelial program'),
('GSE152048__BC11','1','Chondrogenic_pending','COL2A1,MATN4','chondrogenic program; malignancy pending'),
('GSE152048__BC16','7','Myeloid','ADAM8,SLC11A1','myeloid program with lipid/inflammatory state'),
('GSE152048__BC17','14','Epi','AGR3,KRT19','epithelial program'),
('GSE152048__BC2','2','DC','FCER1A,CD1C','dendritic lineage'),
('GSE152048__BC2','7','DC','LILRA4,LRRC26','plasmacytoid dendritic program'),
('GSE152048__BC2','12','DC','CLEC9A,LAMP3','dendritic lineage'),
('GSE152048__BC2','13','Mixed','TPSAB1,AHSP','mast and erythroid overlap; not a reference'),
('GSE152048__BC20','14','Myeloid','CXCR2,FCGR3B','granulocytic myeloid program'),
('GSE152048__BC20','6','MesAmb','MFAP5,FAP','fibroblast/mural overlap; not resolved as mural'),
('GSE152048__BC21','15','Mixed','CD3E,SIGLEC15','T-cell and osteoclast overlap'),
('GSE152048__BC22','9','Chondrogenic_pending','COL2A1,COL9A1','chondrogenic program; malignancy pending'),
('GSE152048__BC22','2','MesAmb','NOTCH3,COL16A1','fibroblast/mural overlap'),
('GSE152048__BC3','10','Chondrogenic_pending','COL2A1,CNMD','chondrogenic program; malignancy pending'),
('GSE152048__BC5','1','FibLike','MFAP5,FGF7','fibroblast-associated program; CNV required'),
('GSE152048__BC5','13','Myogenic_pending','MYF5,PAX7','myogenic program; not CAF')]
manifest=pd.read_csv(p5/'manifest.tsv',sep='\t').fillna('');packets=pd.read_csv(p5/'review/cluster_review_packets.tsv',sep='\t',dtype={'cluster':str}).fillna('')
clinical=pd.read_csv(r/'phase4_download_20260914/clinical_source_audit.tsv',sep='\t',encoding='utf-8-sig').fillna('')
clinical.to_csv(p5/'review/clinical_context.tsv',sep='\t',index=False)
protocol={'seed':270231,'reference_selection_time':'before any new cohort CNV result','reference_lineages':['T_NK','B_Plasma','Myeloid','DC','Endo'],'min_genes':500,'min_UMI':750,'min_lineage_detected_genes':3,'min_own_lineage_UMI_percent':0.5,'max_combined_Fib_Osteo_Mural_UMI_percent':0.2,'max_competing_Epi_Myogenic_UMI_percent':0.2,'max_competing_immune_or_endo_UMI_percent':0.2,'max_OC_specific_UMI_percent':0.1,'exclude_doublet_score_above_sample_90th_percentile':True,'immune_cap':400,'endothelial_cap':160,'minimum_reference_lineage_eligible':20,'minimum_total_train':30,'minimum_total_holdout':10,'holdout_fraction':0.2,'holdout_min_per_eligible_lineage':10,'previous_detection_count_rule_retained':True,'ambient_RNA_claim':'relative competing expression screen only; no ambient correction claimed'}
(p6/'resources/reference_protocol.json').write_text(json.dumps(protocol,indent=2))
audits=[];summaries=[];labelcounts=[]
for row in manifest.to_dict('records'):
 sample=row['sample'];source=p5/'results'/sample;o=pd.read_csv(p5/'review'/(sample+'_cell_evidence.tsv.gz'),sep='\t',dtype={'cluster':str});a=ad.read_h5ad(source/'annotated.h5ad');x=a.layers['counts']
 assert o.cell_id.tolist()==a.obs_names.tolist()
 mk=pd.read_csv(source/'cluster_markers.tsv',sep='\t',dtype={'group':str})
 o['lineage_reviewed']=o.lineage_review_proposed.copy();o['expert_review_note']='marker_packet_reviewed; identity pending CNV'
 for sm,c,label,genes,reason in OVERRIDES:
  if sm!=sample:continue
  hits=set(genes.split(','))&set(mk.loc[mk.group==c,'symbol'])
  assert len(hits)>=2,(sm,c,genes)
  old=o.loc[o.cluster==c,'lineage_reviewed'].unique().tolist();o.loc[o.cluster==c,'lineage_reviewed']=label;o.loc[o.cluster==c,'expert_review_note']=reason
  audits.append({'sample':sm,'cluster':c,'old_label':';'.join(old),'reviewed_label':label,'supporting_top50_markers':';'.join(sorted(hits)),'reason':reason})
 pp=packets[packets['sample']==sample].set_index('cluster')
 cycle=o.cluster.map(pp.cycle_context).astype(bool)
 o.loc[(o.lineage_reviewed=='Unresolved')&cycle,'lineage_reviewed']='Cycling_pending'
 o['lineage_phase3']=o.lineage_reviewed
 for context in ['CTA_context','cycle_context']:o[context]=o.cluster.map(pp[context]).astype(bool)
 symbol=a.var.symbol.astype(str).to_numpy()
 def umi(genes):
  ix=np.flatnonzero(np.isin(symbol,genes));return np.asarray(x[:,ix].sum(1)).ravel()
 o['pctUMI_OC_specific']=umi(['OCSTAMP','DCSTAMP','ATP6V0D2','SIGLEC15','SUCNR1','CALCR'])/o.n_counts*100
 o['pctUMI_competing_mes']=o.pctUMI_Fib+o.pctUMI_Osteo+o.pctUMI_Mural
 immune_pct=sum(o['pctUMI_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC'])
 qscore=o.doublet_score.quantile(.9)
 base=(o.n_genes>=500)&(o.n_counts>=750)&o.doublet_score.notna()&(o.doublet_score<=qscore)&(o.pctUMI_competing_mes<=.2)&((o.pctUMI_Epi+o.pctUMI_Myogenic)<=.2)&(o.pctUMI_OC_specific<=.1)&~o.CTA_context
 o['ref_eligible_immune_original_count_rule']=o.ref_eligible_immune;o['ref_eligible_endo_original_count_rule']=o.ref_eligible_endo
 o['ref_eligible_immune']=False
 for k in ['T_NK','B_Plasma','Myeloid','DC']:
  o['ref_eligible_immune']|=(o.lineage_reviewed==k)&(o['det_'+k+'_review']>=3)&(o['pctUMI_'+k]>=.5)&base&(o.pctUMI_Endo<=.2)
 o['ref_eligible_endo']=(o.lineage_reviewed=='Endo')&(o.det_Endo_review>=3)&(o.pctUMI_Endo>=.5)&base&(immune_pct<=.2)
 immune_strong=(o.det_T_NK_review>=3)|(o.det_B_Plasma_review>=3)|(o.det_Myeloid_review>=4)
 o['mixed_lineage_fraction_flag']=(immune_strong&(immune_pct>.2)&(o.pctUMI_competing_mes>.2))|((o.det_Endo_review>=4)&(o.pctUMI_Endo>.2)&(o.pctUMI_Osteo>.2))
 target_labels=['FibLike','OsteoLike','Mural','MesAmb','Mixed','Unresolved','Myogenic_pending','Chondrogenic_pending','Cycling_pending','OC_like']
 o['cnv_target_expanded_original_count_rule']=o.cnv_target_expanded
 o['cnv_target_expanded']=o.lineage_reviewed.isin(target_labels)|(o.pctUMI_competing_mes>.2)|o.CTA_context
 o['cnv_role']='not_selected';o.loc[o.cnv_target_expanded,'cnv_role']='target'
 rng=np.random.default_rng(270231)
 for typ,cap in [('immune',400),('endo',160)]:
  ids=np.flatnonzero(o['ref_eligible_'+typ].to_numpy());rng.shuffle(ids);ids=ids[:cap]
  if len(ids)<20:continue
  nh=max(10,int(len(ids)*.2));o.loc[ids[:nh],'cnv_role']='holdout_'+typ;o.loc[ids[nh:],'cnv_role']='ref_'+typ
 refs=o.cnv_role.str.startswith('ref_');holds=o.cnv_role.str.startswith('holdout_')
 eligible=refs.sum()>=30 and holds.sum()>=10
 safe=sample.replace('.','_');o['tool_id']=[f'{safe}_c{i:06d}' for i in range(len(o))]
 o['patient_mapping_status']='known_within_dataset';o['sample']=sample;o['short']=row['short'];o['phase3_note']='lineage reviewed; CAF and malignancy remain pending CNV'
 ctxt=clinical[(clinical.dataset==row['dataset'])&(clinical['sample']==row['short'])].iloc[0]
 o['clinical_context_note']=ctxt['note'];o['clinical_unresolved']=ctxt['unresolved']
 o.to_csv(p6/'review'/(sample+'_cell_evidence.tsv.gz'),sep='\t',index=False)
 # Save review metadata separately; original phase5 expression objects remain immutable.
 labelcounts.extend({'sample':sample,'lineage':k,'cells':int(v)} for k,v in o.lineage_reviewed.value_counts().items())
 summary={'sample':sample,'short':row['short'],'patient_id':row['patient_id'],'site':row['site'],'clinical_group':row['group'],'retained':len(o),'fibroblast_lineage_candidates':int((o.lineage_reviewed=='FibLike').sum()),'target_expanded_total':int(o.cnv_target_expanded.sum()),'ref_immune':int((o.cnv_role=='ref_immune').sum()),'ref_endo':int((o.cnv_role=='ref_endo').sum()),'holdout_immune':int((o.cnv_role=='holdout_immune').sum()),'holdout_endo':int((o.cnv_role=='holdout_endo').sum()),'eligible_immune_count_rule':int(o.ref_eligible_immune_original_count_rule.sum()),'eligible_endo_count_rule':int(o.ref_eligible_endo_original_count_rule.sum()),'eligible_immune_fraction_rule':int(o.ref_eligible_immune.sum()),'eligible_endo_fraction_rule':int(o.ref_eligible_endo.sum()),'cnv_eligible':bool(eligible),'reason':'ready_for_dual_CNV' if eligible else 'insufficient_same_sample_train_or_holdout','CAF_confirmed':False}
 if eligible:
  sel=o.cnv_role!='not_selected';inp=p6/'inputs'/sample;inp.mkdir(exist_ok=True)
  o.loc[sel].to_csv(inp/'cells.tsv',sep='\t',index=False)
  genes=a.var[['symbol']].copy();genes.insert(0,'ensembl_id',a.var_names);genes.to_csv(inp/'genes.tsv',sep='\t',index=False)
  with gzip.open(inp/'counts.mtx.gz','wb') as f:mmwrite(f,x[sel.to_numpy()].T,field='integer')
  ann=pd.DataFrame({'tool_id':o.loc[sel,'tool_id'],'annotation':np.where(refs[sel],o.loc[sel,'cnv_role'],'C'+o.loc[sel,'cluster']+'_'+o.loc[sel,'lineage_reviewed'])})
  ann.to_csv(inp/'annotations.tsv',sep='\t',index=False,header=False)
  summary['run_cells']=int(sel.sum());summary['target_run']=int((o.cnv_role=='target').sum())
  (inp/'INPUT_COMPLETE').write_text('CNV_INPUT_EXPORT_COMPLETE\n')
 summaries.append(summary);pd.DataFrame(summaries).to_csv(p6/'review/input_summary.tsv',sep='\t',index=False)
 print('INPUT_READY',json.dumps(summary),flush=True)
 del a,x,o
pd.DataFrame(audits).to_csv(p6/'review/manual_lineage_review.tsv',sep='\t',index=False)
pd.DataFrame(labelcounts).to_csv(p6/'review/lineage_counts_reviewed.tsv',sep='\t',index=False)
(p6/'INPUTS_READY').write_text('PHASE6_INPUTS_READY\n')
(p6/'STATUS.json').write_text(json.dumps({'state':'INPUTS_READY','samples':len(summaries),'eligible_samples':sum(x['cnv_eligible'] for x in summaries),'cnv_started':False},indent=2))
print('PHASE6_INPUTS_READY',len(summaries),sum(x['cnv_eligible'] for x in summaries),flush=True)
