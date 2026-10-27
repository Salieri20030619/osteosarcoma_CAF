"""SYNTAX_CHECK: python -m py_compile code/compare_lineage_programs.py
RUN_COMMAND: python -u code/compare_lineage_programs.py
EXPECTED_SUCCESS_MARKER: LINEAGE_PROGRAM_COMPARISON_COMPLETE
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='4'
from pathlib import Path
import json,gc
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad
from threadpoolctl import threadpool_limits
threadpool_limits(limits=4)
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=r/'phase7_review_states_20260917';out=p/'analysis'
panels=json.loads((out/'gene_panels.json').read_text());panels['Pericyte_support']=['RGS5','CSPG4','KCNJ8','ABCC9','COX4I2','PDGFRB'];panels['Vascular_SMC_support']=['MYH11','ACTG2','CNN1','LMOD1','MYLK','TAGLN'];genes=sorted(set(sum(panels.values(),[])+['FAP','PDPN','PDGFRA','PDGFRB','THY1','CD34','ITGA11','CD36']))
manifest=pd.read_csv(p/'review/input_summary.tsv',sep='\t');program_rows=[];gene_rows=[];groups=[];mural_rows=[]
for row in manifest.to_dict('records'):
 sample=row['sample'];a=ad.read_h5ad(r/'phase5_nonpdx_qc_20260916/results'/sample/'annotated.h5ad');o=pd.read_csv(p/'identity'/(sample+'_identity.tsv.gz'),sep='\t');assert o.cell_id.tolist()==a.obs_names.tolist()
 cl=o.final_identity.astype(str).copy();cl.loc[o.final_identity=='Other']=o.loc[o.final_identity=='Other','lineage_reviewed'].astype(str);cl.loc[o.final_identity=='Fib_candidate']='Fib_candidate_non_strict';cl.loc[o.caf_strict]='CAF_strict'
 symbol=a.var.symbol.astype(str).to_numpy();mapping_rows=[];mapping_cols=[]
 for j,g in enumerate(genes):
  hits=np.flatnonzero(symbol==g);mapping_rows.extend(hits.tolist());mapping_cols.extend([j]*len(hits))
 mapping=sp.csr_matrix((np.ones(len(mapping_rows)),(mapping_rows,mapping_cols)),shape=(a.n_vars,len(genes)))
 raw=(a.layers['counts']@mapping).toarray();expr=(a.X@mapping).toarray();total=np.asarray(a.layers['counts'].sum(1)).ravel()
 per=(raw[:,[genes.index(g) for g in panels['Pericyte_support']]]>0).sum(1)
 smc=(raw[:,[genes.index(g) for g in panels['Vascular_SMC_support']]]>0).sum(1)
 subtype=np.where((per>=3)&(smc<3),'Pericyte_like',np.where((smc>=3)&(per<3),'Vascular_SMC_like',np.where((per>=3)&(smc>=3),'Mural_mixed','Mural_unresolved')))
 verified=(o.final_identity=='Mural_SMC').to_numpy();cl.loc[verified]=subtype[verified]
 candidate=(o.lineage_reviewed=='Mural').to_numpy()|verified
 for i in np.flatnonzero(candidate):mural_rows.append({'sample':sample,'cell_id':o.cell_id.iloc[i],'final_identity':o.final_identity.iloc[i],'mural_subtype':subtype[i],'pericyte_markers_detected':int(per[i]),'SMC_markers_detected':int(smc[i]),'status':'RNA_marker_supported; original CNV status retained'})
 # Redraw the full-cell atlas using reviewed labels and existing per-sample UMAP coordinates.
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axs=plt.subplots(1,2,figsize=(14,5))
 for ax,values,title in [(axs[0],cl,'Reviewed lineage and identity'),(axs[1],o.final_identity,'CNV-supported identity')]:
  labels=sorted(values.unique());colors=plt.get_cmap('tab20')
  for j,label in enumerate(labels):
   ix=(values==label).to_numpy();ax.scatter(a.obsm['X_umap'][ix,0],a.obsm['X_umap'][ix,1],s=2,color=colors(j%20),label=label,rasterized=True)
  ax.set_title(sample+' '+title,fontsize=9);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(fontsize=5.5,markerscale=3,loc='upper left',bbox_to_anchor=(1.01,1))
 fig.tight_layout();fig.savefig(p/'figures'/(sample+'_reviewed_atlas.png'),dpi=160,bbox_inches='tight');plt.close(fig)
 for label in sorted(cl.unique()):
  mask=(cl==label).to_numpy();n=int(mask.sum());lib=float(total[mask].sum());groups.append({'sample':sample,'patient_id':row['patient_id'],'cell_class':label,'cells':n,'total_UMI':lib})
  for name,panel in panels.items():
   idx=[genes.index(g) for g in panel if g in genes and (symbol==g).any()]
   if not idx:continue
   program_rows.append({'sample':sample,'patient_id':row['patient_id'],'site':row['site'],'cell_class':label,'cells':n,'program':name,'covered_genes':len(idx),'mean_log1p_panel':float(expr[mask][:,idx].mean()),'pct_cells_detecting_at_least_3_genes':float(((raw[mask][:,idx]>0).sum(1)>=3).mean()*100),'metric':'mean log1p normalized expression; not a control-subtracted score'})
  sums=raw[mask].sum(0);means=expr[mask].mean(0);pct=(raw[mask]>0).mean(0)*100
  for j,g in enumerate(genes):gene_rows.append({'sample':sample,'patient_id':row['patient_id'],'site':row['site'],'cell_class':label,'cells':n,'gene':g,'gene_present':bool((symbol==g).any()),'pct_detected':float(pct[j]),'mean_log1p':float(means[j]),'pseudobulk_count':float(sums[j]),'pseudobulk_CPM':float(sums[j]/lib*1e6) if lib else 0.})
 print('LINEAGE_PROGRAMS',sample,n,flush=True);del a,o,raw,expr,mapping;gc.collect()
pd.DataFrame(mural_rows).to_csv(out/'mural_pericyte_SMC_review.tsv.gz',sep='\t',index=False)
pr=pd.DataFrame(program_rows);ge=pd.DataFrame(gene_rows)
pr.to_csv(out/'sample_lineage_program_expression.tsv',sep='\t',index=False);ge.to_csv(out/'sample_lineage_candidate_genes.tsv',sep='\t',index=False);pd.DataFrame(groups).to_csv(out/'sample_lineage_comparison_cell_counts.tsv',sep='\t',index=False)
left=ge[(ge.cell_class=='CAF_strict')&(ge.cells>=20)];right=ge[(ge.cell_class=='Malignant_OS_candidate')&(ge.cells>=20)]
paired=left.merge(right,on=['sample','patient_id','gene'],suffixes=('_CAF','_malignant'));paired['CAF_minus_malignant_log2CPM']=np.log2(paired.pseudobulk_CPM_CAF+1)-np.log2(paired.pseudobulk_CPM_malignant+1);paired.to_csv(out/'patient_CAF_vs_malignant_candidate_genes.tsv',sep='\t',index=False)
if len(paired):
 summary=paired.groupby('gene').agg(paired_patients=('patient_id','nunique'),median_log2CPM_difference=('CAF_minus_malignant_log2CPM','median'),patients_CAF_higher=('CAF_minus_malignant_log2CPM',lambda x:int((x>0).sum())),median_CAF_detection=('pct_detected_CAF','median'),median_malignant_detection=('pct_detected_malignant','median')).reset_index();summary['interpretation']='patient-paired descriptive specificity screen, not protein or clinical validation';summary.to_csv(out/'CAF_vs_malignant_gene_specificity.tsv',sep='\t',index=False)
(out/'LINEAGE_COMPARISON_COMPLETE').write_text('LINEAGE_PROGRAM_COMPARISON_COMPLETE\n');print('LINEAGE_PROGRAM_COMPARISON_COMPLETE',len(groups),len(paired),flush=True)
