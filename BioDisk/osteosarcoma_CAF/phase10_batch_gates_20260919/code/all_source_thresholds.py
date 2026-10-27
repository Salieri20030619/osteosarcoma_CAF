from pathlib import Path
import pandas as pd,numpy as np,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919';e=pd.read_csv(O/'tables/all_cell_unified_evidence.tsv.gz',sep='\t',index_col='cell_id',low_memory=False)
raw=e[[c for c in e if c.startswith('raw_')]].rename(columns=lambda c:c[4:]);cor=e[[c for c in e if c.startswith('corrected_')]].rename(columns=lambda c:c[10:])
src=(O/'code/apply_quality_guard.py').read_text();start=src.index('def conflict(');end=src.index('\nfixed=',start);exec(src[start:end])
def molecular(ng=500,fib=2,risk=.2,frac=.2,osteo=.2):
 fixed=e.dual_low&e.reference_QC_pass&~e.CTA_context&~e.independent_doublet_excluded&(raw.det_CTA<2)&(raw.det_Fib>=fib)&(e.n_genes>=ng)&(e.cluster_both_high_fraction.fillna(1)<risk)&(e.cluster_any_high_fraction.fillna(1)<risk)&e.ambient_checked&(e.corrected_library_counts>0)
 return fixed&~conflict(raw,frac,osteo)&~conflict(cor,frac,osteo)&(cor.det_Fib>=fib)&(cor.det_CTA<2)
base=molecular()&e.core_scope;assert np.array_equal(base,e.unified_strict)
conditions=[('baseline','.2',{})]+[('n_genes',str(v),{'ng':v}) for v in [300,800,1000]]+[('Fib_genes',str(v),{'fib':v}) for v in [1,3]]+[('cluster_CNV_fraction',str(v),{'risk':v}) for v in [.1,.3,.5]]+[('competing_fraction_percent',str(v),{'frac':v,'osteo':v}) for v in [.1,.5]]+[('osteo_fraction_percent',str(v),{'osteo':v}) for v in [.1,.5,1.0]]
rows=[];controls=[]
for name,value,kw in conditions:
 mol=molecular(**kw);mask=mol&e.core_scope
 for ds in ['ALL']+sorted(e.dataset.unique()):
  m=np.ones(len(e),bool) if ds=='ALL' else e.dataset.eq(ds).to_numpy();rows.append({'varied_rule':name,'value':value,'dataset':ds,'candidates':int((mask&m).sum()),'added_vs_unified':int((mask&~base&m).sum()),'removed_vs_unified':int((~mask&base&m).sum()),'old_joint_retained':int((mask&e.old_joint_member&m).sum()),'scope':'diagnostic only; unchanged primary'})
 for role in ['holdout_immune','holdout_endo']:
  rr=e.cnv_role.eq(role);controls.append({'varied_rule':name,'value':value,'role':role,'cells':int(rr.sum()),'molecular_pass_without_scope':int((rr&mol).sum())})
 if name=='cluster_CNV_fraction':e.loc[mask&~base,['dataset','sample','original_cluster','cluster_any_high_fraction','cluster_both_high_fraction','old_joint_member','exclusion_reasons']].to_csv(O/'tables'/('CNV_risk_'+value+'_additional_candidates.tsv.gz'),sep='\t',index_label='cell_id')
pd.DataFrame(rows).to_csv(O/'tables/all_source_threshold_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(controls).to_csv(O/'tables/all_source_threshold_control_checks.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).query('dataset == "ALL"').to_string(index=False));(O/'ALL_SOURCE_THRESHOLD_COMPLETE').write_text('ALL_SOURCE_THRESHOLD_SENSITIVITY_COMPLETE\n')