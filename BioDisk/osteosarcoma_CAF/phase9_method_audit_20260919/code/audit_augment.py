from pathlib import Path
import inspect,json,hashlib
src=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase9_method_audit_20260919.py').read_text()
exec(src.split('baseline=gate()')[0])
base=e.caf_supported_candidate.astype(bool);rows=[]
for cutoff in [.1,.2,.5,1.0,float('inf')]:
 raw=fib&e.joint.eq('both_low')&(e.det_Fib>=2)&~mixed(e)&~e.CTA_context.astype(bool)&(e.det_CTA<2)&~comp(e,osteo=False)&(e.cluster_both_high_fraction.fillna(1)<.2)&(e.cluster_any_high_fraction.fillna(1)<.2)&(e.n_genes>=500)&e.reference_QC_pass.astype(bool)
 cor=(e.corrected_det_Fib>=2)&~comp(e,'corrected_',osteo=False)&~mixed(e,'corrected_')&(e.corrected_det_CTA<2)&e.ambient_checked.astype(bool)
 mask=raw&cor&~((e.det_Osteo>=3)&(e.pct_Osteo>cutoff))&~((e.corrected_det_Osteo>=3)&(e.corrected_pct_Osteo>cutoff))
 if cutoff==.2:assert (mask==base).all()
 for ds in ['ALL']+sorted(e.dataset.unique()):
  m=np.ones(len(e),bool) if ds=='ALL' else e.dataset.eq(ds).to_numpy();rows.append({'osteogenic_percent_threshold':cutoff,'dataset':ds,'candidates':int((mask&m).sum()),'added':int((mask&~base&m).sum()),'removed':int((~mask&base&m).sum()),'status':'diagnostic_only_not_admitted'})
save(pd.DataFrame(rows),'osteogenic_gate_separate_audit.tsv')
(O/'resources/harmonypy_run_harmony_source.py').write_text(inspect.getsource(harmonypy.run_harmony))
old=P/'phase7_review_states_20260917/code/integrate_identity.py';(O/'resources/legacy_GSE152048_identity.py').write_text(old.read_text())
log=(J/'logs/integration.log').read_text();(O/'resources/harmony_convergence_lines.txt').write_text('\n'.join(x for x in log.splitlines() if 'converg' in x.lower()))
for name in ['code/run_cnv.R','result_review_20260919/tables/calibration_by_role.tsv','result_review_20260919/tables/parameters.tsv']:
 p=P/'geo_cnv_20260918'/name;(O/'resources'/p.name).write_bytes(p.read_bytes())
print(pd.DataFrame(rows).query('dataset == "ALL"').to_string(index=False));print('AUGMENT_COMPLETE')