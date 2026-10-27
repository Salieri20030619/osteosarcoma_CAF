from pathlib import Path
import pandas as pd,json,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase10_batch_gates_20260919'
for d in ['code','logs','tables','resources','inputs','decontx','objects','figures','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
for phase in ['phase3_identity_20260914/identity','phase7_review_states_20260917/identity','geo_identity_20260919/tables']:
 f=P/phase/'all_cell_identity_evidence.tsv.gz';d=pd.read_csv(f,sep='\t',nrows=3);print(phase);print(d.columns.tolist());print(d[[c for c in ['sample','short','dataset','lineage_reviewed','lineage_phase3','lineage_reviewed_candidate','cluster','cnv_concordance','joint','infercnv_call','decision_cnv2_source','decision_cnv2_call','copykat_call','reference_QC_pass','technology','technology_reported'] if c in d]].to_string(index=False))
f=P/'OS_CAF_phase2/results/GSM8337654/annotated.h5ad';a=ad.read_h5ad(f,backed='r');print('EARLY_FULL',a.shape,list(a.obs.columns),list(a.var.columns));a.file.close()
for d in [P/'phase7_review_states_20260917/remediation_20260917/decontx',P/'geo_identity_20260919/decontx']:print(str(d),[p.name for p in d.iterdir() if p.is_dir() and (p/'COMPLETE').exists()])
(O/'STATUS.json').write_text(json.dumps({'status':'UNIFIED_GATES_AND_BATCH_REMEDIATION_STARTED','original_stage':'phase9_joint_20260919','target_name_preserved':'SR1004','original_results_preserved':True},indent=2))