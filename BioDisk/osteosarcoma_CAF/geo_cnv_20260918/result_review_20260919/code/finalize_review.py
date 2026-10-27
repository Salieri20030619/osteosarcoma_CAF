from pathlib import Path
import pandas as pd,json,datetime,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918');O=P/'result_review_20260919';T=O/'tables'
def md(cols,rows):return '\n'.join(['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
s=pd.read_csv(T/'sample_summary.tsv',sep='\t');z=pd.read_csv(T/'cell_review.tsv.gz',sep='\t');a=pd.read_csv(T/'calibration_by_role.tsv',sep='\t');r=(O/'code/report_template.md').read_text(encoding='utf-8').replace('?????????','?????????')
cal=[];sam=[]
for x in s.to_dict('records'):
 sm=x['sample'].split('__')[1]
 cal.append([sm,x['holdout'],x['holdout_infer_high'],x['holdout_copy_high'],x['holdout_copy_ne'],f"{x['holdout_copy_ne']/x['holdout']*100:.2f}%"])
 sam.append([sm,x['targets'],x['target_binary_n'],f"{x['target_binary_agreement']*100:.2f}%",x['target_copy_ne'],x['fib_cells'],x['fib_both_low'],x['fib_both_high']])
r=r.replace('{{CALIBRATION_TABLE}}',md(['??','???','infer?','Copy?','Copy??','Copy????'],cal)).replace('{{SAMPLE_TABLE}}',md(['??','????','?????','??????','????','FibLike','Fib??','Fib??'],sam))
t=z[z.role_group=='target'];vc=t.joint.value_counts();agree=int(vc.get('both_low',0)+vc.get('both_high',0));binary=agree+int(vc.get('infer_low_copy_high',0)+vc.get('infer_high_copy_low',0));txt=f"91,352 ????????? {vc['both_low']:,}??? {vc['both_high']:,}?infer?/Copy? {vc['infer_low_copy_high']:,}?infer?/Copy? {vc['infer_high_copy_low']:,}?Copy???? {vc['copy_not_evaluable']:,}?infer ambiguous {vc['infer_ambiguous']:,}???????? {agree:,}/{binary:,}?{agree/binary*100:.2f}%????????????????????? {agree/len(t)*100:.2f}%?????????????????????????????"
r=r.replace('{{TARGET_SUMMARY}}',txt)
# Resolve the fixed author-code path from the downloaded source tree.
tree=json.loads((P/'sources/author_repository_tree.json').read_text());paths=[q['path'] for q in tree.get('tree',[]) if q['path'].endswith('02_infercnv_sample_list.r')]
if len(paths)==1:
 r=r.replace('06_major_celltype_annotation/02_infercnv_sample_list.r',paths[0])
 print('AUTHOR_SOURCE_PATH',paths[0])
assert '{{' not in r
(O/'CAF_GEO_CNV_review_20260919.md').write_text(r,encoding='utf-8')
# Verify conservation and essential review outputs, without rerunning upstream analyses.
assert len(z)==97238 and not z.cell_id.duplicated().any()
assert int(s.fib_cells.sum())==8489 and int(s.fib_both_low.sum())==3028
assert int(s.targets.sum())==91352 and int(s.holdout.sum())==1190
assert (pd.read_csv(T/'integrity_checks.tsv',sep='\t').iloc[:,2:6]==True).all().all()
refs=pd.read_csv(T/'reference_gate_checks.tsv',sep='\t');assert len(refs)==15 and refs.reference_rules_pass.all() and refs.no_holdout_baseline_overlap.all()
assert (O/'logs/profile_export.exit_code').read_text().strip()=='0' and (O/'logs/profile_review.exit_code').read_text().strip()=='0'
assert 'RESULT_TABLE_REVIEW_COMPLETE' in (O/'logs/review.log').read_text()
assert 'PROFILE_REVIEW_COMPLETE' in (O/'logs/profile_review.log').read_text()
for f in ['cnv_review_overview.png','lung_chromosome_profiles.png']:assert (O/'figures'/f).stat().st_size>10000
state=dict(review_status='COMPLETE_IDENTITY_ACCEPTANCE_OPEN',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),samples=15,cells=97238,computational_integrity='PASS',reference_rule_audit='PASS',figure_visual_QA='PASS',biological_identity_acceptance='OPEN',confirmed_new_CAF_count=None,labels_changed=False,original_project_state_changed=False,new_CNV_jobs_started=False,inferential_tests_performed=False,report='CAF_GEO_CNV_review_20260919.md',findings=['GSM9504303 dominant FibLike population withheld from CAF inclusion: cross-method call discordance with shared chromosome signals','GSM9504304: 431 dual-low FibLike candidates, not confirmed CAF','1190 heldout normals: CopyKAT 1070 diploid, 4 aneuploid, 116 not evaluable','FLEX cohort CopyKAT missing calls 11732/42686; no independent sensitivity validation','GSM9504299 remains unassessed due insufficient normals'],log_note='Mean-of-empty-slice warnings in descriptive medians correspond to groups with all missing CopyKAT signal; NA is retained, not imputed.')
(O/'REVIEW_STATUS.json').write_text(json.dumps(state,indent=2)+'\n')
(O/'REVIEW_COMPLETE').write_text('RESULT_REVIEW_COMPLETE_IDENTITY_ACCEPTANCE_OPEN\n')
print('REPORT_AND_REVIEW_ACCEPTANCE_COMPLETE')
