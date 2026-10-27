"""SYNTAX_CHECK: python -m py_compile code/freeze_review_report.py
RUN_COMMAND: python code/freeze_review_report.py
EXPECTED_SUCCESS_MARKER: REVIEW_AUDIT_READY_FOR_DAY5_6_ANALYSIS
"""
from pathlib import Path
import json,datetime,hashlib
import pandas as pd,numpy as np
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=r/'phase7_review_states_20260917';old=r/'phase6_nonpdx_cnv_20260916';delivery=p/'delivery';delivery.mkdir(exist_ok=True)
assert (p/'identity/COMPLETE').read_text().strip()=='PHASE6_IDENTITY_INTEGRATION_COMPLETE'
q=pd.read_csv(p/'identity/tool_reference_QC.tsv',sep='\t').fillna('');s=pd.read_csv(p/'identity/sample_identity_summary.tsv',sep='\t');before=pd.read_csv(old/'identity/sample_identity_summary.tsv',sep='\t');qc=pd.read_csv(r/'phase5_nonpdx_qc_20260916/results/sample_summary.tsv',sep='\t')
allcells=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t');strict=pd.read_csv(p/'identity/caf_strict_cells.tsv',sep='\t');ext=pd.read_csv(p/'identity/caf_extended_cells.tsv',sep='\t')
assert len(allcells)==125212 and allcells.cell_id.is_unique
assert len(strict)==int(s.caf_strict.sum())==int(allcells.caf_strict.sum())
assert len(ext)==int(s.caf_extended.sum())==int(allcells.caf_extended.sum())
assert set(strict.cell_id)<=set(ext.cell_id)
assert set(pd.read_csv(old/'identity/caf_strict_cells.tsv',sep='\t').cell_id)<=set(strict.cell_id)
checks=[]
for x in s.to_dict('records'):
 sm=x['sample'];d=allcells[allcells['sample']==sm];assert len(d)==int(x['retained'])
 assert int(d.caf_strict.sum())==x['caf_strict'] and int(d.caf_extended.sum())==x['caf_extended']
 for row in q[(q['sample']==sm)&(q.state=='complete')].to_dict('records'):
  metrics=json.loads(row['lineage_checks']);expected=bool(metrics) and all(v['evaluable']>=10 and v['normal_fraction']>=.9 for v in metrics.values()) and float(row['heldout_evaluable'])>=10 and float(row['heldout_normal_fraction'])>=.9
  assert bool(row['reference_check_pass'])==bool(expected),(sm,row['method'])
 if d.caf_strict.any():
  z=d[d.caf_strict];assert z.reference_QC_pass.all() and (z.n_genes>=500).all() and (z.det_Fib_review>=2).all() and (z.lineage_reviewed=='FibLike').all() and (z.cnv_concordance=='both_low').all() and (~z.CTA_context).all() and (z.cluster_high_fraction<.2).all()
 checks.append({'sample':sm,'cell_accounting':'PASS','reference_gate_recalculation':'PASS','strict_criteria':'PASS','strict_cells':int(d.caf_strict.sum()),'extended_cells':int(d.caf_extended.sum())})
# Cross-check the identity lists against an independent recalculation, not only completion flags.
(p/'review/final_identity_audit.json').write_text(json.dumps(checks,indent=2))
compare=before[['sample','reference_QC_pass','caf_strict','caf_extended']].merge(s[['sample','reference_QC_pass','caf_strict','caf_extended','selected_secondary_method']],on='sample',suffixes=('_before','_after'))
compare['strict_added']=compare.caf_strict_after-compare.caf_strict_before;compare['extended_change']=compare.caf_extended_after-compare.caf_extended_before
compare.to_csv(delivery/'复核前后候选集合对照.tsv',sep='\t',index=False)
s.to_csv(delivery/'样本身份与CAF数量.tsv',sep='\t',index=False);q.to_csv(delivery/'工具参考检查完整结果.tsv',sep='\t',index=False)
clinical=pd.read_csv(r/'phase4_download_20260914/clinical_source_audit.tsv',sep='\t',encoding='utf-8-sig').fillna('')
clinical.to_csv(delivery/'临床分组依据与缺失字段.tsv',sep='\t',index=False)
now=datetime.datetime.now(datetime.timezone.utc).astimezone(datetime.timezone(datetime.timedelta(hours=8)))
strict_samples=s[s.caf_strict>0];patients=strict.patient_id.unique().tolist();lung=strict[strict.site=='Lung'];nl=int(lung.patient_id.nunique());nb=int(strict.loc[strict.site!='Lung','patient_id'].nunique());npass=int(s.reference_QC_pass.sum())
raw_discord=int(before.CNV_discordant_or_ambiguous.sum());review_discord=int(s.CNV_discordant_or_ambiguous.sum());pending=int(s.Identity_pending.sum())
failed_samples=['GSE299023__OS.C.21.BX.OCT','GSE299023__OS.C.32.BX','GSE299023__OS.C.33.BX','GSE152048__BC3','GSE152048__BC6']
review_decisions=[]
for sm in failed_samples:
 row=s[s['sample']==sm].iloc[0];methods=q[(q['sample']==sm)&(q.state!='not_run')]
 issues=[]
 for z in methods.to_dict('records'):
  if z['state']!='complete':issues.append(z['method']+'运行失败');continue
  m=json.loads(z['lineage_checks'])
  for role,values in m.items():
   if values['evaluable']<10:issues.append(f"{z['method']} {role}仅{values['evaluable']}个可评估细胞")
   elif values['normal_fraction']<.9:issues.append(f"{z['method']} {role}正常比例{values['normal_fraction']:.1%}")
 decision='通过预设参考门槛，采用'+row.selected_secondary_method+'与inferCNV联合判定' if row.reference_QC_pass else '未通过双方法参考门槛，不纳入严格CAF集合'
 review_decisions.append({'sample':sm,'decision':decision,'strict_CAF':int(row.caf_strict),'extended_CAF':int(row.caf_extended),'issues_retained':'; '.join(issues)})
pd.DataFrame(review_decisions).to_csv(delivery/'五个重点样本复核结论.tsv',sep='\t',index=False)

print('REVIEW_AUDIT_READY_FOR_DAY5_6_ANALYSIS',len(strict),len(ext),npass,flush=True)
