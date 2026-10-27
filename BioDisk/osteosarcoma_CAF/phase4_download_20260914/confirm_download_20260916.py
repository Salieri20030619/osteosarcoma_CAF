"""SYNTAX_CHECK: python -m py_compile confirm_download_20260916.py
RUN_COMMAND: python confirm_download_20260916.py
EXPECTED_SUCCESS_MARKER: CAF_DOWNLOAD_AUDIT_PASS
"""
from pathlib import Path
import json,hashlib,gzip,csv,datetime
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914')
assert (r/'resume_20260916.exit_code').read_text().strip()=='0'
assert (r/'COMPLETE').read_text().strip()=='NON_PDX_DOWNLOAD_VALIDATION_COMPLETE'
manifest=json.loads((r/'download_manifest.json').read_text())
checks=json.loads((r/'download_checks.json').read_text())
rows=list(csv.DictReader((r/'source_and_matrix_validation.tsv').open(),delimiter='\t'))
expected={(s['dataset'],u.rsplit('/',1)[-1]) for s in manifest for u in s['urls']}
assert len(checks)==len(expected)==29
assert {(x['dataset'],Path(x['file']).name) for x in checks}==expected
old=json.loads((r/'download_status.json').read_text())
known={x['name']:x['expected_bytes'] for x in old.get('remaining',[]) if x.get('expected_bytes')}
for x in checks:
 p=Path(x['file']);assert p.is_file() and p.stat().st_size==x['bytes']
 if p.name in known:assert p.stat().st_size==known[p.name],p.name
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 assert h.hexdigest()==x['sha256'],p.name
 with gzip.open(p,'rb') as f:
  while f.read(4*1024*1024):pass
 assert x['status']=='download_integrity_pass'
assert len(rows)==len(manifest)==17
assert {(x['dataset'],x['sample']) for x in rows}=={(x['dataset'],x['sample']) for x in manifest}
assert all(x['validation']=='PASS' and x['CAF_confirmed']=='False' for x in rows)
assert all(x['direct_patient_verified'] and x['is_xenograft_verified'] is False for x in manifest)
report={'status':'CAF_DOWNLOAD_AUDIT_PASS','confirmed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files_verified':len(checks),'total_bytes':sum(x['bytes'] for x in checks),'tissue_matrices_validated':len(rows),'input_barcodes':sum(int(x['input_barcodes']) for x in rows),'new_cohort_QC_clustering_CNV_started':False,'next_state':'PAUSED_AFTER_DOWNLOAD_VALIDATION','cohorts':{d:sum(x['dataset']==d for x in rows) for d in sorted({x['dataset'] for x in rows})}}
(r/'completion_audit_20260916.json').write_text(json.dumps(report,indent=2))
status={'files_total':29,'files_verified':29,'verified_bytes':report['total_bytes'],'remaining':[],'process_exists':False,'status':'PAUSED_AFTER_DOWNLOAD_VALIDATION','matrix_validation_complete':True,'audit_file':'completion_audit_20260916.json','updated_at_utc':report['confirmed_at_utc']}
(r/'download_status.json').write_text(json.dumps(status,indent=2))
(r/'下载完成确认_20260916.md').write_text('# CAF 非 PDX 队列下载完成确认\n\n确认时间：'+report['confirmed_at_utc']+'。\n\n29/29 个文件已完成，合计 '+str(report['total_bytes'])+' 字节。全部文件重新核对 SHA256 并完整读取 gzip 流；17/17 个组织矩阵通过维度、非负整数计数、条形码和基因标识检查。\n\nGSE299023：6 个组织；GSE152048：11 个组织。组织数及条形码数不作为独立患者数。新队列 QC、聚类、CAF 注释与 CNV 未启动。下载与验证后暂停。\n\n完成证据：completion_audit_20260916.json、download_checks.json、source_and_matrix_validation.tsv、download_resume_20260916.log、COMPLETE、PAUSED_AFTER_DOWNLOAD_VALIDATION。9 月 14 日暂停记录保留为历史；当前状态以本次确认和 download_status.json 为准。\n',encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
