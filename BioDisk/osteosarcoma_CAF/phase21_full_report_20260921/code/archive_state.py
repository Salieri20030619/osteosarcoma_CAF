from pathlib import Path
import json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase21_full_report_20260921'
S=P/'phase22_no_UMI_sensitivity_20260921'
pdf=next((O/'delivery').glob('*.pdf'))
assert (S/'CNV611/ANALYSIS_COMPLETE').exists()
cnv=json.loads((S/'CNV611/delivery/QA.json').read_text())
assert cnv['all_methods_complete']
qc=json.loads((S/'delivery/QA.json').read_text())
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
qa={'completed_utc':now,'PDF':pdf.name,'PDF_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'PDF_pages':15,'embedded_figures':12,'visual_QA':'all 14 original pages inspected; final changed pages 14 and 15 inspected; unchanged pages pixel-verified','main_analysis_unchanged':True,'no_UMI_scope':'all 52 QC recovery audit; SCPCS000611 dual CNV and DecontX rerun; no global reclustering','no_UMI_QC':qc,'no_UMI_CNV611':cnv,'client_message_sent':False}
(O/'delivery/QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
(S/'delivery/FINAL_QA.json').write_text(json.dumps({'QC_audit':qc,'CNV611':cnv,'completed_utc':now},indent=2))
state=P/'PROJECT_STATE.json'
backup=P/('PROJECT_STATE.before_full_report_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S')+'.json')
shutil.copy2(state,backup)
j=json.loads(state.read_text());j['phase21_full_report_20260921']=qa
j['phase22_no_UMI_sensitivity_20260921']={'status':'QC_AUDIT_AND_CNV611_COMPLETE','scope':qa['no_UMI_scope'],'QC':qc,'CNV611':cnv,'legacy_added_doublets':'not_assessed; not counted as confirmed singlets','global_reclustering':'not_performed'}
j['current_user_priority']='Full main-line report with embedded figures and no-UMI supplementary analysis; PDF-only delivery; small-scale analytical flow'
tmp=state.with_suffix('.json.tmp');tmp.write_text(json.dumps(j,ensure_ascii=False,indent=2));tmp.replace(state)
print(json.dumps({'PDF_sha256':qa['PDF_sha256'],'project_state_backup':str(backup)}))
