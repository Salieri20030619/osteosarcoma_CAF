from pathlib import Path
import json,hashlib,datetime,shutil
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'report_update_allclusters_GSEA_20261006';R=O/'resources';R.mkdir(exist_ok=True)
receipt=json.loads((O/'delivery/DELIVERY_RECEIPT.json').read_text())
report=next((O/'delivery').glob('*.docx'));z=next((O/'delivery').glob('*.zip'))
assert hashlib.sha256(report.read_bytes()).hexdigest()==receipt['report_sha256']
assert hashlib.sha256(z.read_bytes()).hexdigest()==receipt['zip_sha256']
task={'title':'胶原与基质报告全亚群记录及GSEA整合','status':'completed_verified_delivered','report':str(report),'local_report':receipt['report'],'results_zip':str(z),'local_results_zip':receipt['zip'],'pages':receipt['pages'],'figures':7,'completed_at':receipt['completed_at'],'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'remaining_required_work':[],'scope':'33亚群 全来源54条名义上调记录 39基因 477条既有GSEA跟进记录 保留模型统计与患者统计边界 未新增推断性检验'}
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','PROJECT_STATE.json']:
 p=B/name;backup=R/('BEFORE_REPORT_UPDATE_'+name)
 if not backup.exists():shutil.copy2(p,backup)
 d=json.loads(p.read_text());d['latest_report_revision']=task;d['latest_report']=str(report);d['report']=str(report);d['local_report']=receipt['report'];d['latest_report_revision_completed_at']=receipt['completed_at'];d['last_activity_at']=receipt['completed_at'];d['all_subgroup_GSEA_report_summary']={'nominal_gene_source_rows':54,'unique_genes':39,'nominal_subgroups':7,'GSEA_rows':477,'GSEA_all_contexts':33,'GSEA_Fibro_contexts':20,'new_inferential_tests':False,'historical_phase68_gene_statistics_unchanged':True}
 tmp=p.with_name(p.name+'.report_update.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n');tmp.replace(p)
(O/'REPORT_DELIVERED.json').write_text(json.dumps(task,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'passed':True,'report':str(report),'zip':str(z),'progress_updated':3},ensure_ascii=False))
