from pathlib import Path
import json,hashlib,shutil
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'report_merge_consistency_collagen_GSEA_20261006';R=O/'resources';R.mkdir(exist_ok=True)
r=json.loads((O/'delivery/FINAL_ACCEPTANCE.json').read_text())
report=next((O/'delivery').glob('*.docx'))
assert hashlib.sha256(report.read_bytes()).hexdigest()==r['report_sha256'] and r['passed']
task={'title':'CAF功能分类一致性报告整合胶原基质与GSEA证据','status':'completed_verified_delivered','report':str(report),'local_report':r['report'],'base_report':r['base_report'],'incorporated_report':r['added_report'],'pages':26,'figures':16,'tables':14,'completed_at':r['completed_at'],'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'remaining_required_work':[],'scope':'以用户明确指定的分类一致性报告为主 保留UMAP组成及33亚群分类附表 整合54条胶原基质记录和GSEA 统一摘要图表编号结论 两份原报告保留','statistics_archive':str(B/'report_update_allclusters_GSEA_20261006/delivery/CAF全亚群与GSEA整合报告_配套统计_20261006.zip')}
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','PROJECT_STATE.json']:
 p=B/name;back=R/('BEFORE_CONSISTENCY_MERGE_'+name)
 if not back.exists():shutil.copy2(p,back)
 d=json.loads(p.read_text());d.update({'latest_report_revision':task,'latest_report':str(report),'report':str(report),'local_report':r['report'],'latest_report_revision_completed_at':r['completed_at'],'last_activity_at':r['completed_at']})
 tmp=p.with_name(p.name+'.merge_report.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n');tmp.replace(p)
(O/'REPORT_DELIVERED.json').write_text(json.dumps(task,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'passed':True,'report':str(report),'pages':26,'progress_updated':3},ensure_ascii=False))
