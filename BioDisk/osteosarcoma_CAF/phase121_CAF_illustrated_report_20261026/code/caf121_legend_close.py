from pathlib import Path
import json,hashlib,shutil,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase121_CAF_illustrated_report_20261026';q=json.loads((O/'resources/LEGEND_QA.json').read_text());p=O/'delivery/CAF_illustrated_report_legends.pdf';assert hashlib.sha256(p.read_bytes()).hexdigest()==q['pdf_sha256'];assert q['changed_pages']==[3,4,11] and q['other_16_pages_pixel_identical']
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
q.update(visual_QA='PASS: pages3,4,11 rendered with Poppler and visually checked; complete categorical legends, no clipping',status='completed_verified',checked_at=datetime.datetime.now(datetime.timezone.utc).isoformat());dump(O/'resources/LEGEND_QA.json',q)
f=O/'FINAL_STATUS.json';shutil.copyfile(f,O/'resources/before_legend_revision/FINAL_STATUS.json');x=json.loads(f.read_text());x.update(report=str(p),local_report='C:/Users/31937/Downloads/CAF_illustrated_report_20261026/骨肺CAF亚群_图文分析与Lung04推荐报告_图例补充版.pdf',pdf_sha256=q['pdf_sha256'],revision='colour_and_symbol_legends_added',revised_at=q['checked_at'],changed_pages=[3,4,11]);dump(f,x);dump(O/'STATUS.json',x)
for n in ['CAF_ILLUSTRATED_REPORT_CURRENT.json','CAF_CANDIDATE_RECOMMENDATION_CURRENT.json']:dump(R/n,x)
r=R/'CAF_RESUME_STATUS_CURRENT.json';y=json.loads(r.read_text());y.update(latest_task=x,last_completed_task=x,updated_at=q['checked_at']);dump(r,y)
fa=O/'FIGURES_COMPLETE.json';g=json.loads(fa.read_text());rev=json.loads((O/'resources/LEGEND_REVISION.json').read_text());g['figure_hashes'].update(rev['hashes']);g['legend_revision']=rev;dump(fa,g)
ac=O/'ACCEPTANCE.json';a=json.loads(ac.read_text());shutil.copyfile(ac,O/'resources/before_legend_revision/ACCEPTANCE.json');a.update(pdf_sha256=q['pdf_sha256'],local_pdf_sha256=q['pdf_sha256'],legend_QA=q);dump(ac,a)
dump(O/'LEGEND_DELIVERY_RECEIPT.json',x)
print(json.dumps(dict(status='completed_verified',pages=19,changed_pages=[3,4,11],pdf_sha256=x['pdf_sha256']),ensure_ascii=False))
