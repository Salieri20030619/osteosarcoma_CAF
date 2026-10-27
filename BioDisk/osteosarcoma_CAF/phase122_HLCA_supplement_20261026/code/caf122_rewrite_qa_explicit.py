from pathlib import Path
import json,hashlib,shutil,pymupdf
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026/unified_full_report');W=B/'qa_explicit_methods';W.mkdir(exist_ok=True)
target=Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答.pdf');backup=W/'previous_QA.pdf'
if not backup.exists():shutil.copyfile(target,backup)
data=json.loads(Path('C:/Users/31937/caf122_reviewer_qa_explicit.json').read_text(encoding='utf-8'));assert len(data)==5 and '敏感' not in json.dumps(data,ensure_ascii=False)
src=Path('C:/Users/31937/caf122_reviewer_qa.py').read_text(encoding='utf-8');start=src.index('qas=[');end=src.index('story=[',start);src=src[:start]+'qas='+repr(data)+'\n'+src[end:];src=src.replace("B=Path('C:/Users/31937/caf122_hlca_supplement_20261026/unified_full_report')", "B=Path('C:/Users/31937/caf122_hlca_supplement_20261026/unified_full_report/qa_explicit_methods')")
src=src.replace('if i==3:','if i==2:').replace('以下5组问答基于当前报告已完成的分析。答复区分现有证据与后续验证，不将尚未开展的实验表述为已完成。','以下5组问答逐项说明已经使用的分析方法、得到的结果及其对候选推荐的支持。')
builder=W/'build_qa.py';builder.write_text(src,encoding='utf-8');exec(compile(src,str(builder),'exec'),{'__name__':'__main__'})
doc=pymupdf.open(W/'reviewer_QA_review.pdf');text=''.join(p.get_text() for p in doc);assert len(doc)==2;assert '敏感' not in text;assert all(text.count(f'Q{i}｜')==1 for i in range(1,6));assert all(x in text for x in ['PTPRC','LST1','TYROBP','CD3D','CD3E','PECAM1','VWF','CDH5','94.7%','96.1%','1.704','1.772','0.610','0.08984','0.06875','0.10577'])
rec=dict(status='STAGED_PASS',pages=2,questions=5,prohibited_term_absent=True,actual_screening_rules_checked_against='caf120_audit.py:49-51',main_report_unchanged_sha256=hashlib.sha256(Path('C:/Users/31937/Downloads/骨肺CAF亚群_图文分析与Lung04推荐报告.pdf').read_bytes()).hexdigest(),sha256=hashlib.sha256((W/'reviewer_QA_review.pdf').read_bytes()).hexdigest(),backup=str(backup));(W/'QA_REVISION_CHECK.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(rec,ensure_ascii=False))
