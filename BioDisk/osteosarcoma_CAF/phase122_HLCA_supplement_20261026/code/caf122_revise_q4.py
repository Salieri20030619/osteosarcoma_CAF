from pathlib import Path
import json,shutil,hashlib
W=Path('C:/Users/31937/caf122_hlca_supplement_20261026/unified_full_report/qa_explicit_methods/q4_revision');W.mkdir(exist_ok=True)
path=Path('C:/Users/31937/caf122_reviewer_qa_explicit.json');old=json.loads(path.read_text(encoding='utf-8'));new=json.loads(json.dumps(old));pdf=Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答.pdf')
if not (W/'previous_QA.pdf').exists():shutil.copyfile(pdf,W/'previous_QA.pdf');shutil.copyfile(path,W/'previous_answers.json')
new[3][1]=[
'本研究目前只有7位有效配对患者，属于小样本探索研究，主要目标是识别值得优先验证的候选。Lung_04在基质、胶原交联、纤维型胶原和FNDC1程序上均呈7/7同向增强，说明其推荐依据来自多个患者一致的功能方向，而不是单个患者的突出表现。',
'基质和FNDC1面板的平均效应分别为+0.849和+0.318，原始P均为0.015625，对应FDR为0.08984和0.06875。小样本限制了效应估计的精度和检验能力；现有数据的校正后证据尚未达到0.05阈值。因此，我们将其作为有一致方向支持的名义显著、探索性候选推进，以增加独立患者样本进一步检验其重复性，而不因FDR未过就直接舍弃这条线索。'
]
assert all(new[i]==old[i] for i in [0,1,2,4]);assert '敏感' not in json.dumps(new,ensure_ascii=False);path.write_text(json.dumps(new,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');(W/'Q4_EDIT_RECORD.json').write_text(json.dumps(dict(changed_question='Q4 only',unchanged_questions=['Q1','Q2','Q3','Q5'],old=old[3],new=new[3],previous_QA_sha256=hashlib.sha256((W/'previous_QA.pdf').read_bytes()).hexdigest()),ensure_ascii=False,indent=2),encoding='utf-8');print('Q4_ONLY_UPDATED')
