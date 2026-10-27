from pathlib import Path
import json,hashlib,pymupdf
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026/unified_full_report')
for name,file in [('YH','msyh.ttc'),('YHB','msyhbd.ttc')]:pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
styles={k:ParagraphStyle(k,fontName='YHB' if k!='body' else 'YH',fontSize=s,leading=l,spaceAfter=g,wordWrap='CJK',textColor=colors.HexColor('#243b49')) for k,s,l,g in [('title',18,27,14),('q',12,19,9),('body',10.5,17.5,10),('cite',8.5,13,14)]}
qas=[
('Q1｜为什么用7位患者推荐Lung_04？11、7和5这些人数为何不同？',
 ['11是跨来源至少检出1个Lung_04细胞的患者记录数（NCI 10条、SCPCP000017 1条）；7是NCI内主分析目标与同患者对照各至少5细胞的有效配对患者数；5是严格身份层下仍满足相同配对条件的人数。细胞数和文库数不作为独立生物学重复。',
 '本项目将6位设为操作性下限，依据之一是双侧精确符号翻转检验：N=5的最小P为0.0625，N=6为0.03125。该下限不是充分功效保证。Lung_04的7位患者在多个核心程序上同向，故适合作为当前队列的优先验证候选；证据仍属小样本探索，后续需要独立外部患者验证。'],
 '对应证据：正文第二节、图8–11、表4–5。'),
('Q2｜如何排除周细胞、平滑肌或其他细胞混入？diploid能证明CAF身份吗？',
 ['本研究联合标志物表达、成纤维支持层、NCI作者CopyKAT标注和明确免疫/内皮组合共表达排除进行敏感性分析，并使用HLCA健康核心图谱作外部身份参照。Lung_04在可评估患者中的成纤维细胞样比例平均为94.7%，冻结参考参数延长训练后身份大类一致率为96.1%。',
 '这些结果支持候选以成纤维细胞样为主，但不能证明绝对纯度。diploid并非CAF身份或非恶性身份的充分条件，HLCA映射也不是实验鉴定。Bone_07及其子群仍保留平滑肌样混合信息。后续应结合组织原位共定位或独立分选，确认成纤维身份并检查竞争谱系标志物。'],
 '对应证据：图3–5、图11–12，正文第六、七节。'),
('Q3｜正常肺本来就表达胶原，如何说明候选超出正常背景？会不会只是批次差异？',
 ['主分析首先在同来源、同病灶、同患者内部比较目标亚群与其余候选细胞，提供队列内相对增强的证据。HLCA再以健康肺泡、外膜和肿瘤邻近肺成纤维细胞提供独立组织背景；其细胞不参与原队列分群，供者数也不增加候选的有效配对人数。',
 'Lung_04的基质、纤维型胶原和交联分数高于健康肺参考，方向在肿瘤邻近肺及细胞门槛敏感性下保持。跨数据集仍存在平台、测序深度和过滤差异，因此外部比较仅报告分数差与百分位，不作跨数据集显著性检验，也不能据此单独认定肿瘤特异性。推荐基础仍是队列内的患者配对证据。'],
 '对应证据：图9、图13及表4–5；正文第六、七节。'),
('Q4｜为什么原始P显著而FDR未通过，仍然保留并推荐这个群？',
 ['本研究按预设目标程序识别原始P<0.05的名义显著线索，同时完整报告BH-FDR，并结合效应方向、幅度和患者一致性决定后续验证优先级。Lung_04基质程序的效应为+0.849，7/7同向，P=0.015625、FDR=0.08984；FNDC1面板效应+0.318，7/7同向，P=0.015625、FDR=0.06875。它们应称为名义显著或探索性信号，不能称为校正后显著。',
 'FDR在预先定义的同类终点家族内计算；文献audit面板与program等不同家族分开，不能将最有利的P值替代目标效应。当前推荐表示“优先验证”，不代表机制已确证。严格身份层N=5时，即使5/5同向，精确检验也只能达到P=0.0625，这一结果连同效应和小样本定位保留。'],
 '对应证据：表4、正文第七节统计方法。'),
('Q5｜这是否只是COL1A1等胶原基因升高？能直接称为促纤维沉积或典型myCAF吗？',
 ['证据不仅包括结构胶原，还包括基质、交联、FNDC1及不含COL结构基因的胶原加工面板。后者在主分析中6/7患者同向，平均效应+0.228，P=0.03125、FDR=0.10577；严格身份层5/5同向。它提供了结构胶原以外的补充线索，但仍属于转录层证据。',
 '单基因并非全部同向增强：Lung_04的ACTA2、POSTN和LOXL2患者配对结果没有显示一致增强，且典型myCAF程序与健康肺泡背景差值仅+0.036。因此采用“ECM/胶原生成相关CAF候选”，不直接命名为典型收缩型myCAF，也不将RNA高分等同于真实纤维沉积。胶原分泌、组织沉积和交联机制需分别通过蛋白、原位形态及相应功能实验验证。'],
 '对应证据：表4、正文第八节、补充表S1；探索性疾病状态比较见S2–S3。')]
story=[Paragraph('骨肺CAF研究<br/>审稿人可能问题与建议答复',styles['title']),Paragraph('以下5组问答基于当前报告已完成的分析。答复区分现有证据与后续验证，不将尚未开展的实验表述为已完成。',styles['body'])]
for i,(q,answers,cite) in enumerate(qas):
 if i==3:story.append(PageBreak())
 story.append(Paragraph(q,styles['q']))
 for answer in answers:story.append(Paragraph(answer.replace('<','&lt;'),styles['body']))
 story.append(Paragraph(cite,styles['cite']))
def footer(c,d):
 c.setStrokeColor(colors.HexColor('#176c7a'));c.line(42,810,553,810);c.setFont('YH',8);c.setFillColor(colors.HexColor('#176c7a'));c.drawString(42,819,'骨肺CAF研究｜独立审稿人问答');c.drawString(42,25,'2026-10-27');c.drawRightString(553,25,str(d.page))
out=B/'reviewer_QA_review.pdf';SimpleDocTemplate(str(out),pagesize=(595,842),leftMargin=42,rightMargin=42,topMargin=47,bottomMargin=43).build(story,onFirstPage=footer,onLaterPages=footer)
doc=pymupdf.open(out)
for i,p in enumerate(doc):p.get_pixmap(matrix=pymupdf.Matrix(1.15,1.15)).save(B/f'qa_page_{i+1}.png')
text=''.join(p.get_text() for p in doc);assert all(text.count(f'Q{i}｜')==1 for i in range(1,6));(B/'REVIEWER_QA_CONTENT.json').write_text(json.dumps(qas,ensure_ascii=False,indent=2),encoding='utf-8');print('QA_PAGES',len(doc),hashlib.sha256(out.read_bytes()).hexdigest())
