from pathlib import Path
import json,hashlib,shutil,io,re
import pymupdf as fitz
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.lib import colors
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'expanded_QA_revision';W.mkdir(exist_ok=True);target=Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答_扩充版.pdf');backup=W/'original_expanded_QA.pdf'
if not backup.exists():shutil.copyfile(target,backup)
assert hashlib.sha256(backup.read_bytes()).hexdigest()=='821ca7507da6adc83ad9da03530be19ae4d213740fadb015d08b1f052eab8918'
doc=fitz.open(backup);plans={};records=[]
pdfmetrics.registerFont(TTFont('YH','C:/Windows/Fonts/msyh.ttc'));pdfmetrics.registerFont(TTFont('YHB','C:/Windows/Fonts/msyhbd.ttc'));pdfmetrics.registerFontFamily('YH',normal='YH',bold='YHB',italic='YH',boldItalic='YHB')
def oldtext(page,indices):return ''.join(doc[page-1].get_text('blocks')[i][4] for i in indices)
def edit(page,indices,prefix,text,issue,size=9.5):
 bs=doc[page-1].get_text('blocks');old=oldtext(page,indices);assert re.sub(r'\s','',old).startswith(re.sub(r'\s','',prefix)),(page,indices,old[:100]);r=fitz.Rect(bs[indices[0]][:4])
 for i in indices[1:]:r|=fitz.Rect(bs[i][:4])
 r.x1=550;r.y0-=.7;r.y1+=3
 plans.setdefault(page,[]).append((r,text,size));records.append(dict(issue=issue,page=page,indices=indices,before=old,after=text,rect=list(r)))
approved=json.loads(Path('C:/Users/31937/caf122_reviewer_qa_explicit.json').read_text(encoding='utf-8'))
edit(2,list(range(1,7)),'基质和 FNDC1',approved[3][1][1],'QA-03: restore approved Q4 small-sample explanation')
q6='Lung_04表现为成纤维细胞样身份并伴成骨相关转录特征。成骨程序在主分析（+0.328，7/7同向，P=0.015625）、筛选后细胞（+0.409，5/5同向）及健康肺泡参照（+0.800）中均升高；矿化程序相对健康肺泡仅+0.072，Bone_07_s1为−0.039，因此不直接命名为已证实矿化的细胞群。<br/>这里使用两项独立计算：按作者diploid及免疫/内皮共检出规则筛选细胞后，重新计算患者配对效应，仍有5位有效患者，基质、纤维型胶原和FNDC1均5/5同向；另将6,423个主分析细胞映射至HLCA健康核心图谱（584,944细胞），在候选至少5细胞的7位可评估患者中，平均成纤维细胞样比例为94.7%，训练200至500轮后标签一致率为96.1%。<br/>健康图谱用于比较已知细胞身份，不能凭其未包含恶性类型排除恶性身份。文献[3]（GSE162454）也报道CAF的骨化相关程序。综合结果支持“伴成骨相关转录特征”的定位，进一步通过独立样本原位共定位确认成纤维身份及基质蛋白表达。'
edit(2,list(range(20,31)),'外部健康肺图谱',q6,'QA-01: separate filtered-cell contrasts from full-query mapping')
q8='可能的原因是，不同父群的患者构成和表达背景不同，使跨父群Z-score与同患者配对差值反映不同层面的差异。例如，Lung_00、Lung_01的部分标志物较低，可能使Lung_04的图内相对颜色偏红；而采用同患者自身对照时，上述基因未呈一致增强。气泡图用于描述表达，患者内检验用于判断相对效应。'
edit(3,[13,14,15],'由于',q8,'QA-05: state patient-composition explanation as possible')
q10='一个可能的原因是各来源的患者覆盖及细胞构成不同：NCI覆盖更多患者，其他来源的细胞较集中于少数患者。骨、肺候选再依据各自患者内配对效应和HLCA身份构成分别定位。'
edit(3,[32,33],'Lung_04 主要检出于',q10,'QA-04: possible sampling explanation')
q10detail=oldtext(3,list(range(34,40))).replace('\n','').replace('限制了对其他跨患者亚群的覆盖','这种分布可能影响其他跨患者亚群的检出')
assert '这种分布可能' in q10detail
edit(3,list(range(34,40)),'由表 2 与图 7',q10detail,'QA-04: qualify downstream causal wording')
q11='这些核查分别提供训练长度改变后的标签一致性、共同特征集合中的程序基因覆盖，以及正常肺表达背景证据。它们补充支持候选解释，但标签一致性和数值核验不能视为缺失基因及跨平台影响已被排除。<br/>第一，在HLCA健康核心图谱（584,944细胞）的参考空间中，NCI匹配1,631/2,000个模型特征，缺失特征仅为模型结构补零。固定参考参数从200轮训练至500轮后，整体身份大类一致率97.4%，Lung_04、Bone_07、Bone_07_s1、Lung_04_s2分别为96.1%、96.9%、91.4%和100.0%；这项比较检验的是训练长度改变后的标签一致性。<br/>第二，正常表达比较采用378个共同唯一基因。按该特征集合，基质面板7/7、纤维型胶原6/6、FNDC1面板6/6及非COL加工13/13均完整覆盖，其他评估程序覆盖率均≥80%。健康肺参考包含140位供者/17项研究，肿瘤邻近肺为33位/5项研究；47项核验确认了计数、汇总和评分等计算的一致性。<br/>第三，Lung_04相对健康肺泡与外膜参照的基质（+1.704、+1.76）、纤维型胶原（+1.772、+1.54）和交联（+0.610、+0.53）分数均较高；非COL加工相对外膜差值为+0.13、中位百分位为0.759，提示部分加工表达也属于正常外膜成纤维细胞背景。跨数据集结果按分数差和百分位作描述，不使用统计显著性表述，推荐仍以原队列的患者内配对证据为主。'
edit(4,list(range(7,19)),'延长冻结训练',q11,'QA-02: limit conclusions to tested factors; distinguish descriptive elevation')
# Carry forward the user's standing terminology preference without deciding whether Q7 is retained.
q7=oldtext(2,[33,34]).replace('\n','').replace('去除 COL1A1 的敏感性分析','去除COL1A1后重新计算评分')
edit(2,[33,34],'去除 COL1A1',q7,'Standing wording rule only; Q7 retained pending decision')
edit(1,[1],'以下 11 组问答','以下11组问答依据17页主报告及补充材料整理，区分已观察到的结果与可能的原因，逐项回应方法学与生物学问题。','Introduction aligned with qualified explanations',9)
edit(1,[2],'第一部分','第一部分：既有核心审稿问答（Q1–Q5）','Remove obsolete Q4 calculation-detail label',10)
for page,items in plans.items():
 p=doc[page-1]
 for r,txt,size in items:p.add_redact_annot(r,fill=None)
 p.apply_redactions(images=0,graphics=0)
 buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(p.rect.width,p.rect.height))
 for r,txt,size in items:
  fitted=False
  for fs in [size,size-.15,size-.3,size-.45]:
   style=ParagraphStyle('body',fontName='YH',fontSize=fs,leading=fs*1.52,wordWrap='CJK',textColor=colors.HexColor('#243b49'));para=Paragraph(txt,style);_,height=para.wrap(r.width,1000)
   if height<=r.height+1:fitted=True;break
  assert fitted,('overflow',page,r.height,height,txt[:60])
  para.drawOn(c,r.x0,p.rect.height-r.y0-height)
 c.save();overlay=fitz.open(stream=buf.getvalue(),filetype='pdf');p.show_pdf_page(p.rect,overlay,0)
out=W/'expanded_QA_review.pdf';doc.save(out,garbage=4,deflate=True);doc.close();check=fitz.open(out);text=''.join(p.get_text(sort=True) for p in check)
assert len(check)==4 and '敏感性' not in text and '均显著更高' not in text and '均证实上述' not in text
assert '可能的原因是' in text and '一个可能的原因是' in text
assert '6,423' in text and '7/7' in text and '13/13' in text
for i,p in enumerate(check):p.get_pixmap(matrix=fitz.Matrix(1.2,1.2)).save(W/f'page_{i+1}.png')
result=dict(status='STAGED_PASS',source_sha256=hashlib.sha256(backup.read_bytes()).hexdigest(),sha256=hashlib.sha256(out.read_bytes()).hexdigest(),pages=4,questions_retained=11,issues_addressed=['QA-01','QA-02','QA-03','QA-04','QA-05'],question_deletions=0,pending_decision='Whether to retain supplementary-material questions, especially Q7; Q9 denominator precision remains outside this edit scope',edits=records)
(W/'EXPANDED_REVISION_QA.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print('PASS',result['sha256'],len(records))
