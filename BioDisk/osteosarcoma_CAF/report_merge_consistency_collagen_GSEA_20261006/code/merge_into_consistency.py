from pathlib import Path
from io import BytesIO
from copy import deepcopy
import re,json,hashlib
from docx import Document
from docx.shared import RGBColor,Pt
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.enum.text import WD_BREAK

B=Path(__file__).parent;O=B/'consistency_merge';O.mkdir(exist_ok=True)
A=B.parent/'骨肉瘤CAF功能分类一致性分析报告_20261005_图文增强版_P值修订版.docx'
S=B.parent/'骨肉瘤胶原与基质相关基因上调分析报告_20261006_全亚群与GSEA整合版.docx'
DEST=B.parent/'骨肉瘤CAF功能分类一致性分析报告_20261006_胶原基质与GSEA合并版.docx'
d=Document(A);s=Document(S);ap=list(d.paragraphs);sp=list(s.paragraphs)
def text_replace(p,text):
 p.clear();p.add_run(text)
def pagebreak():
 p=OxmlElement('w:p');r=OxmlElement('w:r');br=OxmlElement('w:br');br.set(qn('w:type'),'page');r.append(br);p.append(r);return p
text_replace(ap[0],'骨肉瘤CAF功能分类与胶原基质联合分析报告')
text_replace(ap[1],'分类一致性  细胞组成  全亚群胶原基质证据  GSEA')
text_replace(ap[2],'2026年10月6日  分类一致性报告与胶原基质证据整合版')
text_replace(ap[5],'功能解释置信度为高11群、中12群、低10群，支持保留原编号并附加可重叠的功能注释。胶原基质分析进一步得到7个亚群、39个基因的54条全来源名义上调记录；F2-C3与F2-C1分别构成纤维胶原及IV型胶原交联相关的并行候选。GSEA与患者级效应按各自比较背景联合解释。[1–5]')
text_replace(ap[8],'组成按总体、来源、患者记录和F类内部四种分母计算。功能评分及GSEA沿用主要来源内同F类对照，使用全转录组排序与1000次基因集置换。胶原基因补充同患者同部位配对检验。原始P<0.05用于名义筛选，FDR保留原同类家族；模型P与患者级P分开记录。[3–5]')
# Split the base at its references and appendices; preserve every result and original table.
body=d.element.body
ref_start=ap[79]._p;appendix_start=ap[84]._p
children=list(body);i_ref=children.index(ref_start);i_app=children.index(appendix_start)
refs=[deepcopy(e) for e in children[i_ref:i_app] if e.tag!=qn('w:sectPr')]
base_appendices=[deepcopy(e) for e in children[i_app:] if e.tag!=qn('w:sectPr')]
for e in children[i_ref:]:
 if e.tag!=qn('w:sectPr'):body.remove(e)
# The repair paragraph is consolidated with the later GSEA version description.
for p in [ap[65],ap[66]]:body.remove(p._p)

# Reuse base layout and styles. Matching layout was verified before merge.
for sn in ['Normal','Title','Subtitle','Heading 1','Heading 2','Caption']:
 print('STYLE',sn,'base',d.styles[sn].font.name,d.styles[sn].font.size,'source',s.styles[sn].font.name,s.styles[sn].font.size)
print('NUMBERED_PARAGRAPHS',len(s.element.body.findall('.//'+qn('w:numPr'))))
assert len(s.element.body.findall('.//'+qn('w:numPr')))==0
assert len(s.sections)==len(d.sections)==1
assert s.sections[0].page_width==d.sections[0].page_width and s.sections[0].page_height==d.sections[0].page_height
def append(el):
 sect=body.find(qn('w:sectPr'));sect.addprevious(el)
relmap={}
for r in s.part.rels.values():
 if r.reltype.endswith('/image'):
  rid,img=d.part.get_or_add_image(BytesIO(r.target_part.blob));relmap[r.rId]=rid
 elif r.is_external:relmap[r.rId]=d.part.relate_to(r.target_ref,r.reltype,is_external=True)
def import_element(el):
 e=deepcopy(el)
 for node in e.iter():
  for attr in [qn('r:embed'),qn('r:link'),qn('r:id')]:
   if attr in node.attrib:
    old=node.attrib[attr]
    if old not in relmap:raise ValueError('Unmapped relationship '+old)
    node.attrib[attr]=relmap[old]
 return e

# Recast the added report as a continuation, then renumber only its references.
text_replace(sp[0],'胶原基质相关基因与GSEA联合证据');sp[0].style=s.styles['Heading 1']
text_replace(sp[3],'胶原与基质分析概览')
text_replace(sp[8],'胶原基因的患者配对分析方法')
text_replace(sp[56],'功能比较覆盖33个亚群的全范围结果及20个亚群的Fibro限定结果。分类、坐标及历史GSEA沿用8,174个冻结细胞；双联体证据修复后标记26个模型阳性，排除视图保留8,148个细胞。主要候选机制结果在修复核对中保持一致；失去有效重复的辅助比较停止患者级推断。历史GSEA未全部重算，目标及参考供者数保存在配套统计表，不作为有效配对N。[4,5]')
text_replace(sp[79],'分类与胶原基质证据共同支持多亚群推进')
text_replace(sp[80],'分类一致性分析支持保留F1至F5及33个内部编号，以允许重叠的功能注释解释现有结构。结合基因与GSEA，F2-C3和F2-C1分别作为纤维胶原及IV型胶原交联相关候选；F2-C0、F4-C9补充组分线索，F1-C10及F4-C8保留单患者功能候选。')
text_replace(sp[86],'整合结果与资料来源')
text_replace(sp[87],'本报告保留分类一致性分析的9组图及胶原基质分析的7组图，共16组。附表一与二保留33亚群的工作名称、组成和置信度，附表三列出54条全来源名义上调记录。完整GSEA及患者级统计保留于配套结果。')
# The new references replace the added document's repeated provenance paragraph.
text_replace(sp[88],'GSEA功能富集与患者级效应共同用于解释亚群差异。保留未通过校正、来源相反及Fibro限定后减弱的结果；所有分析沿用原检验家族，不因整合报告重新筛选或改变统计值。')
text_replace(sp[90],'附表三 全来源54条名义上调记录')
text_replace(sp[93],'附表三 全来源名义上调记录续')
text_replace(sp[95],'附表三 全来源名义上调记录续')
for p in s.paragraphs:
 if p.text:
  t=re.sub(r'图\s*(\d+)',lambda m:'图'+str(int(m.group(1))+9),p.text)
  t=re.sub(r'(?<!附)表\s*(\d+)',lambda m:'表'+str(int(m.group(1))+2),t)
  if t!=p.text:text_replace(p,t)

# Split supplemental appendices from its narrative and insert base appendices first.
schildren=list(s.element.body);s_app=schildren.index(sp[89]._p)
skip={sp[1]._p,sp[2]._p}
append(pagebreak())
for e in schildren[:s_app]:
 if e.tag==qn('w:sectPr') or e in skip:continue
 append(import_element(e))
 # Add a caption to the priority table, which had no caption in the source.
 if e is sp[80]._p:
  p=deepcopy(sp[7]._p)
  for r in list(p):
   if r.tag!=qn('w:pPr'):p.remove(r)
  r=OxmlElement('w:r');t=OxmlElement('w:t');t.text='表9  分类与胶原基质联合候选优先级';r.append(t);p.append(r);append(p)
# Original source list retained and extended with the exact incorporated report.
for e in refs:
 if e.tag==qn('w:p'):
  text=''.join(e.itertext())
  if '[1]' not in text and '[2]' not in text and '[3]' not in text and '[4]' not in text:continue
 append(e)
p=deepcopy(ap[83]._p)
for e in list(p):
 if e.tag!=qn('w:pPr'):p.remove(e)
r=OxmlElement('w:r');t=OxmlElement('w:t');t.text='[5] 骨肉瘤胶原与基质相关基因上调分析报告 全亚群与GSEA整合版及配套统计，2026年10月6日。';r.append(t);p.append(r);append(p)
for e in base_appendices:append(e)
for e in schildren[s_app:]:
 if e.tag!=qn('w:sectPr'):append(import_element(e))

# Use paragraph page breaks so an otherwise full page cannot spill a blank break page.
for el in list(body):
 if el.tag!=qn('w:p'):continue
 if el.findall('.//'+qn('w:br')) and not el.findall('.//'+qn('w:t')) and not el.findall('.//'+qn('w:drawing')):
  nxt=el.getnext()
  while nxt is not None and nxt.tag==qn('w:p') and not nxt.findall('.//'+qn('w:t')) and not nxt.findall('.//'+qn('w:drawing')):nxt=nxt.getnext()
  if nxt is not None and nxt.tag==qn('w:p'):
   pr=nxt.find(qn('w:pPr'))
   if pr is None:pr=OxmlElement('w:pPr');nxt.insert(0,pr)
   if pr.find(qn('w:pageBreakBefore')) is None:pr.append(OxmlElement('w:pageBreakBefore'))
   body.remove(el)
# Match compact original appendix intent while keeping every label and value.
for table in d.tables:
 if len(table.columns)==7 and table.rows[0].cells[1].text=='细胞数':
  for row in table.rows:
   trpr=row._tr.get_or_add_trPr()
   for el in list(trpr):
    if el.tag==qn('w:trHeight'):trpr.remove(el)
   for cell in row.cells:
    cp=cell._tc.get_or_add_tcPr();mar=cp.find(qn('w:tcMar'))
    if mar is None:mar=OxmlElement('w:tcMar');cp.append(mar)
    for edge in ['top','bottom','left','right']:
     el=mar.find(qn('w:'+edge))
     if el is None:el=OxmlElement('w:'+edge);mar.append(el)
     el.set(qn('w:w'),'40' if edge in ['top','bottom'] else '60');el.set(qn('w:type'),'dxa')
    for p in cell.paragraphs:
     p.paragraph_format.space_before=Pt(0);p.paragraph_format.space_after=Pt(0);p.paragraph_format.line_spacing=1.0
     for run in p.runs:run.font.size=Pt(8.5)

# Ensure unique drawing and bookmark IDs after merging binary image relationships.
for i,el in enumerate(body.findall('.//'+qn('wp:docPr')),1):el.set('id',str(i))
for i,el in enumerate(body.findall('.//'+qn('pic:cNvPr')),1):el.set('id',str(i))
for p in d.paragraphs:
 if p.style.name in ['Title','Subtitle','Heading 1','Heading 2']:
  for run in p.runs:run.font.color.rgb=RGBColor(0,0,0)
assert len(d.inline_shapes)==16 and len(d.tables)==14
captions=[p.text for p in d.paragraphs if re.match(r'^图\d+\s',p.text)]
assert [int(re.match(r'^图(\d+)',x).group(1)) for x in captions]==list(range(1,17))
# Verify every original image's bytes survive, independent of renamed media parts.
hashes=lambda doc:[hashlib.sha256(rel.target_part.blob).hexdigest() for rel in doc.part.rels.values() if rel.reltype.endswith('/image')]
assert set(hashes(Document(A))+hashes(Document(S)))<=set(hashes(d))
d.save(DEST)
record={'report':str(DEST),'base_report':str(A),'added_report':str(S),'base_sha256':hashlib.sha256(A.read_bytes()).hexdigest(),'added_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'report_sha256':hashlib.sha256(DEST.read_bytes()).hexdigest(),'figures':16,'tables':14,'all_original_images_preserved':True,'nominal_records':54,'subgroups':33,'source_reports_preserved':True,'visual_review_pending':True}
(O/'MERGE_ACCEPTANCE.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(record,ensure_ascii=False))
