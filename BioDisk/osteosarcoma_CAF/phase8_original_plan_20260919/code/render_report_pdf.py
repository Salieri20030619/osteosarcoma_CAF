from pathlib import Path
import re,html,io,json
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4,A3,landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfgen import canvas
from PIL import Image
from pypdf import PdfReader,PdfWriter
import pypdfium2 as pdfium
P=Path(r'C:\Users\31937\caf_original_plan_20260919');src=P/'delivery/original_plan_report.md'
out=P/'delivery/骨肉瘤CAF_原计划首版报告_20260919.pdf';qa=P/'pdf_qa';qa.mkdir(exist_ok=True)
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
styles={'body':ParagraphStyle('body',fontName='STSong-Light',fontSize=10,leading=15,spaceAfter=7,wordWrap='CJK'),
 'title':ParagraphStyle('title',fontName='STSong-Light',fontSize=19,leading=27,spaceAfter=15,wordWrap='CJK',keepWithNext=True),
 'h2':ParagraphStyle('h2',fontName='STSong-Light',fontSize=14,leading=20,spaceBefore=12,spaceAfter=8,keepWithNext=True,wordWrap='CJK'),
 'h3':ParagraphStyle('h3',fontName='STSong-Light',fontSize=11.5,leading=17,spaceBefore=10,spaceAfter=6,keepWithNext=True,wordWrap='CJK'),
 'cell':ParagraphStyle('cell',fontName='STSong-Light',fontSize=8,leading=11,wordWrap='CJK')}
names={'sample':'样本','dataset':'队列','site':'部位','strict_CAF':'严格CAF数','subcluster':'子簇','cells':'细胞数','patients_with_ge20_cells':'达到20细胞的患者数','patients_with_ge10_cells':'达到10细胞的患者数','state':'状态','gene':'基因','patients':'配对患者数','median_log2CPM_difference':'中位log2CPM差','positive_patients':'正向患者','negative_patients':'负向患者','program':'表达程序','OR_per_SD':'每SD的OR','CI_low':'95%CI下限','CI_high':'95%CI上限','P_value':'P值','FDR_BH':'BH FDR','two_gene_positive_cells':'至少2基因阳性数','two_gene_positive_pct':'阳性比例%','positive_genes':'阳性组合','donors':'供者数','median_RNA_target_coverage':'RNA目标覆盖中位数','median_RNA_known_CAF_fraction_gate':'RNA门内已知CAF占比中位数','worst_RNA_CNV_high_fraction_gate':'RNA门内CNV高细胞最大占比','Frozen_H0_top30':'H0冻结签名','Frozen_H2_top30':'H2冻结签名','Contractile':'收缩程序','Inflammatory':'炎症程序','Antigen_presentation':'抗原呈递程序','Matrix_remodeling':'基质重塑程序','Osteogenic_full_customer':'成骨完整面板','Calcium_regulation_full_customer':'钙调节完整面板','Gong_human_PTGS2_coexpression_reference':'Gong人肺共表达参考','Alveolar_reference':'肺泡参考','Adventitial_reference':'外膜参考'}
def clean(t):return names.get(t,t)
def fmt(t):
 t=html.escape(t)
 t=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#225277">\1</link>',t)
 return t
def para(t,sty='body'):return Paragraph(fmt(t),styles[sty])
lines=src.read_text(encoding='utf-8').splitlines();story=[];images=[];i=0
while i<len(lines):
 l=lines[i].strip()
 if not l:i+=1;continue
 if l.startswith('|'):
  rows=[]
  while i<len(lines) and lines[i].strip().startswith('|'):
   vals=[v.strip() for v in lines[i].strip().strip('|').split('|')]
   if not all(re.fullmatch(r'[-: ]+',v) for v in vals):rows.append([para(clean(v),'cell') for v in vals])
   i+=1
  n=len(rows[0]);width=487
  if n>=5:widths=[width*.25]+[width*.75/(n-1)]*(n-1)
  elif n==4:widths=[width*.15,width*.22,width*.22,width*.41] if len(rows)>12 else [width/n]*n
  else:widths=[width/n]*n
  table=Table(rows,colWidths=widths,repeatRows=1,hAlign='LEFT')
  table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eef2')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#728191')),('LINEBELOW',(0,1),(-1,-1),.2,colors.HexColor('#d6dce1')),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
  story += [table,Spacer(1,9)];continue
 m=re.match(r'!\[([^\]]*)\]\(([^)]+)\)',l)
 if m:
  images.append((m.group(1),(src.parent/m.group(2)).resolve()));story.append(para('图件 '+m.group(1)+'，见文末全尺寸附图 A'+str(len(images))+'。'));i+=1;continue
 if l.startswith('# '):story.append(para(l[2:],'title'))
 elif l.startswith('### '):story.append(para(l[4:],'h3'))
 elif l.startswith('## '):story.append(para(l[3:],'h2'))
 else:story.append(para(l))
 i+=1
body=qa/'body.pdf'
def footer(c,d):
 c.setFont('STSong-Light',8);c.setFillColor(colors.HexColor('#657180'));c.drawRightString(A4[0]-54,28,str(d.page))
doc=SimpleDocTemplate(str(body),pagesize=A4,rightMargin=54,leftMargin=54,topMargin=48,bottomMargin=45,title='骨肉瘤肺转移相关CAF原计划首版分析报告')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
appendix=qa/'figures.pdf';c=canvas.Canvas(str(appendix),pagesize=landscape(A3));w,h=landscape(A3)
for j,(title,path) in enumerate(images,1):
 c.setFont('STSong-Light',15);c.drawString(35,h-35,'附图 A'+str(j)+' '+title)
 im=Image.open(path);iw,ih=im.size;scale=min((w-70)/iw,(h-100)/ih)
 c.drawImage(str(path),(w-iw*scale)/2,(h-ih*scale)/2-10,iw*scale,ih*scale,preserveAspectRatio=True)
 c.showPage()
c.save();writer=PdfWriter()
for f in [body,appendix]:writer.append(PdfReader(str(f)))
with open(out,'wb') as f:writer.write(f)
d=pdfium.PdfDocument(str(out))
for j in range(len(d)):d[j].render(scale=1.3).to_pil().save(qa/f'page-{j+1}.png')
print('PDF_CREATED',out,'PAGES',len(d))
