import pathlib,csv,json,hashlib
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
P=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase104_full_gene_four_models_20261020')
pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/arphic/uming.ttc'))
pdfmetrics.registerFontFamily('CJK',normal='CJK',bold='CJK',italic='CJK',boldItalic='CJK')
style=ParagraphStyle('body',fontName='CJK',fontSize=8.7,leading=12.7,wordWrap='CJK',textColor=colors.HexColor('#253543'),spaceAfter=5)
small=ParagraphStyle('small',parent=style,fontSize=7.6,leading=10.6)
head=ParagraphStyle('head',parent=style,fontSize=11.5,leading=16,spaceBefore=7,spaceAfter=5,textColor=colors.HexColor('#155C69'))
title=ParagraphStyle('title',parent=style,fontSize=17,leading=23,spaceAfter=7)
def para(t,st=style):return Paragraph(t,st)
with (P/'tables/original14_full_family_results.tsv').open() as f: cand=list(csv.DictReader(f,delimiter='\t'))
with (P/'tables/four_analysis_summary.tsv').open() as f: sums={r['model']:r for r in csv.DictReader(f,delimiter='\t')}
models=['pooled_tumor','source_adjusted_tumor','pooled_normal_min20','source_adjusted_normal_min20']
lookup={(r['gene'],r['model']):r for r in cand}
genes=['NEBL','RTKN','NDUFS7','PNKD','TOR1AIP1','GRID2','TAF10','ZNF526','BEST1','CALHM5','PEF1-AS1','RGS2-AS1','AL096854.1','EBF1']
def result(g,m):
 r=lookup[g,m]
 if r['P']=='NA':return 'NA'
 return ('↑' if float(r['effect'])>0 else '↓')+(' FDR' if float(r['FDR'])<.05 else ' P' if float(r['P'])<.05 else ' NS')
def formatp(v):
 if v=='NA':return 'NA'
 x=float(v);return f'{x:.2g}'
def maketable(rows,widths,fontsize=8):
 t=Table([[para(str(v),ParagraphStyle('cell',parent=small,fontSize=fontsize,leading=fontsize+2)) for v in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAF2F5')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('BOX',(0,0),(-1,-1),.4,colors.HexColor('#B8C7D0')),('INNERGRID',(0,0),(-1,-1),.25,colors.HexColor('#D6DFE5')),('TOPPADDING',(0,0),(-1,-1),2),('BOTTOMPADDING',(0,0),(-1,-1),2)]))
 return t
story=[para('骨肉瘤CAF候选：四项全基因分析简报',title),para('2026-10-20  |  骨16位 / 肺23位  |  6,633个纳入细胞  |  全有效基因分别进行BH校正',small)]
story += [para('主要结果',head),para('联合非配对比较中，原14个候选全部保持原效应方向，9个原始P显著，其中7个通过全基因FDR。来源调整后，12个候选保持同向名义显著，但全检验集合没有FDR显著基因。正常背景扣除后，两类模型分别有6个和4个原候选同向通过全基因FDR。')]
rows=[['分析','有效检验基因','P&lt;0.05（含FDR）','FDR&lt;0.05']]
labels=['A 联合非配对','B 来源调整','C 联合 + 正常背景','D 来源调整 + 正常背景']
for lab,m in zip(labels,models):
 s=sums[m];rows.append([lab,f"{int(s['valid_tests']):,}",f"{int(s['nominal']):,}",f"{int(s['FDR_significant']):,}"])
story.append(maketable(rows,[185,110,130,100]))
story += [para('原14个候选在完整检验家族中的表现',head)]
rows=[['基因','A 联合','B 来源调整','C 正常背景','D 来源+背景']]
for g in genes:rows.append([g]+[result(g,m) for m in models])
story.append(maketable(rows,[125,100,100,100,100],fontsize=7.8))
story.append(para('↑/↓：肺相对骨的效应方向；FDR：全家族FDR&lt;0.05；P：仅原始P&lt;0.05；NS：原始P≥0.05；NA：正常参考不可评估。EBF1原方向为骨高，背景扣除后的肺高不计作原结果支持。',small))
story.append(para('A中通过全基因FDR的7个候选为NEBL、PNKD、TOR1AIP1、GRID2、BEST1、RGS2-AS1、AL096854.1；NDUFS7与TAF10仅名义显著。D中PNKD、TOR1AIP1、TAF10、RGS2-AS1同向通过全基因FDR。',small))
story += [para('方法与解释',head),para('以患者为单位汇总原始计数，计算log2(CP10k+0.1)。按稳定ID或唯一符号对齐，纳入肿瘤三来源共同覆盖且至少1位患者有表达的基因，正常参考缺失不填零。A为Welch非配对比较；B为来源固定效应模型，HC3稳健方差、有效来源自由度19。C/D再减去正常肺与骨的均值差，并合并参考方差。正常骨5位、肺3位，每位≥20个参考细胞。各项在全部有效基因内分别计算BH-FDR；来源模型无可检验变化的基因不计入有效集合；无效检验不补造P值，原14个仅作标注。',small),para('NCI仅有肺样本，来源调整的部位信息来自原两队列，不能称为NCI独立验证。分析对象为原扩展成纤维候选与NCI作者CAF；来源、单核测序及注释差异仍需考虑。正常背景来自外部参考，不是严格配对组织校正。本报告未采用严格子集作为筛选门槛，结论不宣称对细胞身份筛选稳健。',small)]
captions={
'pooled_tumor':'A：骨16位、肺23位。14个原候选中9个原始P显著，7个全基因FDR显著。RGS2-AS1效应2.104，P=1.21e-7，FDR约6.0e-5。',
'source_adjusted_tumor':'B：调整研究来源后，14个候选中12个同向原始P显著；20,161个有效检验中无FDR显著基因。图中没有红点是实际结果。',
'pooled_normal_min20':'C：扣除正常肺减骨背景后，PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1同向FDR显著；EBF1反向FDR显著。',
'source_adjusted_normal_min20':'D：来源调整并扣除正常背景后，PNKD、TOR1AIP1、TAF10、RGS2-AS1同向FDR显著。PEF1-AS1、AL096854.1因参考缺失未绘入C/D。'}
for page_models,page_title in [(models[:2],'全基因火山图：肿瘤部位比较'),(models[2:],'全基因火山图：正常背景敏感性')]:
 story += [PageBreak(),para(page_title,title),para('红：FDR&lt;0.05；橙：仅P&lt;0.05；灰：不显著。深色圈及基因名标出原14个候选。横轴为肺减骨的log2(CP10k+0.1)效应，纵轴为原始P的-log10；虚线为P=0.05，未设置效应量筛选门槛。',small)]
 for m in page_models:
  story += [Image(str(P/'figures'/(m+'_volcano.png')),width=525,height=267.3),para(captions[m],small),Spacer(1,9)]
def footer(c,doc):
 c.saveState();c.setStrokeColor(colors.HexColor('#CCD7DE'));c.line(35,29,A4[0]-35,29)
 c.setFont('CJK',7);c.setFillColor(colors.HexColor('#61727E'));c.drawString(35,17,'数据版本：phase104 / 2026-10-20  ·  全基因FDR，不沿用14候选内部FDR')
 c.drawRightString(A4[0]-35,17,str(doc.page));c.restoreState()
out=P/'report/CAF_four_fullgene_analyses_brief_20261020.pdf'
doc=SimpleDocTemplate(str(out),pagesize=A4,rightMargin=35,leftMargin=35,topMargin=30,bottomMargin=40,title='骨肉瘤CAF候选四项全基因分析简报',author='CAF项目分析',pageCompression=1)
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print('PDF_COMPLETE',out,hashlib.sha256(out.read_bytes()).hexdigest())
