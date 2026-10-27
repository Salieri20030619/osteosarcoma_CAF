import pathlib,hashlib,json,subprocess,html
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
P=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase110_flow_marker_append_20261020')
source=P/'resources/original_user_report.pdf'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='c90188b91fe40ca47fcfe82a22c95011b56302e97d6b9ee1f621905d240a4ac0'
pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/arphic/uming.ttc'))
body=ParagraphStyle('body',fontName='CJK',fontSize=10,leading=14.5,wordWrap='CJK',spaceAfter=7,textColor=colors.HexColor('#243744'))
small=ParagraphStyle('small',parent=body,fontSize=8.5,leading=11.5,spaceAfter=4)
head=ParagraphStyle('head',parent=body,fontSize=12,leading=17,spaceBefore=7,spaceAfter=6,textColor=colors.HexColor('#155C69'))
title=ParagraphStyle('title',parent=body,fontSize=17,leading=24,spaceAfter=10,textColor=colors.black)
def pa(t,s=body):return Paragraph(t,s)
def table(rows,widths):
 st=ParagraphStyle('cell',parent=small,fontSize=9,leading=12.5,spaceAfter=0)
 t=Table([[pa(str(v),st) for v in r] for r in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAF2F5')),('GRID',(0,0),(-1,-1),.35,colors.HexColor('#D1DCE2')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]))
 return t
story=[pa('附录 六个候选的流式marker建议',title),
pa('建议优先对GRID2开展活细胞表面染色预实验；BEST1作为第二顺位的定位与方法探索。PNKD、TOR1AIP1、TAF10不列入常规活细胞表面分选panel；RGS2-AS1改用RNA检测。当前证据不足以将六者中任何一个直接定义为已验证的人骨肉瘤CAF分选marker。'),
pa('六个候选的用途分级',head)]
rows=[['候选','活细胞表面分选建议','更合适的检测路线及依据'],
['GRID2','优先级1 条件性探索','胞外N端抗体有活细胞成像证据；先验证人样本表面信号和流式区分度，再考虑分选。尚无本项目CAF流式验证。[11]'],
['BEST1','优先级2 暂不作主分选标志','膜相关通道不等于存在可用的胞外流式表位。先做定位检测；固定通透流式也须单独验证克隆。[12,13]'],
['TOR1AIP1','不推荐常规活细胞表面染色','核内膜LAP1。优先固定通透后的定位或胞内检测；核膜定位不代表细胞表面可及。[14]'],
['TAF10','不推荐常规活细胞表面染色','核内转录相关蛋白。考虑固定通透表达验证，不能仅凭阳性定义CAF身份或特定亚群。[15]'],
['PNKD','不推荐作为当前表面分选标志','亚型和定位需核对，膜相关注释不证明抗体能从细胞外结合。优先亚型特异RNA/蛋白验证；不能误用MR1抗体。'],
['RGS2-AS1','不适用抗体表面分选','非编码RNA，无相应蛋白表位。可考虑转录本特异RNA流式或原位检测；常规RNA流式需固定通透，不能据此获得可继续培养的活细胞。[16]']]
story += [table(rows,[70,148,307]),pa('优先级表示验证顺序，不代表统计显著性排序、CAF特异性或现成可用的分选性能。前文14→8→6及全部P、FDR结果保持不变。',small),
pa('GRID2为何可以先试',head),
pa('可核对的试剂实例为Alomone AGC-039，Thermo Fisher目录号AGC-039-200UL。产品页标注人、鼠、大鼠反应性，免疫原位于大鼠GluD2胞外N端206-218位；展示的是完整活大鼠PC12细胞的表面成像。该页列出的应用为WB和ICC/IF，没有将流式列为已验证应用。因此它是试剂评估线索，不是已通过人CAF流式验证的采购结论。[11]',small),
pa('若用于活细胞分选后培养，应核对具体配方：上述规格列有0.05%叠氮钠，应优先向厂家确认适合活细胞用途的无叠氮钠规格及相容性。荧光偶联版本或“胞外”名称本身不能替代人细胞流式验证。',small),
pa('BEST1需要先解决什么',head),
pa('Abcam ab259836（EPR23926-112）当前产品资料将Flow Cyt标为不推荐；其固定通透ICC证据不能转述为活细胞流式证据。ab14927当前页面列出ICC、ICC/IF、IP、WB，未列流式；其免疫原在人体BEST1的1-100位范围，具体表位未公开。这两款均不应作为已经验证的活细胞BEST1分选抗体推荐。上述判断针对具体产品，不等于BEST1永远无法做流式。[12,13]',small)]
story += [PageBreak(),pa('建议的验证路径与panel安排',title),
pa('先建立候选成纤维细胞背景',head),
pa('活细胞预实验可先用活性染料、单细胞门及CD45/CD31排除明显免疫和内皮群；如有红细胞或上皮混入，再按样本情况评估CD235a或EpCAM。FAP、CD140a（PDGFRα）、CD90可作为候选富集轴，先看分布，不预设必须三阳，也不把任一单标志等同于CAF。抗体克隆、通道和染色条件需按仪器与样本验证后确定。[17]'),
pa('骨肉瘤肿瘤细胞本身可具有间质表达，CD45-/CD31-/EpCAM-或FAP阳性不能自动排除恶性细胞。GRID2高低群应在同一候选背景中比较；分选后结合成纤维标志表达、形态/组织定位及可获得的肿瘤遗传证据复核身份。单一CNV算法的二倍体预测也不能单独确证非恶性CAF。',small),
pa('把表面可及性和抗体特异性分开验证',head),
pa('同一样本平行比较不通透表面染色与固定通透总蛋白信号；仅通透后出现信号的组合不进入活细胞表面分选。检查组织消化对胞外表位的影响，并使用确认表达的人阳性/阴性材料或敲低/敲除对照评估特异性。未染、单染补偿及FMO用于背景、补偿和设门；同型或二抗对照不能替代靶标特异性验证。[18]'),
pa('进入分选的条件应是：在人样本中有可重复的非通透特异信号、可解释的高低群分离，并通过分选后身份和目标表达复核。首次分选保留GRID2高、低/阴性两组作对照，不预设GRID2阳性就是CAF。若GRID2未通过验证，先用已验证的候选富集组合获得细胞，再用RNA或蛋白读数研究这6个基因，不勉强将六者都放入表面panel。',small),
pa('固定样本分析与活细胞分选分别设计',head),
pa('BEST1、TOR1AIP1、TAF10、PNKD可按克隆和固定通透相容性评估胞内流式，但本次不承诺现有抗体可直接用于该应用。RGS2-AS1需核对ENSG00000285280及旧名AL390957.1/LINC01032，确认探针对应转录本；RNA流式同时设置阳性、阴性和RNA质量对照。固定通透方案用于分析，不与活细胞培养分选的目标混用。[16,18]',small),
pa('新增参考与试剂资料',head)]
refs=[
(11,'Thermo Fisher / Alomone GRID2 extracellular AGC-039-200UL 产品页','https://www.thermofisher.com/antibody/product/GRID2-extracellular-Antibody-Polyclonal/AGC-039-200UL'),
(12,'Abcam BEST1 EPR23926-112 ab259836 产品应用与验证资料','https://www.abcam.com/en-us/products/primary-antibodies/bestrophin-best1-antibody-epr23926-112-ab259836'),
(13,'Abcam BEST1 ab14927 产品应用与免疫原资料','https://www.abcam.com/en-us/products/primary-antibodies/bestrophin-best1-antibody-ab14927'),
(14,'Human Protein Atlas TOR1AIP1 亚细胞定位','https://www.proteinatlas.org/ENSG00000143337-TOR1AIP1/subcellular'),
(15,'Human Protein Atlas TAF10 核内定位与功能注释','https://www.proteinatlas.org/ENSG00000166337-TAF10'),
(16,'Thermo Fisher PrimeFlow RNA Assay 用户指南 MAN0019788','https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0019788_PrimeFlowRNAAssay_UG.pdf'),
(17,'Sahai E et al. A framework for advancing our understanding of cancer-associated fibroblasts. Nature Reviews Cancer. 2020.','https://www.nature.com/articles/s41568-019-0238-1'),
(18,'Abcam Flow cytometry protocol for intracellular and extracellular targets','https://www.abcam.com/en-us/technical-resources/protocols/flow-cytometry-for-intracellular-and-extracellular-targets')]
for n,label,url in refs:story.append(pa(f'[{n}] <link href="{url}" color="#155C69">'+html.escape(label)+'</link>',small))
def footer(c,doc):
 c.saveState();c.setStrokeColor(colors.HexColor('#CFD8DE'));c.line(35,28,A4[0]-35,28);c.setFont('CJK',7);c.setFillColor(colors.HexColor('#647681'))
 c.drawString(35,16,'追加附录  流式marker用途与验证建议  |  原统计结果及前8页保持不变')
 c.drawRightString(A4[0]-35,16,str(8+doc.page));c.restoreState()
appendix=P/'report/flow_marker_appendix.pdf'
SimpleDocTemplate(str(appendix),pagesize=A4,leftMargin=35,rightMargin=35,topMargin=30,bottomMargin=40,title='六个候选的流式marker建议').build(story,onFirstPage=footer,onLaterPages=footer)
out=P/'report/user_report_with_flow_markers.pdf'
subprocess.run(['pdfunite',str(source),str(appendix),str(out)],check=True)
(P/'resources/SOURCES_AND_RECOMMENDATIONS.json').write_text(json.dumps(dict(references=[dict(number=n,title=l,url=u) for n,l,u in refs],recommendation='GRID2 first conditional surface pilot; BEST1 second exploratory; no validated live human osteosarcoma CAF marker among six',GRID2_flow_validated=False,GRID2_live_imaging_species='rat PC12',BEST1_ab259836_flow='manufacturer not recommended',original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),merged_sha256=hashlib.sha256(out.read_bytes()).hexdigest()),ensure_ascii=False,indent=2),encoding='utf-8')
print('MERGED',hashlib.sha256(out.read_bytes()).hexdigest())
print(subprocess.run(['pdfinfo',str(out)],capture_output=True,text=True).stdout)
