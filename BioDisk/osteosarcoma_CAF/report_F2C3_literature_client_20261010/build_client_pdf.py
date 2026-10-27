from pathlib import Path
import json,html,hashlib
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
import pymupdf as fitz
B=Path('C:/Users/31937/CAF_F2C3_client_literature_20261010');O=B/'output/pdf'
refs=json.loads(Path('C:/Users/31937/CAF_F2C3_subtype_review_20261010/inputs/verified_reference_registry.json').read_text(encoding='utf-8'));refs.append({'number':11,'citation':'Liu et al. Front Oncol. 2021;11:709210.','title':'Single-Cell Transcriptomics Reveals the Complexity of the Tumor Microenvironment of Treatment-Naive Osteosarcoma.','doi':'10.3389/fonc.2021.709210','url':'https://doi.org/10.3389/fonc.2021.709210','origin':'补充原始研究'})
pdfmetrics.registerFont(TTFont('YH','C:/Windows/Fonts/msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('YHB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0));pdfmetrics.registerFontFamily('YH',normal='YH',bold='YHB',italic='YH',boldItalic='YHB');C=colors.HexColor;W,H=A4;FW=W-88
styles={k:ParagraphStyle(k,fontName='YHB' if k in ['title','h2','h3'] else 'YH',fontSize=fs,leading=ld,spaceAfter=gap,wordWrap='CJK',textColor=C('#17485e') if k in ['title','h2','h3'] else C('#263b49')) for k,fs,ld,gap in [('title',21,29,15),('h2',15,22,12),('h3',11.5,18,8),('body',10.3,17.8,10),('small',8.8,13.6,7),('cell',9.1,14.6,0),('ref',8.5,12.5,9)]};story=[]
def P(s,k='body'):return Paragraph(s,styles[k])
def add(s,k='body'):story.append(P(s,k))
def table(head,rows,widths):
 t=Table([[P('<b>'+x+'</b>','cell') for x in head]]+[[P(str(x),'cell') for x in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),C('#e9f2f5')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),9),('RIGHTPADDING',(0,0),(-1,-1),9),('TOPPADDING',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),9),('LINEBELOW',(0,0),(-1,0),.7,C('#adc8d4')),('LINEBELOW',(0,1),(-1,-1),.4,C('#d9e5eb'))]));story.extend([t,Spacer(1,13)])
def box(s):
 t=Table([[P(s)]],colWidths=[FW]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),C('#edf5f7')),('BOX',(0,0),(-1,-1),.6,C('#bfd6e0')),('LEFTPADDING',(0,0),(-1,-1),14),('RIGHTPADDING',(0,0),(-1,-1),14),('TOPPADDING',(0,0),(-1,-1),12),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,14)])
def newpage():story.append(PageBreak())
# 1
add('骨肉瘤F2-C3亚群\n文献交叉对比报告'.replace('\n','<br/>'),'title');add('对照基准：2026年10月6日报告  |  文献核查日期：2026年10月10日','small')
add('核心问题：F2-C3此前是否已在骨肉瘤中被研究？','h2')
box('<b>骨肉瘤中已有明确的胶原／基质相关亚群报道。</b>与F2-C3最接近的既往先例是COL5A1+间质细胞和matrix CAF；FNDC1+ CAF等也属于需要比较的相关群体。<br/><b>本项目适合强调该方向的进一步细化描述，目前不宜以“首次发现独立新亚群”作为核心创新主张。</b>')
add('本次对比使用的F2-C3定义','h3')
table(['原报告层次','用于文献对照的特征'],[
['亚群名称','F2-C3为本项目F2内部的工作聚类编号，不是跨研究通用的细胞类型名称。'],
['表型特征','CD248、POSTN、VCAN及胶原相关表达。'],
['胶原及加工相关组合','COL3A1、COL5A1／COL5A2、PCOLCE、LOXL2等，体现原报告关注的III／V型胶原及加工相关表达。']], [116,FW-116])
add('交叉对比范围','h3')
add('对所提供清单中的6篇文章逐一核查，并补充5篇相关原始研究或图谱，共11篇。逐篇比较文献所定义的细胞群、展示标志及功能描述，再与10月6日报告的F2-C3定义对应。')
add('主要判断','h3')
add('已有文献覆盖了F2-C3的主要胶原／基质背景；原报告较具体的工作内容，是描述这一背景下的标志组合与胶原加工相关表达侧重。<b>组合侧重可以作为细化研究的内容，但不等于已经建立新的亚群身份。</b>')
add('以下各页以文献交叉对照为主体。文献中的原有亚群名称均按作者用法保留；基因名称相同或不同，均不单独用于裁定两个亚群是否完全等同。','small')
# 2
newpage();add('一、最直接的三项既往亚群对照','h2')
table(['文献及已报道亚群','与10月6日F2-C3的共同点与差异','交叉对比判断'],[
['<b>Zheng等，2024 [2]</b><br/>J Cancer Res Clin Oncol<br/><br/><b>COL5A1+ MSC</b><br/>展示COL5A1、COL16A1、EMILIN1、ADAMTS2等基质重塑相关基因。',
'<b>共同点：</b>F2-C3原报告包含COL5A1、COL16A1等，胶原及ECM表达方向与该群明确相交。<br/><br/><b>侧重点：</b>本报告以CD248／POSTN／VCAN表型及III／V型胶原、PCOLCE／LOXL2组合进一步描述C3。',
'<b>最直接的胶原型先例。</b><br/><br/>“发现骨肉瘤胶原高表达细胞群”本身已有先例；C3的表述应突出组合特征的细化。'],
['<b>Zheng等，2024 [1]</b><br/>J Hematol Oncol<br/><br/><b>matrix CAF（mCAF）</b><br/>补充图展示FAP、MMP11、COL1A1、COL11A1、POSTN、FN1。',
'<b>共同点：</b>POSTN、FN1及胶原／ECM相关表达与C3相交。<br/><br/><b>侧重点：</b>文献建立的是基质型CAF框架；10月6日报告侧重描述C3的胶原组成和加工相关成员。',
'<b>最直接的既有基质型框架。</b><br/><br/>将C3定位为相关框架中的候选细分状态较为合适，名称不同不构成独立创新证据。'],
['<b>Zeng等，2026卷期 [3]</b><br/>Molecular Cancer<br/><br/><b>FNDC1+ CAF</b><br/>展示FNDC1, COMP, MMP2,<br/>COL11A1, COL10A1, P4HA3。<br/>涉及ECM合成与代谢。',
'<b>共同点：</b>同属骨肉瘤基质相关细胞研究，功能背景存在重叠。<br/><br/><b>侧重点：</b>其展示组合偏FNDC1／COMP／COL10A1等；C3原报告突出III／V型胶原及PCOLCE／LOXL2。',
'<b>相关方向，展示组合不同。</b><br/><br/>应保留为近邻对照，现有资料不足以将C3直接等同于该群。']], [151,215,FW-366])
add('如何理解这些差异','h3')
add('表中的“侧重点不同”是原报告与文献展示内容的比较，不表示文献群不表达PCOLCE／LOXL2，也不表示C3不表达文献中的其他基因。完整表达状态可能重叠，不能用一份短marker表的缺项证明新亚群。')
add('[3]发表于2026年卷期，在线发表时间为2025年12月，属于10月6日报告之前的研究。MSC为间质细胞，ECM为细胞外基质。','small')
#3
newpage();add('二、其他骨肉瘤亚群研究的交叉对照','h2')
table(['文献与细胞群','与F2-C3的对应关系','对新颖性判断的意义'],[
['<b>He等，2024 [4]</b><br/>J Bone Oncology<br/><br/>SPP1、IGFBP4、COL10A成纤维细胞群。',
'该研究采用另一套CAF分群命名；C3以CD248／POSTN／VCAN及胶原加工相关组合描述。SPP1、IGFBP4和COL10A三个群均应纳入对照。',
'已存在CAF内部精细分群先例。不能只比较其中的SPP1或COL10A群，就排除与全部既有状态的联系。'],
['<b>Liu等，2024 [5]</b><br/>BMC Medicine<br/><br/>myCAF、OstCAF、apCAF。',
'文献分别展示POSTN／TAGLN、ALPL／RUNX2和抗原呈递相关特征。C3与其中部分基质表达相交；研究组织环境和定义方式不同。',
'属于相关的既有CAF异质性框架。共享POSTN或采用不同名称，都不能直接决定是否为同一亚群。'],
['<b>Liu等，2021 [11]</b><br/>Front Oncol<br/><br/>C1_CAFs、C2_CAFs、C3_CAFs。',
'该文C1偏MMP9／MCAM及血管相关特征，C2偏OMD／OGN及骨化，C3以TOP2A／MKI67增殖特征定义。其C3与本项目C3的定义不同。',
'<b>相同“C3”编号不代表重复。</b>应比较具体基因组合与状态，不能按聚类数字认定同群。'],
['<b>Budhathoki等 [6]</b><br/>2026年<br/>OsteoCAR图谱预印本v1<br/><br/>含mCAF相关状态。',
'作者公开代码列出SPARC、THBS2、FN1、POSTN相关签名，与C3部分特征相交；图谱包含本项目主要数据来源SCPCP000017。',
'<b>优先核对的同来源图谱。</b>存在数据重用关系，不能当作独立的新发现证据；作者标签与C3的一一对应尚未建立。']], [152,215,FW-367])
add('综合解读','h3')
add('骨肉瘤CAF内部异质性及不同表达状态已有多种分类方案。F2-C3的文献定位应建立在这些既有框架之上，应以具体表达组合及其与既有亚群的关系作为主要说明内容。')
#4
newpage();add('三、其余参考文献与最终查新判断','h2')
table(['参考文献','实际研究层次','是否构成直接亚群重复证据'],[
['Huang等，2022 [7]<br/>Theranostics','CAF／Div-CAF等分层及复发相关CAF研究。','证明骨肉瘤CAF及其异质性已有研究；未提供与本C3完整定义直接对应的证据。'],
['Hu等，2026 [8]<br/>Cancer Drug Resistance','总体CAF变化；所述九基因签名来自恶性细胞。','相应比较未细分CAF。不能把恶性细胞签名用于裁定C3是否为已报道CAF亚群。'],
['Hu等，2024 [9]<br/>J Immunother Cancer','骨肉瘤PDX中的总体CAF和ECM特征。','提供共同的基质研究背景，不是与本C3一一对应的精细CAF亚群报告。'],
['Rossow等，2018 [10]<br/>Oncogene','非骨肉瘤的胶原／基质机制研究，关键过表达对象为4T1肿瘤细胞。','不属于骨肉瘤CAF亚群发现研究，不能作为F2-C3已经被报道的直接依据。']], [137,173,FW-310])
add('最终判断','h3')
box('<b>1. 是否已有相近研究？有。</b>COL5A1+间质细胞和matrix CAF是最直接先例，骨肉瘤胶原／基质相关细胞方向并非研究空白。<br/><br/><b>2. 本项目可以强调什么？</b>在这一研究方向中进一步划分并描述F2-C3的表型及胶原组成、加工相关表达组合，作为候选细分状态开展研究。<br/><br/><b>3. 目前是否适合称为首次发现的新亚群？不建议。</b>现有文献对照尚未建立F2-C3与直接相近文献群之间的明确鉴别边界。')
add('建议采用的表述','h3')
add('<b>F2-C3：本项目识别的胶原／基质相关候选细分状态。</b><br/>其主要工作内容是对已有基质相关细胞异质性的进一步描述，具体关注CD248／POSTN／VCAN表型与胶原组成及加工相关表达。')
add('说明：本次定向查新以原报告和原始文献公开的亚群定义、展示标志及作者代码为依据，未完成跨研究原始细胞的统一注释。因此，“已有相近先例”不等于已证明C3与某一群完全相同；展示组合不同也不等于已经证明是新群。未检索到同名群不能作为首次发现的依据。','small')
#5
newpage();add('参考文献','h2')
for r in refs:
 add(f'<b>[{r["number"]}] {html.escape(r["citation"])}</b><br/>{html.escape(r["title"])}<br/><link href="{r["url"]}" color="#146985">DOI: {r["doi"]}</link>','ref')
add('资料与引用说明','h3')
add('所提供清单对应[1]、[4]、[7]-[10]；补充对照为[2]、[3]、[5]、[6]、[11]。[1]标志核至补充图S2I，[2]核至正文与S3F，[3]核至CAF亚群定义，[4]保留原文三群命名，[5]与[11]核至CAF分群正文，[6]依据v1预印本及<link href="https://github.com/kidcancerlab/Osteo_atlas" color="#146985">作者公开代码</link>。本项目特征依据2026年10月6日胶原与基质分析报告及同日整合报告。','small')

def footer(c,doc):
 c.saveState();c.setFont('YH',8);c.setFillColor(C('#647b89'));c.drawString(44,H-29,'骨肉瘤F2-C3 | 文献交叉对比');c.drawRightString(W-44,H-29,'2026-10-10');c.setStrokeColor(C('#cbdde5'));c.line(44,H-40,W-44,H-40);c.line(44,40,W-44,40);c.drawString(44,26,'既往亚群 · 标志组合 · 共同点与差异 · 查新判断');c.drawRightString(W-44,26,f'{doc.page} / 5');c.restoreState()
pdf=O/'骨肉瘤F2-C3亚群_文献交叉对比报告_客户版_20261010.pdf';doc=SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=44,rightMargin=44,topMargin=58,bottomMargin=53,title='骨肉瘤F2-C3亚群文献交叉对比报告',author='CAF项目');doc.build(story,onFirstPage=footer,onLaterPages=footer)
d=fitz.open(pdf);tx=[p.get_text() for p in d];assert len(d)==5,('pages',len(d));assert all(len(t)>450 for t in tx);joined='\n'.join(tx);assert not any(s in joined for s in ['ΔR','留出','0.611','Phase79','阶段79','FDR','匹配后','可直接回复客户']);assert not '\ufffd' in joined;assert sum(len(p.get_links()) for p in d)==12;(B/'tmp/pdfs/extracted.txt').write_text(joined,encoding='utf-8');(B/'references.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2),encoding='utf-8');qa={'status':'text_passed_pending_visual','pages':5,'references':11,'customer_supplied_references':6,'additional_references':5,'links':12,'new_analysis_results_in_customer_body':False,'page_text_lengths':[len(t) for t in tx],'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()};(B/'PDF_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print('PDF',pdf);print(json.dumps(qa,ensure_ascii=False))
