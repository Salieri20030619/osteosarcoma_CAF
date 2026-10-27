from pathlib import Path
import pandas as pd,json,html,hashlib,shutil
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
import pymupdf as fitz
B=Path('C:/Users/31937/CAF_C3_subdivision_value_20261010');R=B/'results';D=B/'delivery';T=R/'tables';(D/'output/pdf').mkdir(parents=True,exist_ok=True);(B/'tmp/pdfs').mkdir(parents=True,exist_ok=True)
for name in ['tables','resources','logs']:shutil.copytree(R/name,D/name,dirs_exist_ok=True)
(D/'code').mkdir(exist_ok=True)
for name in ['analyze_subdivision.py','analyze_subdivision_initial.py','spatial_assay_limited.py','verify_subdivision.py']:shutil.copy2(B/name,D/'code'/name)
shutil.copy2(B/'literature_update.json',D/'resources/literature_update.json');shutil.copy2(R/'COMPLETE.json',D/'COMPLETE.json')
s=pd.read_csv(T/'composition_statistics.tsv',sep='\t');m=pd.read_csv(T/'matching_statistics.tsv',sep='\t');a=pd.read_csv(T/'assay_limited_spatial_association_statistics.tsv',sep='\t');v=pd.read_csv(T/'assay_limited_spatial_prediction_summary.tsv',sep='\t')
pdfmetrics.registerFont(TTFont('YH','C:/Windows/Fonts/msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('YHB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0));pdfmetrics.registerFontFamily('YH',normal='YH',bold='YHB',italic='YH',boldItalic='YHB');C=colors.HexColor;W,H=A4;FW=W-92
styles={k:ParagraphStyle(k,fontName='YHB' if k in ['title','h2','h3'] else 'YH',fontSize=sz,leading=lead,spaceAfter=sp,wordWrap='CJK',textColor=C('#175d77') if k in ['title','h2','h3'] else C('#2a3c4a')) for k,sz,lead,sp in [('title',20,28,14),('h2',14.5,22,11),('h3',11.3,17,8),('body',10.2,17,9),('small',8.6,13,7),('cell',9.1,14.2,0)]};story=[]
def P(t,k='body'):return Paragraph(t,styles[k])
def add(t,k='body'):story.append(P(t,k))
def tab(headers,rows,widths):
 t=Table([[P('<b>'+str(x)+'</b>','cell') for x in headers]]+[[P(str(x),'cell') for x in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),C('#eaf3f6')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),('LINEBELOW',(0,0),(-1,-1),.4,C('#d5e3e9'))]));story.extend([t,Spacer(1,11)])
def box(t):
 z=Table([[P(t)]],colWidths=[FW]);z.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),C('#eef5f7')),('BOX',(0,0),(-1,-1),.5,C('#c5dae3')),('LEFTPADDING',(0,0),(-1,-1),13),('RIGHTPADDING',(0,0),(-1,-1),13),('TOPPADDING',(0,0),(-1,-1),12),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([z,Spacer(1,13)])
def page():story.append(PageBreak())
# 1
add('F2-C3细分价值检验与客户答复','title');add('2026-10-10  |  回答10月6日F2-C3是否已有骨肉瘤研究先例','small')
box('<b>核心判断：</b>骨肉瘤中相近的胶原／基质相关群体此前已有报道。我们的F2-C3可强调的是<b>进一步细分得到的表达差异线索</b>：胶原组成有不同侧重，在一部分可严格匹配的细胞中仍保留加工偏向。目前尚不足以宣称首次发现一个独立新亚群，也未证实更强的组织结构定位能力。')
add('三个方向的实际检验结果','h2')
tab(['方向','本轮结果','当前可强调的范围'],[
['胶原组成','C3对C1：III／V型平均差+0.493，3/4患者正向；IV型平均差-0.540，4/4较低。','有不同表达侧重的探索性线索，最一致的是C1的IV型评分较高。'],
['加工偏向','严格匹配COL5程序及文库量后，4名患者均保留PCOLCE／LOXL2较高表达，平均差+0.611。','支持部分C3细胞在相当胶原背景下仍有加工增量；尚不是普遍属性。'],
['组织信息增量','19名患者的空间RNA中存在小幅条件关联；胶原组成的留出预测平均ΔR²仅+0.0022，原发9例为-0.0198。','存在关联，未建立稳定的组成定位优势；不能称为已验证结构分工。']], [78,237,FW-315])
add('可直接回复客户','h2')
add('查新后，骨肉瘤中类似的胶原／基质相关亚群已有研究，最接近的是COL5A1+间质细胞和matrix CAF。因此，胶原型CAF这个方向本身不是空白。我们10月6日的F2-C3目前可以强调进一步细分的意义：它与C1的胶原表达组成存在不同侧重，且在胶原程序相当的一部分细胞中仍保留加工相关表达增量。现阶段适合将其表述为有数据线索的候选细分状态；还不能说它就是某篇文章的同一个群，也不宜直接宣称首次发现独立新亚群。')
add('本报告给出细分的实际证据及限制，作为此前亚群查新报告的增量核查。全部比较采用既有标签，不重新定义C3以追求更好结果。','small')
#2
page();add('胶原组成与加工偏向：有线索，证据仍属探索性','h2')
add('同患者直接比较C3与C1','h3')
q=s[(s.matrix=='corrected')&(s.minimum==1)&(s.reference=='F2_C1')];names={'III_V':'III／V型胶原','IV':'IV型胶原','balance':'III／V减IV的偏向'}
tab(['指标','平均差','正向/N','P','FDR'],[[names[r.panel],f'{r.mean_effect:+.3f}',f'{int(r.positive)}/{int(r.N)}',f'{r.P:.3f}',f'{r.FDR:.3f}'] for r in q.itertuples()],[175,84,74,82,FW-415])
add('差值均为C3减C1。III／V程序为COL3A1、COL5A1、COL5A2、COL5A3；IV程序为COL4A1-A6。评分为基因平均ln(1+CP10k)，两程序之差表达转录侧重，不代表胶原质量或沉积比例。')
add('<b>关键敏感性：</b>当每组至少5个细胞时，仅剩1名患者，其III／V减IV差为-0.062，未保留正向，N=1不检验。因此，目前不能把“C3偏III／V、C1偏IV”写成普遍、已确立的生物学分工。')
add('控制共同胶原背景后的加工比较','h3')
tab(['比较方案','可匹配C3细胞／34','患者N及方向','加工差','P / FDR'],[
['严格匹配：两项变量均在0.5个标准差范围内','12 / 34','4例，4/4正向','+0.611','0.125 / 0.125'],
['无距离限制的最近邻敏感性','34 / 34','7例，7/7正向','+0.940','0.0156 / 0.0313']], [179,90,94,62,FW-425])
add('严格匹配在同患者、同部位的其他Fibro中，为每个C3细胞选择3个相近对照，允许复用对照。匹配变量为已发表COL5四基因程序和文库量，该四基因程序与PCOLCE／LOXL2没有重叠。匹配后两变量的总体标准化均值差约-0.003、-0.056；原始矩阵结果同向，加工差+0.606。')
add('无距离限制方案仍存在COL5程序及文库量不平衡（标准化均值差约+0.416、+0.394），其FDR显著不能代替严格匹配结果。严格匹配只覆盖12/34个C3细胞，结论应限于可比较的部分。')
add('4名患者的双侧精确检验最小可达P为0.125。不能因未达到0.05而抹去观察到的方向，但当前结果不属于名义显著。C3仍是既有筛选得到的候选，本轮未新增显著性筛选。','small')
#3
page();add('空间检验：小幅关联未转化为稳定组成定位优势','h2')
add('数据为GSE293065的26张空间切片、19名患者；原发9名另作分析。以CD248／POSTN／VCAN平均表达代表C3相关RNA特征，在控制总体matrix CAF、其他细胞背景、文库量及空间坐标后，检验其与胶原程序的关联。')
add('全部19名患者的条件相关','h3')
q=a[(a.scope=='all_patients')&(a.control=='broad_CAF_available')]
tab(['胶原读出','患者平均ρ','同向/N','P','FDR'],[[names[r.endpoint],f'{r.mean_effect:.3f}',f'{int(r.positive)}/{int(r.N)}',f'{r.P:.3g}',f'{r.FDR:.3g}'] for r in q.itertuples()],[153,95,75,90,FW-413])
add('关联方向可重复且通过本轮FDR校正，但效应较小，C3特征同时关联III／V和IV。不能把小幅相关直接解释成特定胶原结构的准确定位。原发9名患者的组成偏向平均ρ=0.033，8/9正向，P=0.0234、FDR=0.0234，仍按小样本探索性结果呈现。')
add('整名患者留出后的预测增量','h3')
q=v[v.control=='broad_CAF_available'];rows=[]
for scope in ['all_patients','primary_only']:
 for ep in ['III_V','IV','balance']:
  r=q[(q.scope==scope)&(q.endpoint==ep)].iloc[0];rows.append(['全部19例' if scope=='all_patients' else '原发9例',names[ep],f'{r.mean_delta_R2:+.4f}',f'{int(r.positive)}/{int(r.N)}'])
tab(['评估范围','读出','平均ΔR²','改善患者数'],rows,[100,192,109,FW-401])
add('ΔR²比较“总体CAF及背景模型”与“加入C3相关分数后的模型”，按患者等权汇总。组成读出的增量在总体接近零、在原发病例平均为负。再加入可测收缩背景后，总体组成ΔR²约+0.00006，原发为-0.0125，结论不变。留出训练折相互重叠，这些预测增量只作描述，不附加P值。')
add('<b>判断：</b>目前没有证据支持“细分C3后就能更准确定位特定胶原组成”的稳定优势。这里评估的是RNA代理分数，尚不是对真实C3细胞或胶原蛋白结构的直接定位。')
#4
page();add('此前是否有人做过：已有直接相近报道','h2')
tab(['骨肉瘤先例','与本C3的关系'],[
['Zheng等，2024，COL5A1+ MSC [1]','COL5A1、COL16A1、EMILIN1、ADAMTS2是最直接的胶原型参照。我们的C3具有这些共同特征，本轮匹配考察的是它们相当时是否仍有加工差异。'],
['Zheng等，2024，matrix CAF [2]','已报道基质相关CAF框架，POSTN、FN1及多种胶原与C3相交。共有功能背景不能作为首次发现的依据。'],
['Zeng等，2026卷期，FNDC1+ CAF [3]','既有ECM合成相关群，代表组合与本C3不完全吻合；作为相近对照保留，未证明与C3完全同一。'],
['Liu等，2021，C1/C2/C3 CAF [4]','该文C3_CAFs以TOP2A／MKI67增殖特征定义，与本项目C3定义不同。聚类编号相同不能证明重复。']], [171,FW-171])
add('可采用的定位','h3')
box('<b>在既有胶原／基质相关细胞框架中识别的候选细分状态，具有胶原组成差异及局部加工偏向的表达线索。</b><br/>这一定位有具体数据依据。当前不能进一步写成首次发现独立新亚群、已经证实功能分工，或已经形成明确的结构定位优势。')
add('主要证据边界','h3')
add('本轮组成和匹配结果来自发现队列内的小样本比较；此前GSE237070两例原发病例的加工模块一升一降，不能被本轮同队列匹配结果覆盖。完整空间参照缺少未测得的RUNX2、MYH11，已保留不可评估结果；本页前的空间结果来自明确标记的可测基因敏感性，成骨背景用既有五基因，收缩敏感性用ACTA2／TAGLN。缺失基因未按零处理。')
add('参考文献与可复核记录','h3')
for i,citation,doi in [(1,'Zheng et al. J Cancer Res Clin Oncol. 2024. COL5A1+间质细胞研究。','10.1007/s00432-024-05787-2'),(2,'Zheng et al. J Hematol Oncol. 2024. 骨肉瘤单细胞与空间图谱。','10.1186/s13045-024-01598-7'),(3,'Zeng et al. Molecular Cancer. 2026;25:14（2025-12在线发表）。','10.1186/s12943-025-02516-2'),(4,'Liu et al. Front Oncol. 2021. 未治疗骨肉瘤单细胞图谱。','10.3389/fonc.2021.709210')]:add(f'[{i}] {citation} <link href="https://doi.org/{doi}" color="#176c86">{doi}</link>','small')
add('空间数据：<link href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293065" color="#176c86">GSE293065</link>。完整文献交叉表见前一版三项核查报告。新增52条统计、132条组成患者效应、92条匹配细胞对比、156条空间相关和168条留出预测记录均已复核。原假设、缺测修订、原始P／FDR及未改善结果全部保留于证据包。','small')

def footer(c,doc):
 c.saveState();c.setFont('YH',8);c.setFillColor(C('#607784'));c.drawString(46,H-29,'CAF研究 | C3细分价值检验');c.drawRightString(W-46,H-29,'2026-10-10');c.setStrokeColor(C('#ccdce3'));c.line(46,H-40,W-46,H-40);c.line(46,40,W-46,40);c.drawString(46,26,'胶原组成 · 加工偏向 · 空间信息 · 文献先例');c.drawRightString(W-46,26,f'{doc.page} / 4');c.restoreState()
pdf=D/'output/pdf/F2-C3细分价值检验与客户答复_20261010.pdf';doc=SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=46,rightMargin=46,topMargin=58,bottomMargin=54,title='F2-C3细分价值检验与客户答复',author='CAF项目分析');doc.build(story,onFirstPage=footer,onLaterPages=footer)
doc=fitz.open(pdf);texts=[p.get_text() for p in doc];assert len(doc)==4,('Unexpected pages',len(doc));assert all(len(x)>500 for x in texts);assert not any('\ufffd' in x for x in texts);(B/'tmp/pdfs/extracted.txt').write_text('\n\n'.join(texts),encoding='utf-8');qa={'status':'text_checks_passed_pending_visual_review','pages':len(doc),'page_text_lengths':[len(x) for x in texts],'links':sum(len(p.get_links()) for p in doc),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()};(D/'PDF_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print('PDF',pdf);print(json.dumps(qa))
