from pathlib import Path
import re,json,hashlib,shutil,html
import pandas as pd,numpy as np,pymupdf
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'unified_full_report';W.mkdir(exist_ok=True);D=B/'review_bundle_final';F=W/'figures';F.mkdir(exist_ok=True)
target=Path('C:/Users/31937/Downloads/骨肺CAF亚群_图文分析与Lung04推荐报告.pdf');backup=W/'previous_17page_report.pdf'
if not backup.exists():shutil.copyfile(target,backup)
original=pymupdf.open(B/'full_report_integration/original_customer_report.pdf');extracted=[]
for i in range(1,7):
 for info in sorted(original[i].get_image_info(xrefs=True),key=lambda x:x['bbox'][1]):
  x=info['xref'];assert x>0;pix=pymupdf.Pixmap(original,x);path=F/f'figure_{len(extracted)+1:02d}.png';pix.save(path);extracted.append(dict(original_page=i+1,xref=x,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
assert len(extracted)==11,len(extracted)
notes=json.loads((B/'REPORT_INTERPRETATION.json').read_text(encoding='utf-8'));stats=pd.read_csv('C:/Users/31937/caf121_report_work/resources/focus_functional_audit.tsv',sep='\t');dec=pd.read_csv('C:/Users/31937/caf121_report_work/resources/node_decisions_ZH.tsv',sep='\t').set_index('父群或子群')
for name,file in [('YH','msyh.ttc'),('YHB','msyhbd.ttc')]:pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
pdfmetrics.registerFontFamily('YH',normal='YH',bold='YHB',italic='YH',boldItalic='YHB')
ink=colors.HexColor('#243b49');teal=colors.HexColor('#176c7a')
styles={k:ParagraphStyle(k,fontName='YHB' if k in ['title','h','sub'] else 'YH',fontSize=size,leading=lead,spaceAfter=gap,wordWrap='CJK',textColor=ink,keepWithNext=k in ['h','sub']) for k,size,lead,gap in [('title',20,29,14),('h',14,21,10),('sub',11,17,7),('body',10,16,8),('small',8.4,12.5,6),('ref',8.2,12.2,7)]}
story=[]
def p(t,k='body'):story.append(Paragraph(t,styles[k]))
def h(t):p(t,'h')
def page(t):
 if story:story.append(PageBreak())
 h(t)
def table(rows,widths,fs=8.5):
 st=ParagraphStyle('cell',parent=styles['small'],fontSize=fs,leading=fs*1.45,spaceAfter=0);t=Table([[Paragraph(html.escape(str(v)),st) for v in r] for r in rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5f0f2')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.65,teal),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#d6e2e6')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.append(t);story.append(Spacer(1,8))
def figure(path,title,caption,maxheight=245):
 from PIL import Image as PI
 iw,ih=PI.open(path).size;scale=min(511/iw,maxheight/ih);story.append(KeepTogether([Paragraph(title,styles['sub']),Image(str(path),width=iw*scale,height=ih*scale),Spacer(1,5),Paragraph(caption,styles['small']),Spacer(1,7)]))
captions={1:'骨病灶2,788细胞、15父群；肺病灶3,635细胞、14父群。各病灶在PCA/Harmony表示上独立计算UMAP，随机种子20261026，n_neighbors=15，min_dist=0.3。坐标不跨图比较，分离不等同于功能特异或患者重复。',2:'左列为来源，右列为互斥类别：浅灰为其他父群（骨2,657、肺3,402），青色为目标父群其余细胞，红/黄为指定子群。Bone_07=66+35(s1)+30(s2)=131；Lung_04=197+36(s2)=233。',3:'展示27个预设功能及谱系标志物。点大小为RNA检出比例，颜色为患者等权平均log1p(CP10k)的逐基因Z-score，范围[-2,2]。先在患者内平均，再对来源内患者记录等权；图3–5分别标准化，颜色不跨图比较，也不表示检验显著性。',4:'肺病灶14个父群，计算规则同图3。Lung_04呈胶原及基质表达，并伴收缩和成骨标记；个别标志物的患者配对结果见补充表S1。',5:'全部13个保留信号群的标志物表达，计算规则同图3。Lung_05_s2、Lung_05_s4等同时检出PTPRC、PECAM1等竞争谱系标记。',6:'行是15个骨病灶父群，列是16条来源内患者记录（SCPCP000017：12；SCPCP000023：4）。数字为细胞数，颜色为行内患者占比。单患者主导例：Bone_00为517/521，Bone_01为492/492，Bone_03为239/241，Bone_04为209/210。',7:'14个肺病灶父群、23条来源内患者记录：NCI 17条，另外两个来源各3条。Lung_00为1,027/1,055、Lung_01为372/372、Lung_02为361/362；Lung_07、Lung_11、Lung_12主要由单患者贡献。Lung_04的233细胞覆盖NCI 10条和SCPCP000017 1条记录。',8:'灰色为目标和对照各≥1细胞，青蓝为各≥5细胞，橙色为严格身份层各≥5细胞。均是选定来源内有效配对患者数；红线为主分析建议下限6位。Lung_04的主分析和身份层分别为7、5位。',9:'采用min5配对。点为患者效应，菱形为等权均值（Lung_04为红色）；Lung_11无有效配对，无点。父群对同病灶其余候选，子群对同父群其余细胞；效应为患者内目标减对照的log1p(CP10k)程序均分差，不是倍数变化，不跨层级直接排名。',10:'骨肉瘤文献四类面板[1]的min5患者平均效应。Lung_04的FNDC1为+0.32，7/7同向，P=0.015625、FDR=0.06875；MYH11(-0.03)、CD36(-0.04)、C7(-0.01)未见增强。Bone_07偏MYH11方向（+0.94、4/4同向、P=0.125）。NA为无有效min5配对。',11:'Lung_04严格身份层有效配对由7位降至5位。患者6546和7978的目标细胞仍充足（34→24、14→13），但对照由16→4和5→3，低于5细胞条件；剩余5位的基质方向一致。'}
titles={1:'骨、肺病灶父群UMAP',2:'来源构成与重点父子候选映射',3:'骨病灶15个父群标志物气泡图',4:'肺病灶14个父群标志物气泡图',5:'全部13个保留信号群标志物气泡图',6:'骨病灶逐患者细胞贡献',7:'肺病灶逐患者细胞贡献',8:'不同门槛下的有效配对患者支持',9:'重点信号群逐患者功能效应',10:'文献功能面板交叉比较',11:'Lung_04身份敏感性配对减少归因'}
def f(n,maxheight=245):figure(F/f'figure_{n:02d}.png',f'图{n}｜'+titles[n],captions[n],maxheight)
p('骨肉瘤骨肺病灶CAF亚群<br/>基质相关信号与Lung_04推荐报告','title');p('固定队列分析 · 患者级效应 · HLCA外部参照','sub');p('更新日期：2026-10-27　｜　完整队列6,633细胞，主分析6,423细胞','small')
h('一、执行摘要与核心推荐')
p('优先推荐肺病灶父群Lung_04作为后续验证的ECM/胶原生成相关CAF候选，并分层保留其余12个信号群。Lung_04在NCI内有7位有效配对患者，基质、胶原交联、纤维型胶原及FNDC1面板均7/7同向增强；严格身份层仍有5位有效配对患者，基质与FNDC1方向一致。')
p('HLCA补充支持其以成纤维细胞样身份为主，且基质程序高于正常肺参考背景。候选推荐仍以原骨肺队列的患者支持和功能一致性为核心，外部图谱用于解释身份与组织背景。')
p('表1｜重点候选及主分析定位','sub');table([['候选','层级','细胞','有效配对N','定位'],['Lung_04','父群',233,7,'优先：基质/胶原程序一致增强，伴成骨相关转录特征'],['Bone_07','父群',131,4,'保留：骨侧基质与收缩相关线索'],['Bone_07_s1','子群',35,3,'保留：同父群内胶原相关线索'],['Lung_04_s2','子群',36,1,'保留：同父群内基质偏向，配对患者少']],[86,40,40,70,275])
p('阅读顺序：研究队列与准则 → UMAP及标志物 → 患者支持与功能效应 → HLCA身份与正常肺参照 → 方法与综合推荐 → 参考文献和补充材料。','small')
page('二、研究队列与患者级评估准则')
p('本研究在骨、肺病灶内分别沿用15个父群（Bone_00–Bone_14）和14个父群（Lung_00–Lung_13），并评估其子群。NCI中210个双细胞阳性标记细胞仅保留供溯源，不进入主分析UMAP和标志物图。本固定队列不含NCC及正常组织参考细胞；HLCA作为独立外部参照使用。')
p('表2｜队列构成与来源内患者记录','sub');table([['病灶/来源','完整细胞','主分析细胞','患者记录'],['Bone / SCPCP000017','2,477','2,477',12],['Bone / SCPCP000023',311,311,4],['Lung / SCPCP000017',792,792,3],['Lung / SCPCP000023','1,523','1,523',3],['Lung / NCI','1,530','1,320',17],['合计','6,633','6,423','Bone 16；Lung 23']],[230,85,90,106])
p('不同来源独立编码，以上是来源内患者记录，不能直接当作全球去重人数；样本和建库重复按原donor关系合并。')
p('推荐规则优先考虑同来源至少6位有效配对患者，目标与对照每位患者各≥5细胞（min5），至少75%患者同向，目标功能程序原始双侧P&lt;0.05。达到10位及以上可增加支持，仍结合效应、来源及模型适用性综合判断。每侧≥10细胞（min10）、严格身份层（identity_min5）及留一患者分析用于敏感性评估。')
p('6位是本项目操作性下限，不是普适文献标准或充分功效保证。双侧精确符号翻转检验中，N=5时最小P为0.0625，N=6时为0.03125；增加单患者细胞数不能替代独立患者重复[8,9]。严格身份层的患者数与效应作为稳健性信息保留。')
page('三、病灶内分群结构与来源构成');f(1);f(2)
page('三、标志物图谱：骨、肺父群');f(3);f(4)
page('三、全部信号群与分层保留');f(5,220)
sig=['Bone_00','Bone_07','Bone_07_s1','Bone_07_s2','Bone_13','Lung_04','Lung_04_s2','Lung_05_s2','Lung_05_s4','Lung_06_s2','Lung_07','Lung_11','Lung_12']
roles=['单患者主导的基质线索','骨侧基质与收缩','胶原线索；加工不一致','来源间方向不一致','胶原/成骨；单患者','优先推荐','父群内基质偏向','免疫/内皮标记需复核','免疫/内皮标记需复核','少患者；稳定性有限','收缩/血管壁线索','同患者对照不足','胶原加工；单患者']
rows=[['亚群','细胞/覆盖','评估来源','min1','min5','min10','身份N','保留定位']]
for node,role in zip(sig,roles):
 a=dec.loc[node];rows.append([node,f"{int(a['细胞'])}/{int(a['覆盖记录'])}",str(a['评估来源']).replace('SCPCP000017','017').replace('SCPCP000023','023'),int(a['最少1细胞有效患者']),int(a['最少5细胞有效患者']),int(a['最少10细胞有效患者']),int(a['身份敏感性有效患者']),role])
p('表3｜全部13个信号群（017、023分别为对应SCPCP来源）','small');table(rows,[82,48,55,34,34,38,36,184],7.4)
p('覆盖指跨来源至少1个目标细胞的记录；min1/min5/min10及身份N均为选定来源的有效配对患者数。带_s的编号为父群内子群，其余为父群；父子群有包含关系，不可相加。','small')
page('四、患者贡献构成');f(6);f(7)
page('四、患者支持与核心功能方向');f(8);f(9)
page('五、文献交叉与身份敏感性');f(10);f(11);p('Lung_04呈FNDC1/ECM相关增强，伴成骨程序重叠；Bone_07更偏MYH11收缩/血管相关参考方向。文献面板提供功能参照，不等同于完整复现原文亚群。','small')
page('五、Lung_04核心效应与敏感性统计')
p('表4｜NCI内患者配对结果；效应为目标减对照程序均分差','sub');rows=[['分析层/程序','N','同向','平均效应','原始P','BH-FDR']];auditrows=[]
for mode,title in [('main_min5','主分析'),('identity_min5','身份层')]:
 for ep,label in [('program__mCAF_matrix','基质'),('program__collagen_crosslinking','交联'),('program__fibrillar_I_III_V','纤维型胶原'),('audit__FNDC1_paper','FNDC1'),('audit__collagen_processing_nonCOL','非COL加工')]:
  a=stats[(stats.node_id=='Lung_04')&(stats.source=='NCI')&(stats['mode']==mode)&(stats.endpoint==ep)].iloc[0];rows.append([title+'/'+label,int(a.N),f'{int(a.positive_patients)}/{int(a.N)}',f'{a.mean_effect:+.3f}',f'{a.P:.5f}',f'{a.FDR:.5f}']);auditrows.append(a.to_dict())
table(rows,[170,40,55,76,85,85],9)
p('主分析为名义显著，以上BH-FDR均未达到0.05。基质、交联、纤维型胶原及FNDC1面板均7/7同向，原始P=0.015625；非COL胶原加工为6/7同向，效应+0.228，P=0.03125、FDR=0.10577。身份层中部分程序5/5升高，P=0.06250受精确检验分辨率限制。')
p('成骨程序主分析效应+0.328，7/7同向，P=0.015625；身份层+0.409，5/5同向。当前结果不支持实际矿化增强。非COL加工结果提供了结构胶原以外的补充线索，但功能终点仍需实验检验。')
page('六、HLCA参照：确认候选的主要身份')
p('HLCA（Human Lung Cell Atlas，人类肺细胞图谱）是整合多个研究的肺组织单细胞参考。完成本队列功能评估后，先利用其健康核心图谱判断候选更接近成纤维细胞、周细胞还是平滑肌细胞，再与正常肺成纤维细胞比较基质表达。外部细胞不并入原分群和UMAP。')
figure(B/'figures/identity.png','图12｜HLCA健康核心图谱身份参照','每位患者先计算候选细胞的身份比例，再对患者等权平均；候选每患者至少5细胞。颜色图例列出全部身份类别。',270)
p('Lung_04的7位可评估患者均来自NCI，平均成纤维细胞样比例94.7%；Lung_04_s2的2位为100.0%。Bone_07与Bone_07_s1各4位可评估患者均来自SCPCP000017，成纤维细胞样比例82.2%和75.6%，平滑肌样比例17.8%和24.4%。')
p('保持参考参数冻结，从200轮续训至500轮，整体身份大类一致率97.4%；Lung_04、Bone_07、Bone_07_s1、Lung_04_s2分别为96.1%、96.9%、91.4%、100.0%。骨侧候选同时保留基质生成和收缩/平滑肌样身份信息。')
page('六、正常肺参照：基质程序超过组织固有背景')
p('肺组织本身含有胶原和基质表达。因此，这一步直接衔接前述功能结果，比较候选信号是否高于正常肺泡和外膜成纤维细胞的背景。')
figure(B/'full_report_integration/healthy_background_main.png','图13｜候选与健康肺成纤维细胞的程序分数差','候选患者中位分数减去健康参考的研究等权均值，单位为平均log1p(CP10k)。红色表示候选更高，蓝色表示参考更高。候选≥5细胞，参考每供者每类≥20细胞。',240)
p('四项候选的基质、纤维型胶原和交联程序均高于两类健康参照。Lung_04相对健康肺泡参照差值为+1.704、+1.772和+0.610；非COL加工在外膜参照中的中位百分位为0.759，说明部分加工表达也属于正常外膜成纤维细胞特征。肿瘤邻近肺及10/50细胞门槛下主要方向保持较高。')
p('表5｜外部参照患者数与原主分析配对人数的区别','small');table([['候选','原min5有效配对N','外部参照：候选≥5细胞N'],['Lung_04',7,7],['Bone_07',4,4],['Bone_07_s1',3,4],['Lung_04_s2',1,2]],[150,170,191])
p('外部比较只要求候选细胞数充足，不要求同患者对照；原主分析要求两侧均满足门槛。Lung_04在外部比较≥10细胞时仍有7位患者，原主分析min10则为5位有效配对，两者不能互换。','small')
page('七、统计方法、外部参照方法与质量控制')
p('<b>统计单位：</b>同来源、同病灶、同患者的目标与对照均值差，每位有效患者等权。父群对同病灶其余候选细胞，子群对同父群其余细胞；同患者样本及建库重复按原donor关系合并。采用未整合log1p(CP10k)及预设基因集评分，FNDC1等面板按原文基因列表复核。')
p('<b>检验与多重校正：</b>N≤16采用双侧精确符号翻转检验，更大N使用65,536次蒙特卡洛模拟；N&lt;2或全零不补造P。原始P&lt;0.05用于名义线索，BH-FDR按病灶、来源、父/子层级、比较背景、分析层和终点家族六个键分别计算。program子集1,817条有效检验分属13个非空BH家族；Lung/NCI/parent/within_site/singlet_min5/program含92条检验。新增audit文献面板属独立家族，不与program、rank或不同层级混池。')
p('<b>身份敏感性：</b>沿用成纤维支持层；NCI限定作者CopyKAT标注为diploid的细胞，并对明确免疫/内皮组合共表达作敏感性排除。')
p('<b>HLCA映射：</b>使用v1.1预训练模型和健康核心584,944个细胞的潜在空间，2000基因原始计数输入，查询批次按来源设置，以50近邻加权投票。主映射200轮，保持参考参数冻结追加300轮复核；NCI匹配1,631/2,000个模型特征，缺失特征仅作模型补零。')
p('<b>正常表达参照：</b>健康肺140位供者、17项研究；肿瘤邻近肺33位供者、5项研究。每研究/供者/条件/类型固定种子最多抽100细胞，采用study::donor标识。共同唯一基因378个，各程序覆盖率≥80%；按全库原始计数归一后计算分数，先按供者汇总，再对参考研究等权。主设置候选≥5、参考每类≥20细胞，10/50为敏感性；跨数据集比较报告分数差与百分位，不给出跨数据集P值。')
p('<b>数值核验：</b>既有15,792条功能比较重算一致，UMAP含6,423细胞，13群细胞数及min1/min5/min10人数一致；1,817条program检验的FDR复算最大误差1.11×10<super>-16</super>。HLCA补充的47项核验通过，包含原始计数、患者汇总、程序分数、单基因及冻结参数检查。')
page('八、综合推荐与后续验证')
p('优先推进Lung_04作为ECM/胶原生成相关CAF候选，其余12个信号群全部保留并分层跟进。原队列的患者内配对效应与患者一致性构成推荐基础；HLCA补充主要身份和正常肺背景解释。')
p('表6｜Lung_04推荐依据','sub');table([['证据层','结果'],['本队列患者支持','NCI主分析7位有效配对，严格身份层5位；跨来源共覆盖11条记录。'],['功能一致性','基质+0.849、交联+0.306、纤维型胶原+0.904、FNDC1+0.318，均7/7同向。'],['非COL加工','主分析+0.228，6/7同向；身份层+0.307，5/5同向。'],['HLCA身份','可评估患者的成纤维细胞样比例平均94.7%；延长冻结训练后身份大类一致率96.1%。'],['正常肺背景','主要基质程序高于肺泡、外膜及肿瘤邻近肺参考；提供外部组织背景解释。']],[100,411])
p('Lung_04相对健康肺泡参照的典型myCAF程序差值仅+0.036，成骨程序+0.800、矿化+0.072；Bone_07_s1矿化差值-0.039。因此采用ECM/胶原生成相关定位，保留成骨相关转录特征，不直接命名为典型收缩型myCAF或已证实矿化的细胞群。')
p('<b>证据范围：</b>当前尚无独立外部患者队列验证候选结论；同批数据的多面板及文献相似性提供互补支持。HLCA供者不计入候选患者数，跨数据集技术差异影响绝对分数，骨侧属于跨器官参照。CopyKAT diploid、标志物或映射投票均不能单独证明纯CAF或非恶性身份。')
p('<b>验证路径：</b>补充独立患者或来源样本，在组织原位确认成纤维身份及基质蛋白表达；检验胶原分泌、纤维沉积与组织结构，可结合染色或SHG。若研究化学交联机制，再补充生化定量或力学读出。RNA结果不直接证明矿化、耐药或转移机制。')
p('探索性疾病状态参照及影响解释的敏感性结果完整保留在补充材料S2–S3。','small')
page('参考文献')
raw=original[9].get_text();parts=re.split(r'(?=\[\d+\])',raw[raw.index('[1]'):])
for part in parts:
 if not part.strip():continue
 clean=' '.join(part.split());clean=re.sub(r'(?<=[\u4e00-\u9fff]) (?=[\u4e00-\u9fff])','',clean);p(html.escape(clean),'ref')
p('[10] Sikkema L, et al. An integrated cell atlas of the lung in health and disease. Nature Medicine, 2023. DOI: 10.1038/s41591-023-02327-2。HLCA模型：Zenodo 7599104；完整表达图谱：CELLxGENE数据集9f222629-9e39-47d0-b83f-e08d610c7479。','ref')
page('补充材料S1｜Lung_04标志物配对结果')
p('沿用原报告补充表S1，评估来源NCI，采用min5患者配对口径。');table([['基因','有效N','正向患者','平均效应','原始P','FDR'],['POSTN',7,'4/7','-0.163982','0.625000','0.772059'],['ACTA2',7,'2/7','-0.023889','0.859375','0.984375'],['TAGLN',7,'2/7','+0.029591','0.921875','1.000000'],['LOXL2',7,'1/7','-0.093160','0.546875','0.710374']],[75,45,70,107,107,107])
h('S2–S3：探索性纤维化状态参照的用途')
p('IPF（Idiopathic Pulmonary Fibrosis，特发性肺纤维化）具有肺组织瘢痕化及细胞外基质沉积特征。HLCA中的IPF成纤维细胞用于探索候选是否伴随更广泛的纤维化相关程序，不作为CAF筛选或Lung_04推荐门槛，也不表示本项目患者患有IPF。')
p('HLCA论文Fig.6h列出CCL2、COL1A1、CTHRC1、MMP19、SERPINE1、HIF1A等相关基因。S2比较六基因均分，S3去除COL1A1并展示单基因证据，以检查高分是否主要由胶原推动。')
p('完整抽样参照共24,983个基质细胞，含健康、肿瘤邻近肺及IPF；IPF覆盖52位供者、4项研究。归一化及供者汇总沿用正文。健康与IPF检验限定同研究、同类型且每侧≥3供者，采用双侧Mann-Whitney检验；BH-FDR在同研究、同细胞门槛下的程序/类型对比家族内计算。候选与外部参考仅作描述比较。')
page('补充材料S2｜IPF六基因功能参照');figure(B/'figures/ipf.png','补充图S1｜参考与候选的IPF六基因分数','点为患者/供者，黑线为中位数；健康/IPF点为肺泡成纤维细胞。参考每供者每类≥20细胞，候选每患者≥5细胞，括号给出人数。单位：平均log1p(CP10k)。',280)
for text in notes['ipf']:p(text)
page('补充材料S3｜单基因及COL1A1敏感性');figure(B/'figures/markers.png','补充图S2｜基因表达与检测比例','颜色为平均log1p(CP10k)，点面积为表达比例。前6列为IPF相关基因，其余为基质和myCAF标志物；先按供者汇总，参考再按研究等权。',270)
for text in notes['gene_evidence']:p(text)
nc=pd.read_csv(D/'tables/IPF_without_COL1A1_background.tsv',sep='\t');rows=[['候选','去COL1A1后健康百分位','去COL1A1后IPF百分位']]
for node in ['Lung_04','Bone_07','Bone_07_s1','Lung_04_s2']:
 a=nc[nc.node_id==node].set_index('reference_group');rows.append([node,f"{a.loc['Healthy | Alveolar fibroblasts','median_reference_percentile']:.3f}",f"{a.loc['IPF | Alveolar fibroblasts','median_reference_percentile']:.3f}"])
table(rows,[145,183,183]);p('补充表S2｜使用肺泡成纤维细胞作参照。六基因高分主要受COL1A1贡献，当前数据不支持完整典型IPF状态的认定。','small')
def footer(c,d):
 c.setStrokeColor(teal);c.setLineWidth(.6);c.line(42,810,553,810);c.setFillColor(teal);c.setFont('YH',8);c.drawString(42,819,'骨肺CAF亚群研究｜基质信号与候选推荐');c.setFillColor(colors.HexColor('#657780'));c.drawString(42,25,'2026-10-27　|　患者级证据与外部参照');c.drawRightString(553,25,str(d.page))
staged=W/'full_report_review.pdf';SimpleDocTemplate(str(staged),pagesize=(595,842),leftMargin=42,rightMargin=42,topMargin=47,bottomMargin=43).build(story,onFirstPage=footer,onLaterPages=footer)
doc=pymupdf.open(staged);doc.set_toc([[1,'执行摘要',1],[1,'研究队列与准则',2],[1,'分群与标志物图谱',3],[1,'患者支持与功能结果',6],[1,'文献与敏感性',8],[1,'Lung_04核心统计',9],[1,'HLCA身份与正常肺参照',10],[2,'正常肺背景比较',11],[1,'方法与质量控制',12],[1,'综合推荐',13],[1,'参考文献',14],[1,'补充材料',15],[2,'S1 标志物配对结果',15],[2,'S2 IPF参照',16],[2,'S3 COL1A1敏感性',17]]);doc.set_metadata({'title':'骨肺CAF亚群_图文分析与Lung04推荐报告','author':'CAF研究分析','subject':'患者级功能效应、HLCA正常肺参照与候选推荐'});doc.saveIncr();texts=[p.get_text() for p in doc]
for i,pag in enumerate(doc):pag.get_pixmap(matrix=pymupdf.Matrix(1.1,1.1)).save(W/f'page_{i+1:02d}.png')
qa=dict(status='STAGED',pages=len(doc),source_figures=extracted,body_font='Microsoft YaHei',statistics_recalculated=False,primary_statistic_rows=auditrows,all13_retained=all(n in ''.join(texts) for n in sig),reviewer_QA_separate=True,section_order=[t[:95] for t in texts],sha256=hashlib.sha256(staged.read_bytes()).hexdigest())
(W/'UNIFIED_BUILD_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print('PAGES',len(doc));print('\n'.join(f'{i+1}: '+t[:70].replace('\n',' ') for i,t in enumerate(texts)))
