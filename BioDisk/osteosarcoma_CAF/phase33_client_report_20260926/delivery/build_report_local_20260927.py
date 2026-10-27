from pathlib import Path
import json,math,sys,datetime
import pandas as pd
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT,TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent;IN=ROOT/'inputs';OUT=ROOT/'output';OUT.mkdir(exist_ok=True)
D=json.loads((IN/'report_data.json').read_text(encoding='utf-8'))
pdfmetrics.registerFont(TTFont('Body','C:/Windows/Fonts/simsun.ttc',subfontIndex=0))
pdfmetrics.registerFont(TTFont('Head','C:/Windows/Fonts/simhei.ttf'))
pdfmetrics.registerFontFamily('Body',normal='Body',bold='Head',italic='Body',boldItalic='Head')
W,H=A4;CW=W-88
ST={
 'title':ParagraphStyle('title',fontName='Head',fontSize=23,leading=32,textColor=colors.black,spaceAfter=13),
 'h1':ParagraphStyle('h1',fontName='Head',fontSize=17,leading=24,textColor=colors.black,spaceAfter=12),
 'h2':ParagraphStyle('h2',fontName='Head',fontSize=12.3,leading=18,textColor=colors.black,spaceBefore=8,spaceAfter=6),
 'body':ParagraphStyle('body',fontName='Body',fontSize=10.5,leading=17,spaceAfter=8,wordWrap='CJK',textColor=colors.black),
 'small':ParagraphStyle('small',fontName='Body',fontSize=8.4,leading=12,spaceAfter=6,wordWrap='CJK',textColor=colors.HexColor('#42484D')),
 'cell':ParagraphStyle('cell',fontName='Body',fontSize=9,leading=13,wordWrap='CJK'),
 'cellhead':ParagraphStyle('cellhead',fontName='Head',fontSize=9,leading=13,textColor=colors.white,wordWrap='CJK'),
 'num':ParagraphStyle('num',fontName='Body',fontSize=9,leading=13,alignment=TA_CENTER),
}
story=[];md=[];headings=[]
def text(s,style='body'):
 story.append(Paragraph(escape(str(s)).replace('\n','<br/>'),ST[style]));md.append(str(s)+'\n')
def h(s,level=1):
 story.append(Paragraph(escape(s),ST['h1' if level==1 else 'h2']));md.append(('#'*(level+1))+' '+s+'\n');headings.append(s)
def page(s):
 if story:story.append(PageBreak())
 h(s)
def table(headers,rows,widths=None,numeric=()):
 widths=widths or [CW/len(headers)]*len(headers)
 data=[[Paragraph(escape(str(x)),ST['cellhead']) for x in headers]]
 for row in rows:
  data.append([Paragraph(escape(str(x)).replace('\n','<br/>'),ST['num' if j in numeric else 'cell']) for j,x in enumerate(row)])
 t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#284D63')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#D9D9D9')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F2F6F8')])]))
 story.extend([t,Spacer(1,9)])
 md.extend(['| '+' | '.join(map(str,headers))+' |','|'+'---|'*len(headers)]+['| '+' | '.join(str(x).replace('\n',' ') for x in r)+' |' for r in rows]+[''])
def fig(name,maxh=255,caption=''):
 f=IN/name
 if not f.exists():raise FileNotFoundError(f)
 im=Image(str(f));scale=min(CW/im.imageWidth,maxh/im.imageHeight);im.drawWidth=im.imageWidth*scale;im.drawHeight=im.imageHeight*scale
 story.append(im);story.append(Spacer(1,5))
 if caption:text(caption,'small')
 md.append('图件 '+name+'\n')
def fmt(x,d=3):
 if x is None or (isinstance(x,float) and math.isnan(x)):return 'NA'
 return f'{x:.{d}f}'
def pv(x):
 if x is None or (isinstance(x,float) and math.isnan(x)):return 'NA'
 return f'{x:.5g}'
def row(key,**kw):
 vals=[r for r in D[key] if all(r.get(k)==v for k,v in kw.items())]
 if len(vals)!=1:raise ValueError((key,kw,len(vals)))
 return vals[0]
def statrow(layer,cls,feature):
 return row('main_program_statistics',layer=layer,class_id=cls,feature=feature)
def surfrow(layer,cls,gene):
 return row('main_surface_statistics',layer=layer,class_id=cls,symbol=gene)
def selection(layer,target,symbol,contrast):
 vals=[r for r in D['selected_comparisons'] if r['layer']==layer and r['target']==target and r['symbol']==symbol and r['contrast']==contrast]
 if len(vals)==1:return vals[0]
 vals=[r for r in vals if not r['paired']]
 if len(vals)!=1:raise ValueError((layer,target,symbol,contrast,len(vals)))
 return vals[0]
def foot(c,doc):
 c.saveState();c.setFont('Body',8);c.setFillColor(colors.HexColor('#5B6268'))
 c.drawString(44,27,'骨肉瘤 CAF 分类比较与分析型流式报告');c.drawRightString(W-44,27,str(doc.page));c.restoreState()
# 1
text('骨肉瘤 CAF 分类比较与\n分析型流式报告','title')
text('冻结分析范围：2026 年 9 月 26 日  |  复核更新：9 月 27 日','h2')
text('本研究形成了 CAF 候选表达分类。PDGFRB、MCAM 支持收缩候选表型的分析；在明确肺取材的患者比较中，GPC3 在整体 CAF 候选层升高并达到名义显著。这两类证据构成分析型流式的候选依据，分别用于表达分类和部位关联观察。')
text('核心交付是可追溯的分类、患者级比较、身份与技术证据核查，以及分析型流式抗体和必要对照。细胞集合保持冻结，未新增 UMI、基因数、线粒体比例、签名覆盖率或参考细胞数准入门槛。')
table(['结果','证据','研究解释'],[
 ['候选大类','支持层 1,224；扩展层 8,174；F1-F5 分类','表达类型与 CNV 身份等级分开呈现'],
 ['肺部位关联标志','扩展层整体 GPC3：骨 12、肺 3；P=0.008791，FDR=1','明确肺取材比较中的名义显著线索'],
 ['候选表面标志','MCAM、PDGFRB：内部比较均 P=0.03125，FDR=1；6/6 患者同向','保留为探索性组合，不因 FDR 未通过而剔除'],
 ['质量与身份核查','107 份去污染；512 基因应激；S/G2M；20 组敏感性；肺图谱映射','分群与表达证据可复核，关键身份重叠保留'],
 ],[91,211,CW-302])
h('主要结论',2)
text('扩展层大类的分群对应关系较稳定，myCAF 内部计算群的边界则明显受处理方式影响。肺参考显示平滑肌与成纤维表型重叠；完整 CNV 重跑仍存在候选与正常参考的非中性区段信号。当前数据不支持将 E-F2-C1 命名为已经确认的独立 CAF 亚型。')
text('分析型流式比较 PDGFRB 与 MCAM 联合表达，同时记录 GPC3 与 THY1，并配置免疫、内皮、上皮排除及壁细胞相关对照。E-F2-C1 的肺侧比例方向保留，但排除肺/胸壁混合取材后 P=0.07273，当前数据不支持其达到名义显著。')
text('统计口径：具备有效检验时，以原始 P<0.05 筛选名义显著候选；完整报告效应量、独立患者数及同类假设家族 FDR。描述性结果不补造 P 值。','small')
# 2
page('研究纳入九个来源的完整捕获记录')
table(['来源','样本记录','患者记录','完整捕获细胞'],[[r['dataset'],r['sample_records'],r['donor_records'],f"{r['full_called_cells']:,}"] for r in D['datasets']],[157,85,85,CW-327],numeric=(1,2,3))
text('合计 107 份来源记录、876,648 个完整 called cells，其中包括 67 份 ScPCA 记录及 40 份 GEO 记录。跨研究的患者记录数不直接相加作为独立统计样本量；推断性比较在来源内部进行。','small')
h('临床分组按原始证据保留',2)
cg={r['clinical_group']:r['sample_records'] for r in D['clinical_groups']}
table(['临床标签','记录数'],[
 ['B2 肺转移灶',cg.get('B2',0)],
 ['B1 有肺转移病史支持但时点未明',cg.get('B1_history_supported_timing_unresolved',0)],
 ['A 报告无转移',cg.get('A_reported_nonmetastasis',0)],
 ['A 候选 报告无转移但时点未明',cg.get('A_candidate_reported_no_metastasis_timing_unknown',0)],
 ['未知分组',cg.get('Unknown',0)],
 ['其他转移或复发标签',107-sum(cg.get(k,0) for k in ['B2','B1_history_supported_timing_unresolved','A_reported_nonmetastasis','A_candidate_reported_no_metastasis_timing_unknown','Unknown'])]
 ],[CW-90,90],numeric=(1,))
text('取材部位、原发或转移、疾病时序、治疗记录分别保存。治疗字段的 Yes/No 按来源原义报告，不自动改写为“化疗/未化疗”。淋巴结不归入 B2。完整逐样本台账见附件 S01。')
# 3
page('候选范围分层后形成五类表达类型')
text('最终候选母集合为 11,379 个细胞；此前 11,351 是来源尚未全部处理时的进度数。明确竞争谱系背景 3,205 个细胞单列，余下 8,174 个组成扩展探索层；其中 1,224 个组成支持候选层，6,950 个身份未决候选保留研究。')
table(['类型','表达解释','支持层','扩展层'],[
 ['F1','炎症及外膜样',478,695],['F2','收缩及活化基质',90,225],['F3','基质重塑',343,251],['F4','成骨与基质重叠',0,5304],['F5','混合及表达未定',313,1699]
 ],[45,245,95,CW-385],numeric=(2,3))
text('两层分别拟合，因此同名类别成员不要求嵌套；扩展层包含支持层细胞，不构成独立验证。支持层无独立 F4 类，不表示其中完全没有成骨基因表达。')
h('表达类型与身份等级分别记录',2)
text('冻结时，支持层中 114 个细胞具有双方法无异常支持，1,110 个为单方法支持。114 个只作证据对照，不单独强行分群。新重跑对其中 89 个单方法支持细胞所在亚群检出非中性区段；此项冲突单列，冻结层不因单次新调用而自动删减，也不称为确证 CAF。')
text('myCAF、iCAF、apCAF、基质程序和客户指定成骨、矿化程序均已评估。F1-F5 是当前候选表达分类，不要求人为制造独立 apCAF 类。SPP1、TAGLN、RUNX2 等不单独承担 CAF 身份判定。')
h('内部状态使用独立编号',2)
table(['范围','细胞数','初始计算群','本报告标记'],[
 ['支持层 F2',90,5,'S-F2-C0 至 S-F2-C4'],
 ['扩展层 F2',225,6,'E-F2-C0 至 E-F2-C5']
 ],[130,65,85,CW-280],numeric=(1,2))
text('S 与 E 编号不互相等同。后文 E-F2-C1 指扩展层内部计算群 1，共 38 个细胞、15 份样本记录、13 个来源内患者记录，按收缩及壁细胞重叠候选状态命名。')
# 4
page('F2 收缩程序及表面候选具有患者级差异')
text('在主要研究 SCPCP000017 内，F2 相对其余候选细胞呈收缩程序增强。扩展层的基质程序也保持患者间一致方向。比较在同患者、同原始部位内进行，再在患者内合并多个匹配部位。')
rows=[]
for layer,cls,f,label in [('supported','F2','program_myCAF_customer','支持 F2 收缩'),('expanded','F2','program_myCAF_customer','扩展 F2 收缩'),('expanded','F2','program_mCAF_matrix','扩展 F2 基质')]:
 r=statrow(layer,cls,f);rows.append([label,r['n_independent_donors'],fmt(r['median_delta']),f"{r['positive_donors']}/{r['n_independent_donors']}",pv(r['P']),pv(r['FDR'])])
table(['比较','患者 N','中位差值','正方向','P','FDR'],rows,[128,48,79,67,91,CW-413],numeric=(1,2,3,4,5))
text('程序效应为平均 log1p(CP10k) 评分差。名义显著与 FDR 显著分别标注；完整检验家族和全部结果保留在附件。','small')
rows=[]
for layer,g,label in [('supported','PDGFRB','支持 PDGFRB'),('expanded','PDGFRB','扩展 PDGFRB'),('supported','MCAM','支持 MCAM'),('expanded','MCAM','扩展 MCAM')]:
 r=surfrow(layer,'F2',g);rows.append([label,r['n_independent_donors'],fmt(r['median_delta']),f"{r['positive_donors']}/{r['n_independent_donors']}",pv(r['P']),pv(r['FDR'])])
table(['表面标志','患者 N','中位差值','正方向','P','FDR'],rows,[128,48,79,67,91,CW-413],numeric=(1,2,3,4,5))
text('单基因效应为 log2(患者伪批量 CP10k + 0.1) 差。MCAM 在扩展层的 FDR 未通过，但原始 P<0.05，按用户规定保留为名义显著候选。')
h('去污染后主要方向保留',2)
text('扩展层 F2 收缩程序的中位差由 +0.625 变为 +0.651，仍为 11/12 患者正方向；基质程序由 +1.352 变为 +1.236，仍为 12/12。MCAM 仍为 10/12 正方向；PDGFRB 从 11/12 降为 9/12。')
text('这些前后对照用于评价效应和方向，不把原始计数检验的 P/FDR 移植成去污染后的显著性。','small')
# 5
page('内部群的部位关联受混合取材记录影响')
text('原宽口径包含一名“Right lung/chest wall”患者，临床分组为 Unknown。该口径下 E-F2-C1 比例为骨 17.3%、肺 63.0%（9/3 名患者），P=0.036364、FDR=0.218182。明确肺取材口径为骨 17.3%、肺 61.1%（9/2 名患者），差 +43.8 个百分点，P=0.072727、FDR=0.436364，未达到名义显著。')
fig('E_F2_C1_bone_lung_proportion.png',230,'图 1 每点为一名患者，三角形为混合取材记录。分母为 F2 候选细胞；零分母不记为零比例。两种部位口径完整保留，混合取材不冒充明确 B2。')
h('该状态的标志证据来自自身内部比较',2)
rows=[]
for g in ['MCAM','PDGFRB','JAG1','THY1','FAP']:
 r=selection('expanded','subgroup_1',g,'subgroup_vs_other_myCAF')
 rows.append([g,r['n_first'],fmt(r['median_effect']),f"{int(r['positive_donors'])}/{r['n_first']}",pv(r['P']),pv(r['FDR'])])
table(['基因','患者 N','中位差值','正方向','P','FDR'],rows,[95,48,95,65,95,CW-398],numeric=(1,2,3,4,5))
text('MCAM、PDGFRB、JAG1 均在 6/6 名可配对患者中高于其余 F2，P=0.03125。THY1 在该内部比较中 P=0.0625，未达到名义筛选标准；其作为辅助通道的依据来自父类及蛋白检测用途，不能借用主类 P 值。')
# 6
page('明确肺取材的整体候选层 GPC3 名义升高')
text('剔除混合取材并重新核对患者和完整检验家族后，整体 CAF 候选层的 GPC3 在肺侧仍达到名义显著。主要研究 SCPCP000017 的独立组结果如下；效应为患者表达中位数的肺减骨。')
rows=[]
for layer,target,label in [('supported','all_CAF_candidates','支持层整体'),('expanded','all_CAF_candidates','扩展层整体'),('supported','F3','支持层 F3')]:
 r=row('strict_site_all_surface',layer=layer,target=target,symbol='GPC3',paired=False)
 rows.append([label,f"{r['n_bone']}/{r['n_lung']}",fmt(r['median_effect']),pv(r['P']),pv(r['FDR'])])
table(['GPC3 比较','骨/肺患者 N','中位差值','P','FDR'],rows,[128,100,93,93,CW-414],numeric=(1,2,3,4))
text('去污染后方向保留：支持层整体中位差由 +3.934 变为 +3.890，扩展层由 +3.431 变为 +3.437；这些为前后描述性复核，不借用原始 P 值宣称校正后表达显著。')
text('在 F2 内部，明确肺取材的 GPC3 比较 P=0.054545、FDR=1，未达到名义显著。整体候选层结果不能改写为 GPC3 标记某个肺相关 myCAF 亚型；其与 THY1 的联合蛋白表达列为待测组合。')
h('收缩候选的差异与身份重叠同时保留',2)
text('E-F2-C1 的 MCAM、PDGFRB 去污染后仍在 6/6 名患者中保持正方向，中位差为 +6.917、+3.325；壁细胞程序同样为 6/6 正方向。38 个细胞含 4 个双方法、34 个单方法无异常支持，另有 1 个新增双细胞预测。')
text('HLCA 标签为平滑肌 18、支气管周围成纤维 14、肌成纤维 2、周细胞 2、应激平滑肌 1、肺泡成纤维 1。宽部位口径下，E-F2-C1 内 MCAM、JAG1 的肺侧表达反而较低（骨 4/肺 3；均 P=0.028571、FDR=1）；这与状态比例属于不同假设，不能互相替代。')
# 7
page('成骨矿化与临床比较形成可复核的方向证据')
rows=[]
for layer,cls,feature,name in [('expanded','F4','program_Osteogenic_customer','F4 成骨'),('expanded','F4','program_Mineral_customer','F4 矿化'),('expanded','F2','program_Osteogenic_customer','F2 成骨'),('expanded','F2','program_Mineral_customer','F2 矿化'),('supported','F3','program_mCAF_matrix','支持 F3 基质')]:
 r=statrow(layer,cls,feature);rows.append([name,r['n_independent_donors'],fmt(r['median_delta']),f"{r['positive_donors']}/{r['n_independent_donors']}",pv(r['P']),pv(r['FDR'])])
table(['主要研究比较','患者 N','中位差值','正方向','P','FDR'],rows,[128,48,79,67,91,CW-413],numeric=(1,2,3,4,5))
text('扩展层 F4 成骨程序为名义显著；F2 成骨和矿化程序当前没有达到 P<0.05。不能将 F2 的收缩差异替代成骨或矿化结论。完整 16 个成骨基因、15 个矿化调节基因及内部群、临床比较结果均保留。')
h('治疗记录对应的捕获比例差异',2)
rows=[]
for r in D['clinical_proportion_nominal']:
 if r['study']=='SCPCP000023' and r['layer']=='expanded':
  rows.append([r['target'],f"{r['n_first']}/{r['n_second']}",fmt(r['mean_first']*100,3),fmt(r['mean_second']*100,3),pv(r['P']),pv(r['FDR'])])
table(['扩展层指标','No/Yes N','No 均值 %','Yes 均值 %','P','FDR'],rows,[108,63,79,79,87,CW-416],numeric=(1,2,3,4,5))
text('SCPCP000023 的治疗字段 Yes 组捕获到的扩展候选比例较低。这里的分母是完整捕获细胞，结果描述治疗标签关联，不将 Yes 改写为具体化疗，也不解释为治疗导致 CAF 减少。')
text('所有严格临床标签、部位和治疗对比均建立可评估清单。无数据或无有效重复的条目记为不可检验，不将 Unknown 填入 A，也不补造 P 值。','small')
# 8
page('肺参考显示驻留样程序及收缩表型的重叠')
text('HLCA 正式模型映射全部 8,174 个候选细胞，参考核心细胞为 584,944 个。F1 以外膜和肺泡成纤维参考表型为主；F2 更接近支气管周围成纤维与平滑肌相关表型。')
fig('HLCA_main_cohort_patient_weights.png',278,'图 2 SCPCP000017 内先对患者取类别均值，再对患者取平均。数值为参考近邻标签权重，非确证细胞起源或蛋白阳性比例。')
table(['配对患者','骨侧 AdvF+AlvF 权重','肺侧权重','肺减骨'],[
 ['CCG1106_006','0.253','0.443','+0.189'],['CCG1106_009','0.501','0.712','+0.211']
 ],[123,165,106,CW-394],numeric=(1,2,3))
text('两对患者的合计参考权重在肺侧均较高；AlvF 单独方向一对升高、一对下降。配对可能并非同期取材，结果按 N=2 的探索性方向呈现。F2 本身没有骨肺两端都有细胞的有效表达配对，不能借用整个 CAF 候选层的配对数。')
text('Gong 2022 人源 PTGS2 共表达面板为 PTGS2、IL6、CCL2、CXCL1、PTX3、TNFAIP6。已计算原始和去污染后连续评分；未把六基因面板称为完整鼠源差异签名，也未新增阳性基因数门槛。')
# 9
page('技术复核保留大类信号并识别内部群敏感性')
text('107/107 份 DecontX 与 11/11 份补充 scDblFinder 任务完成；全范围校正前后对照覆盖 8,174 个候选细胞。完整 512 基因应激评分及 S43/G2M54 评分已完成。')
fig('sensitivity_ARI_overview.png',240,'图 3 固定原始细胞、原始高变基因集合和分群参数，比较五种技术处理。ARI 仅描述一致性，不设置通过阈值。')
text('扩展层整体 ARI 为 0.816-0.848；支持层为 0.565-0.835。F2 内部计算群明显更敏感：支持层 ARI 0.121-0.317，扩展层 0.308-0.503。当前数据不支持把初始内部群边界作为稳定亚型边界。')
table(['核查','结果'],[
 ['双细胞补算','223 个候选获得新调用，其中 25 个预测为 doublet；含支持层 14 个'],
 ['污染复核','COL1A1/COL1A2、SPP1、IBSP、客户程序及关键表面候选均有前后对照'],
 ['apCAF 相关复核','CD74/HLA-DRA/HLA-DRB1 与 CIITA/CD80/CD86 纳入污染对照；未人为建立独立 apCAF 类'],
 ['批次与样本','Harmony、逐样本着色图、原始表达患者级差异及完整图例均保存'],
 ],[114,CW-114])
text('新增疑似双细胞和参考身份不确定性均单列，不自动缩减已经冻结的支持层或扩展层。')
# 10
page('完整重跑仍未解决关键候选的 CNV 身份冲突')
text('骨样本 SCPCS000614 与肺样本 SCPCS000610 的四项 inferCNV 已从完整原始计数重新计算，全部正常退出并通过独立验收。运行期间脚本保持冻结，未复用旧对象或检查点。原来源参考混有壁细胞，本次按原注释改用免疫加内皮参考，并与淋巴系参考比较；壁细胞仍保留在输入中。')
table(['样本','完整输入细胞','免疫加内皮参考','淋巴系参考'],[
 ['SCPCS000614 骨',3843,1122,26],['SCPCS000610 肺',5335,1321,153]
 ],[151,115,115,CW-381],numeric=(1,2,3))
text('均使用原生 i6 HMM、cutoff=0.1、denoise 和原随机种子。未添加参考数量、区段长度或振幅门槛。以下“有区段”指所在 HMM 计算亚群具有软件保留的非中性区段，不等同逐细胞恶性判定。')
rf=IN/'rerun_candidate_evidence_summary.tsv'
ref=pd.read_csv(rf,sep='\t');rr=ref[ref.layer.eq('expanded')&ref.target.eq('F2')]
table(['样本','扩展 F2 数','免疫加内皮：有区段','淋巴系：有区段'],[
 [sm,int(rr[rr['sample'].eq(sm)].iloc[0]['cells']),int(rr[rr['sample'].eq(sm)&rr.variant.eq('immune_endothelial_reference')].iloc[0]['subcluster_non_neutral_region']),int(rr[rr['sample'].eq(sm)&rr.variant.eq('lymphoid_reference')].iloc[0]['subcluster_non_neutral_region'])] for sm in ['SCPCS000614','SCPCS000610']
 ],[130,87,153,CW-370],numeric=(1,2,3))
h('正常参考的表现限制了二分类解释',2)
text('免疫加内皮方案中，骨参考 1,122/1,122、肺参考 1,266/1,321 个细胞所在亚群也有非中性区段。淋巴系方案中，骨参考为 0/26，肺参考为 68/153。参考行为和候选冲突均未通过改阈值消除。')
text('这次重跑确认计算可复现，但不能证明候选为正常 CAF，也不能把区段检出直接解释为恶性。骨、肺扩展候选分别 80、61 个，均有此类信号；其中支持层为 55、34 个。历史 CopyKAT 支持与新 inferCNV 区段证据分别保存。')
text('原骨样本末尾脚本报错及旧输出保留；本页使用四项正常退出的新结果。核查范围为两个关键样本，不冒称全部 107 份均换参考重跑。SCPCS000525 独立恢复任务尚未退出，未替换冻结分类。','small')
# 11
page('分析型流式以明确的双阳性候选比较为主')
text('首轮面板比较 PDGFRB 与 MCAM 联合表达，并记录 GPC3 与 THY1。前者依据收缩表达分类，GPC3 依据整体候选层的明确肺取材比较；不将不同统计比较拼成已经验证的新亚型。')
table(['用途','靶标与克隆','产品示例','安排'],[
 ['核心候选','PDGFRB/CD140b 28D4','BD 558821 PE [R3]','主面板'],
 ['核心候选','MCAM/CD146 P1H12','BD 566226 BV421 [R4]','主面板'],
 ['辅助表型','THY1/CD90 5E10','BioLegend 328139/328140 BV711 [R5]','记录表达分布'],
 ['部位关联候选','GPC3 307801','R&D FAB2119A APC [R11]','主面板，记录与 THY1 联合分布'],
 ['阴性排除','CD45 HI30','BD 555482 FITC [R7]','免疫排除'],
 ['阴性排除','CD31 WM59','BioLegend 303104 FITC [R8]','内皮排除'],
 ['阴性排除','EpCAM/CD326 9C4','BioLegend 324204 FITC [R9]','上皮排除'],
 ['重叠对照','CSPG4/NG2 LHM 2','Novus NB100-2688JF646 [R10]','另管替换 GPC3 通道观察'],
 ],[65,141,188,CW-394])
h('次级候选保留但不增加首轮必选项',2)
text('FAP 427819、FAB3715A APC [R6] 作为备用观察标志，不设为全部 F2 必须阳性。JAG1 有内部名义差异，但所核对 188331 的 MAB1277/BAM1277 产品页面标记停产 [R12]，不列入首轮现货组合。')
text('GPRC5A 在支持层整体的明确肺取材比较为 P=0.045455、FDR=1。厂家 481906/MAB5239 有转染细胞膜相关染色和胞内染色资料 [R13]，尚不能直接等同本项目未通透活细胞表面检测，列为备用。','small')
text('ACTA2、TAGLN、MYH11、RUNX2 等功能或胞内标志不作为活细胞表面通道。CSPG4 作为观察性身份对照，不新增“阳性一律排除”的规则。')
# 12
page('必要对照服务于分析读数而不扩展实验目标')
h('分析顺序与组合',2)
text('采集稳定区段 → 单细胞 → 活细胞 → CD45/CD31/EpCAM 阴性区 → PDGFRB×MCAM 四象限；在同一排除后母群观察 GPC3×THY1，并交叉记录各象限。保留单阳性与双阴性参照，不只保留双阳性事件。')
table(['组合','分析目的'],[
 ['PDGFRB+MCAM','检验收缩候选状态的核心双阳性表型'],
 ['GPC3+THY1','检验部位关联候选与辅助表型的联合表达'],
 ['THY1+MCAM 或 THY1+PDGFRB','比较辅助表型和收缩候选覆盖'],
 ['CSPG4 独立对照管','观察壁细胞相关重叠，不作为未经验证的硬排除门']
 ],[183,CW-183])
table(['必要对照','解决的问题'],[
 ['未染与活死染对照','自发荧光、背景及死细胞影响'],
 ['各荧光单染','补偿或光谱解混'],
 ['关键候选门 FMO','依据实测背景设置阴阳性边界'],
 ['已知阳性材料','确认具体克隆和批号能产生可解释信号'],
 ['CD45/CD31/EpCAM 单项检查','合并 dump 前确认每个排除标志有效'],
 ['实际消化条件下表位检查','评价表面信号保留，不自拟统一损失百分比']
 ],[171,CW-171])
text('配色草案以仪器具备对应通道为前提：PDGFRB-PE、MCAM-BV421、THY1-BV711、GPC3-APC，另配 FITC dump 和兼容的活死染；对照管用 CSPG4-JF646 替换 GPC3。APC 与 JF646 不在同管作为两个需分辨的信号。实际仪器配置尚未给出，最终通道匹配据实确认。')
text('所有抗体用量依据具体产品与实际样本滴定；不设客户未规定的固定 MFI、阳性比例或最低细胞数标准。方案围绕分析读数及其必要对照组织。')
# 13
page('关键解释边界与客户主要求核对')
text('名义显著支持候选提出，不能独立解决细胞身份。E-F2-C1 的壁细胞程序和正常肺参考相似性，与其收缩、表面标志差异同时存在。当前数据不支持将它解释为已确证的促转移 CAF 谱系。')
text('全背景 RNA 模拟门中，PDGFRB+MCAM 门内 F2 占比中位数约 0.31%，CNV 异常或冲突比例约 88.4%-94.0%。RNA 未检出并不等于蛋白阴性，这些数值不是实际流式纯度预测，但提示需保留身份对照，不能将候选标志当作 CAF 特异证据。')
text('严格 A/B1/B2 受到转移时点与治疗记录缺失的限制；骨肺部位差异不等同转移风险或因果效应。两层相互重叠，多研究患者去重未完全闭合，来源内患者数用于统计。无效重复和全零不可评估结果不补造 P 值。')
text('明确肺取材复核未支持 E-F2-C1 比例达到名义显著，但保留整体候选层 GPC3 的名义差异。整体层差异可包含组成变化；分群和标志筛选来自同一数据，按探索性证据呈现。')
table(['客户要求','本轮交付状态'],[
 ['临床分组和治疗标注','逐样本台账完成；未知及时间未明保留'],
 ['CAF 分类与程序','双层 F1-F5 分类及完整客户程序完成'],
 ['myCAF 内部分群和差异','5/6 初始群、患者级差异及处理敏感性完成'],
 ['临床组间和部位比较','可估计比较完成，名义线索与不可估计条目均交付'],
 ['配对患者','两对扩展层骨肺描述完成；F2 无有效表达配对如实记录'],
 ['成骨与矿化面板','16+15 基因完整保留并纳入比较'],
 ['肺驻留参考','当前冻结范围 HLCA 映射及 Gong 人源面板完成'],
 ['数据来源覆盖','现有 107 记录完成；Budhathoki 人源子集/去重及 GSE237070 矩阵缺口单列'],
 ['分析型流式方案','候选、克隆、组合和必要对照完成；实际湿实验尚未开展']
 ],[152,CW-152])
# 14
page('方法补充保留完整基因集与统计家族')
h('数据处理和分层',2)
text('以完整来源 called-cell 原始计数为输入，保留来源注释、CNV 调用、双细胞证据和临床元数据。初始候选复核后分别拟合支持层与扩展层。原始流程使用 logCP10k、2000 高变基因、PCA、按 study_technology 的 Harmony、邻居图和 Leiden；主分辨率为 1，边界对照为 0.5，随机种子固定。')
h('患者级比较',2)
text('基因分析采用患者和部位内原始计数聚合，归一化分母为对应细胞的原始完整文库总计数，转换为 log2(CP10k+0.1)。类间比较在同患者同部位匹配后对部位平均；临床比较将同时出现在两条件的患者单列配对，排除出独立组比较。')
text('配对采用双侧均值符号翻转：不超过 15 对时精确枚举；非配对在分配数不超过 65,536 时精确置换，其余使用固定种子的 65,536 次 Monte Carlo。细胞比例分别以完整捕获细胞、候选层或 F2 为分母；零分母记录不可评估，零分子保留。')
text('基因、程序及比例分别在预先定义的同类假设家族内进行 BH 校正。原始 P<0.05 用于名义候选筛选，FDR 完整报告。FDR<0.05 才称为校正后显著，其余 P<0.05 结果称为名义显著或探索性信号。')
h('客户指定完整程序',2)
for key,name in [('myCAF_customer','myCAF'),('iCAF_customer','iCAF'),('apCAF_customer','apCAF'),('Osteogenic_customer','成骨 16 基因'),('Mineral_customer','矿化及钙调节 15 基因')]:
 text(name+'：'+'、'.join(D['panels'][key])+'。','small')
h('质量与参考方法',2)
text('DecontX 逐样本原生自动聚类；scDblFinder 补齐缺少直接调用的来源。应激使用原始 512 基因集，细胞周期使用正式 S43/G2M54。敏感性固定原始特征和细胞集合。HLCA 使用原模型转换与 scANVI/scArches 投射、50 近邻距离加权转移，保存连续不确定度。')
# 15
page('防错指南逐项核查和证据索引')
checks=[
 ('G01-G02','恶性身份及双 CNV','双方法及 marker 证据分层；关键骨肺样本同方法换参考核查，原始失败和不可评估项保留'),
 ('G03-G04','SPP1 及破骨背景','SPP1 不单独定身份；破骨与髓系竞争谱系单列'),
 ('G05','壁细胞重叠','RGS5/CSPG4/NOTCH3/DES 程序及 HLCA 参考纳入；F2 不等同纯 myCAF'),
 ('G06','apCAF 假象','污染、双细胞和 CIITA/CD80/CD86 证据可查；不强造独立 apCAF 群'),
 ('G07-G08','肺来源和队列构成','HLCA/Gong 与来源内临床比较完成；部位与转移概念分开'),
 ('G09-G10','整合与周期','Harmony、样本着色、原始表达复核及正式周期敏感性完成'),
 ('G11-G12','推断单位和 bulk','患者作为统计单位；bulk 非本轮任务'),
 ('G13-G16','污染 双细胞 应激 QC','107 去污染、11 补充双细胞、512 应激及来源 QC；未新增入口门槛'),
 ('G17-G19','表面标志 组合 抗体','胞内标志区分、至少双阳性候选、阴性排除和厂家克隆核查完成'),
 ('G20-G21','消化与状态命名','消化表位检查写入必要对照；表达状态命名保留，实际实验不冒充已执行')
]
table(['编号','项目','核查结果'],checks,[61,127,CW-188])
h('附件统计与数据文件',2)
table(['索引','内容'],[
 ['S01','107 来源临床与取材台账'],
 ['S02','F1-F5 主类程序、表面候选和全基因统计'],
 ['S03','myCAF 内部及临床表达比较 全部原始 P/FDR'],
 ['S04','临床比例、明确肺取材全家族复核及患者清单'],
 ['S05','肺参考逐细胞标签、不确定度、患者权重和骨肺配对'],
 ['S06','去污染前后、双细胞及应激周期证据'],
 ['S07','20 组分群敏感性、原始特征与细胞对齐验收'],
 ['S08','针对性换参考全细胞与候选对照'],
 ['S09','E-F2-C1、GPC3 证据、完整核对表及报告原始数据']
 ],[62,CW-62])
# 16
page('文献与试剂来源')
references=[
 ('R1','HLCA 官方模型映射方法','https://docs.scvi-tools.org/en/stable/tutorials/notebooks/hub/query_hlca_knn.html'),
 ('R2','Gong 等 2022 人源 PTGS2 共表达依据','https://pmc.ncbi.nlm.nih.gov/articles/PMC9830653/'),
 ('R3','BD CD140b 28D4 PE 558821','https://www.bdbiosciences.com/en-nl/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/pe-mouse-anti-human-cd140b.558821'),
 ('R4','BD CD146 P1H12 BV421 566226','https://www.bdbiosciences.com/en-us/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/bv421-mouse-anti-human-cd146.566226'),
 ('R5','BioLegend CD90 5E10 BV711','https://www.biolegend.com/en-gb/products/brilliant-violet-711-anti-human-cd90-thy1-antibody-15886?GroupID=BLG5826'),
 ('R6','R&D FAP 427819 APC FAB3715A','https://www.rndsystems.com/products/human-fibroblast-activation-protein-alpha-fap-apc-conjugated-antibody-427819_fab3715a'),
 ('R7','BD CD45 HI30 FITC 555482','https://www.bdbiosciences.com/en-eu/products/reagents/flow-cytometry-reagents/research-reagents/single-color-antibodies-ruo/fitc-mouse-anti-human-cd45.555482'),
 ('R8','BioLegend CD31 WM59 克隆及流式用途','https://www.biolegend.com/fr-fr/products/brilliant-violet-711-anti-human-cd31-antibody-14948'),
 ('R9','BioLegend EpCAM 9C4 FITC 324204','https://sandbox-punchout.biolegend.com/Default.aspx?Clone=9C4&ID=6680'),
 ('R10','Novus CSPG4 LHM 2 JF646','https://resources.rndsystems.com/products/documents/datasheets/NB100-2688JF646.pdf'),
 ('R11','R&D GPC3 307801 APC FAB2119A','https://www.rndsystems.com/products/human-glypican-3-apc-conjugated-antibody-307801_fab2119a'),
 ('R12','R&D JAG1 188331 MAB1277 停产状态','https://www.rndsystems.com/products/human-mouse-rat-jagged-1-antibody-188331_mab1277'),
 ('R13','R&D GPRC5A 481906 染色资料','https://resources.rndsystems.com/pdfs/datasheets/mab5239.pdf')
]
for key,title,url in references:
 story.append(Paragraph(f'{key}　<link href="{escape(url)}" color="#245A7A">{escape(title)}</link>',ST['body']))
 md.append(f'{key} [{title}]({url})')
text('客户依据：骨肉瘤课题相关需求；骨肉瘤转移灶CAF生信分析防错指南。本文采用用户 2026 年 9 月 26 日最终明确的统计和分析型流式范围。','small')
text('来源核查日期：2026 年 9 月 26 日。产品用途依据厂家资料，不表示本项目已完成实际组织染色。原始计算、完整 P/FDR、失败记录及数据完整性证据随统计附件保留。','small')
f=OUT/'骨肉瘤CAF分类比较与分析型流式报告_20260926.pdf'
doc=SimpleDocTemplate(str(f),pagesize=A4,rightMargin=44,leftMargin=44,topMargin=40,bottomMargin=43,title='骨肉瘤 CAF 分类比较与分析型流式报告',author='CAF 项目分析',pageCompression=1)
doc.build(story,onFirstPage=foot,onLaterPages=foot)
(OUT/'骨肉瘤CAF分类比较与分析型流式报告_20260926.md').write_text('\n'.join(md),encoding='utf-8')
from pypdf import PdfReader
reader=PdfReader(str(f));texts=[p.extract_text() or '' for p in reader.pages]
(ROOT/'pdf_text_check.json').write_text(json.dumps({'pages':len(texts),'characters_by_page':[len(t) for t in texts],'page_starts':[t[:160] for t in texts],'reference_results_included':rf.exists()},ensure_ascii=False,indent=2),encoding='utf-8')
print('PDF_WRITTEN',str(f),'PAGES',len(texts))


