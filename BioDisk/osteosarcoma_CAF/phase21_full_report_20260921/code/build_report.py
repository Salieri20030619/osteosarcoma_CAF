from pathlib import Path
import html,re,json,hashlib
import pandas as pd
from PIL import Image as PILImage
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Spacer,PageBreak,Image,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
B=Path('C:/codex_test');O=B/'CAF_full_report_20260921';D=O/'delivery';F=O/'figures'
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('Head',r'C:\Windows\Fonts\simhei.ttf'))
st={'body':ParagraphStyle('body',fontName='CN',fontSize=10.8,leading=17,spaceAfter=8,splitLongWords=False),'title':ParagraphStyle('title',fontName='Head',fontSize=21,leading=28,spaceAfter=12),'head':ParagraphStyle('head',fontName='Head',fontSize=14,leading=20,spaceBefore=9,spaceAfter=9),'sub':ParagraphStyle('sub',fontName='Head',fontSize=11.5,leading=17,spaceBefore=8,spaceAfter=6),'cell':ParagraphStyle('cell',fontName='CN',fontSize=9.3,leading=13.5,splitLongWords=False),'note':ParagraphStyle('note',fontName='CN',fontSize=9,leading=13.5,spaceAfter=8,textColor=colors.HexColor('#46535c'),splitLongWords=False)}
for style in st.values():
 style.wordWrap='CJK'
 style.splitLongWords=True
story=[];md=[];figures=[]
def esc(t):
 t=str(t).replace('\u2212','-');t=re.sub(r'(?<=[\u4e00-\u9fff])(?=[A-Za-z0-9])|(?<=[A-Za-z0-9])(?=[\u4e00-\u9fff])',' ',t)
 return html.escape(t)
def p(t,style='body'):story.append(Paragraph(esc(t),st[style]));md.append(t+'\n')
def h(t):p(t,'head')
def sub(t):p(t,'sub')
def page():story.append(PageBreak())
def table(headers,rows,widths):
 data=[[Paragraph(esc(x),st['cell']) for x in row] for row in [headers]+rows];t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#eaf0f4')),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#58758a')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#d0d9df')),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.append(t);story.append(Spacer(1,8))
 md.append('| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)+'\n')
def fig(path,caption,height=235):
 path=Path(path);w,hh=PILImage.open(path).size;scale=min(505/w,height/hh);im=Image(str(path),width=w*scale,height=hh*scale);im.hAlign='CENTER';story.append(im);p(caption,'note');figures.append(str(path));md.append('!['+caption+']('+str(path)+')\n')
summary=pd.read_csv(O/'inputs/cohort_summary.tsv',sep='\t');qc=pd.read_csv(O/'inputs/ScPCA_QC_loss_decomposition.tsv',sep='\t')
R=pd.read_csv(B/'CAF_discovery_stats_20260921/tables/C7_C9_primary_results.tsv',sep='\t');M=pd.read_csv(B/'CAF_discovery_stats_20260921/tables/within_myCAF_program_tests.tsv',sep='\t');G=pd.read_csv(B/'CAF_discovery_stats_20260921/tables/discovery_marker_tests.tsv',sep='\t')
def n(v):return f'{int(v):,}'
# 1 Overall narrative
p('骨肉瘤CAF候选状态完整分析报告','title');p('数据选择  QC与CNV身份判定  三层探索  患者级统计  分析型流式\n2026年9月21日','note')
p('本项目从人源骨肉瘤组织单细胞数据出发，先处理细胞质量、双细胞和肿瘤/间质身份，再建立严格层、支持层及广义拓展层开展联合探索。发现队列支持C7与C9具有不同的表达程序：C7相对偏myCAF，C9相对偏iCAF。当前结果适合推进小规模分析型流式，确认候选表面标志物的表达和共表达。')
fig(F/'01_analysis_workflow.png','图1  分析主线。CNV与谱系证据先于候选分层；联合分群后按患者比较表达特征，最后转化为分析型流式的候选验证任务。',235)
table(['核心问题','本轮结果'],[['分析范围','8个来源、52条样本记录；299,483个QC后细胞'],['三层候选','严格4,099个；支持54,997个；广义拓展140,135个'],['主发现','同患者直接比较显示C7的myCAF程序更高、C9的iCAF程序更高，均有探索性统计支持'],['实验落点','先开展小规模分析型流式，评估FAP、THY1等候选标志，不将双阳性直接等同于C7/C9']],[ ninety:=90,421])
p('阅读顺序：数据与分组 → QC → CNV方法及输入可评估性 → 三层定义 → 联合分群和程序 → 患者统计 → 流式与证据解释。全部细胞数均为分析记录；推断性统计另按独立来源与患者映射限定样本量。','note')
# 2 datasets
page();h('1 数据集选择与临床分组')
p('队列选择覆盖骨病灶与肺转移病灶，纳入已取得且通过文件验收的人源非PDX表达数据。原始计数、基因标识、条码和临床来源分别留档；同患者多组织、多个建库或技术分片保留关联，不按多个患者计数。')
purpose={'GSE152048':'骨原发、复发和肺转移发现队列；作者组织矩阵','GSE162454':'6例骨原发队列；补充任何转移标签对照','GSE250015':'有无转移标签的骨原发记录；患者复用对应另列','GSE270231':'肺部组织记录；同供者多记录及未知身份分开处理','GSE299023':'补充非PDX骨肉瘤组织；扩充表达状态覆盖','GSE318840':'8例FLEX样本，覆盖初治、治疗后、复发及肺转移','SCPCP000017':'5份骨/肺标本；含同患者多个肺标本','SCPCP000023':'7份补充标本；含骨肺配对线索及其他肺标本'}
table(['来源','记录数','QC细胞','选择理由与用途'],[[r.dataset,str(r.samples),n(r.QC_cells),purpose[r.dataset]] for r in summary.itertuples()],[101,43,69,298])
sub('临床标签按取材和病史分开记录')
table(['分组','定义','当前处理'],[['A','无肺转移患者的骨病灶','需有肺转移史及随访依据；不由“Primary”自动推定'],['B1','已有肺转移患者的骨原发灶','需明确转移器官及相对取样时间；任何转移标签另命名'],['B2','肺转移灶','以来源明确的肺转移取材为依据；保留治疗信息'],['同患者骨肺','同一患者的不同部位标本','合并同患者重复标本时保留部位；不自动视为同期取样']],[58,174,279])
p('本轮候选状态统计重点是患者内群间差异。GSE162454的Yes/No保留为“任何转移标签”，ScPCA保留官方取材部位、病程与治疗原文。主统计使用上海、广西和河内三个独立招募来源，具体纳入说明见统计章节。','note')
# 3 QC
page();h('2 QC先检查输入质量 再保留可分析细胞')
p('表达输入核对维度、条码唯一性、基因映射、有限非负整数计数及文库和。早期7份样本沿用已冻结QC对象；后续17份非PDX组织、16份新增GEO和12份ScPCA按各阶段流程处理。四部分保留细胞依次为28,542、125,212、101,487和44,242，合计299,483。')
p('后续GEO及ScPCA附加QC采用样本内自适应阈值：在log10基因数/UMI上用中位数减3倍未缩放MAD设下限，并设至少200基因、300 UMI；线粒体比例上限按中位数加3MAD并限制于5%至20%。高复杂度先标记，不机械按固定比例删除。')
fig(F/'02_scpca_qc.png','图2  新增12份ScPCA样本的四步保留数量。S602等为SCPCS000602等的缩写；filtered、官方processed、官方singlet及项目附加QC是不同步骤。',205)
table(['样本','filtered','processed','singlet','最终'],[[r['sample'].replace('SCPCS000','S'),n(r['filtered']),n(r['official_processed']),n(r['official_singlets']),n(r['retained'])] for r in qc[qc['sample'].isin(['SCPCS000794','SCPCS000795','SCPCS000812','SCPCS000813','SCPCS000789'])].to_dict('records')],[ ninety:=91,105,105,105,105])
p('S794和S812在官方processed阶段已大幅减少；S795、S789则分别有432和925个官方singlet低于项目300 UMI下限。S813有4,065个官方singlet低于该下限。阈值造成的保留数量与原始组织可用细胞数分开解释，低于门槛不直接代表生物学上没有相关细胞。','note')
p('GEO以Scrublet初筛并复核异常阈值，部分新增样本追加scDblFinder复核；ScPCA要求官方processed成员且scDblFinder为singlet，再应用附加QC。后续所有候选层均保留已有双细胞排除标记。','note')
#4 CNV pipeline
page();h('3 CNV将间质表达与拷贝数证据结合')
p('先按免疫、内皮、上皮、成纤维、成骨、血管壁、破骨等面板标记谱系，再使用同样本标志清楚的免疫/内皮细胞作为正常参考。CNV输入采用原始counts与基因坐标，不采用HVG子集或Harmony校正表达。Fib/Osteo/Mural混合程序和较高双细胞评分参与参考筛选。')
table(['步骤','实际实施'],[['参考训练与校准','后续GEO/ScPCA至少30个训练参考、20个留出校准细胞；ScPCA每条常染色体至少30个符合表达条件的基因。早期阶段保留原参考与方法版本'],['inferCNV','后续实现用GRCh38/Ensembl98坐标，cutoff 0.1、窗口101、denoise开启、HMM关闭；由去噪信号计算平均绝对log2偏离和染色体偏离'],['信号分类','阈值结合留出校准细胞的99%分位数、median+5 MAD与0.05下限；高信号还需至少2条染色体绝对中位偏离≥0.1，或1条≥0.25；其余区分低信号、模糊或未评估'],['CopyKAT','后续实现用正常参考，ngene.chr=5、win.size=25、KS.cut=0.1；输出diploid、aneuploid或未评估，保留基因映射审计'],['历史方法溯源','GSE270231的GSM8337654第二方法为SCEVAN；该记录沿用其历史判定，其余既有样本第二方法为CopyKAT。未将历史方法全部重写成CopyKAT']],[108,403])
fig(F/'03_CNV_states.png','图3  全部QC细胞的两方法判定组合。low/no-high表示至少一种低且没有高；no-low-evidence包含未评估/模糊情形。图中细胞不是全部CAF，比例是样本内CNV证据构成。',220)
#5 actual evidence
page();h('4 真实CNV信号与方法分歧分别保留')
p('以GSE162454的GSM4952365为例展示实际染色体信号。该样本CNV输入6,546个细胞，包括346个训练参考与85个留出校准细胞。下图使用原计算输出，固定随机种子从训练参考、严格候选和双方法高信号组各抽取30个细胞。')
fig(F/'04_actual_CNV_example.png','图4  inferCNV去噪表达的染色体log2中位信号。左为正常训练参考，中为严格层候选，右为两方法均高的细胞；颜色为相对表达偏离，不是直接DNA拷贝数。每组分别有346、260、2,356个可抽取细胞。',290)
table(['inferCNV判定','CopyKAT aneuploid','CopyKAT diploid','CopyKAT未评估'],[['low_signal','3,052','917','29'],['elevated_signal','2,356','0','20'],['ambiguous_signal','169','0','3']],[147,121,121,122])
p('同一示例中，3,052个细胞出现inferCNV低、CopyKAT高的分歧。因此流程保留每个细胞的方法结果，而不是仅采用其中一个工具的阳性/阴性标签。双低还需经过谱系、双细胞及环境RNA等检查，才可进入严格候选集合。')
sub('环境RNA与竞争谱系')
p('既有样本按阶段补齐DecontX；ScPCA有11份具备环境RNA分析条件。严格层要求校正前后均有成纤维证据，并检查成骨、血管壁、免疫/内皮/上皮等竞争程序。支持层和广义层将部分检查改为风险标记，以探索较宽候选范围；原始counts保留用于表达比较。')
p('本示例用于展示CNV证据链和分歧处理，不以单个样本替代全队列结果；完整逐样本判定留在分析记录中。','note')
#6 eligibility
page();h('5 CNV输入不足与QC损失分开记录')
fig(F/'05_reference_calibration.png','图5  新增ScPCA中5份可评估样本的留出正常细胞低信号比例，均通过该阶段90%的参考检查。留出细胞参与信号阈值校准与质量检查，该图不是外部分类准确率。',170)
p('新增12份ScPCA中，5份完成inferCNV和CopyKAT，10项计算成功。其余7份在当时输入规则下未启动双方法CNV，具体如下。')
table(['样本','QC细胞','训练/留出','未满足条件'],[['S611','3,778','39 / 10','留出不足；参考合计49，比最低总需求少1'],['S794','29','0 / 0','参考不足；最低常染色体基因覆盖为6'],['S795','189','0 / 0','未选出满足规则的正常参考'],['S812','59','0 / 0','未选出满足规则的正常参考'],['S813','5,017','0 / 0','未选出满足规则的正常参考'],['S786','4,210','0 / 0','未选出满足规则的正常参考'],['S789','379','0 / 0','未选出满足规则的正常参考']],[61,62,83,305])
p('“未选出参考”表示细胞未通过当时正常谱系、表达复杂度、竞争程序及双细胞等组合条件，不表示样本中必然没有正常细胞。S786的4,210个官方singlet全部通过附加QC，仍因参考条件而不能运行CNV；S813也不是单纯总细胞数不足。')
p('样本记录继续保留，CNV状态记为未评估。记录保留不等于必然有候选细胞：S794在广义候选筛选后仍为0，其余样本按实际分子证据纳入广义层并标记。现阶段没有将未评估当作CNV阴性，也没有为增加数量改写原严格结果。')
#7 layers
page();h('6 三个分析层级对应不同证据强度')
table(['层级','主要纳入条件','用途'],[['严格层','保留冻结FibLike范围；通常≥500基因，原始与DecontX后Fib检测均≥2；双CNV低且参考通过，严格竞争谱系/CTA/簇风险与双细胞排除','作为较保守身份基线'],['支持层','广义候选中至少一种CNV低、另一种无高信号，且既有参考检查通过；不再强制原始与校正后共同通过，也不以原簇风险一票否决','作为本轮候选发现与重点统计层'],['广义拓展层','≥300基因；原始Fib≥2，或Fib≥1且ECM≥3；排除强竞争谱系及既有双细胞。保留低/高冲突与未定细胞；双高或只有高且没有低支持者不纳入','观察边界状态和敏感性，逐细胞保留风险标记']],[66,325,120])
p('广义层的强竞争筛选采用相关标志≥3个且UMI比例>1%（破骨特异标志≥2个）；成骨与血管壁不在该步统一排除。严格集合以既有结果纳入，最终关系为严格层包含于支持层、支持层包含于广义层。')
table(['层级','候选细胞','有候选的样本记录','独立拟合群数','原始面板myCAF'],[['严格','4,099','17','10','2,867'],['支持','54,997','38','18','19,212'],['广义拓展','140,135','51','24','60,180']],[70,87,141,87,126])
fig(F/'06_layers_sources_and_risks.png','图6  左为各来源三层数量（对数型刻度）；右为各层风险标记比例，标记可重叠。广义层有72,151个低/高冲突、8,398个无低信号证据细胞。支持层中5,345个带血管壁程序，3,093个缺少环境RNA检查。',205)
p('表中myCAF统一按原始计数检测规则计算，用于层间比较；不能与旧版要求原始和校正后同时阳性的1,926个候选直接比较增量。','note')
#8 integration
page();h('7 联合分析得到18个支持层表达群')
p('候选筛选后，从完整来源矩阵取15,391个共同实测Ensembl基因，按完整文库总量计算logCP10k。建模固定HVG 2,000、PCA 30、Harmony研究×技术校正（theta=2）、15近邻和Leiden分辨率0.5；严格层、支持层、广义层分别重拟合。客户程序测量使用完整来源基因空间，不把共同交集外的基因当作零表达。')
fig(F/'07_supported_umap.png','图7  支持层联合UMAP。左突出C7和C9，右标记来源；背景包含其他16群。UMAP用于表达结构展示，不将距离或视觉分离作为显著性检验。',225)
fig(F/'08_program_heatmap.png','图8  支持层18群的客户程序概览。先按群内已知供者求均值，再在程序列内标准化用于显示；该色阶不用于计算P值。',285)
p('三个层级的群编号独立产生，不能把其他层的C7或C9直接当作支持层同名群。共享细胞上，支持层与广义层分群ARI为0.739，严格层与广义层为0.153；因此结合对应关系和患者内程序方向解释层间一致性。','note')
#9 candidate programs genes
page();h('8 C7与C9形成可优先复核的候选状态')
table(['特征','支持层C7','支持层C9'],[['细胞与面板','2,509个；2,240个原始myCAF面板阳性','3,754个；2,686个原始myCAF面板阳性'],['严格层共同细胞','374个','2,391个'],['血管壁程序标记','7.3%','0.7%'],['主要表达线索','FN1、COL6A3、TIMP1、CALD1、POSTN、TPM2等','DCN、SFRP2、LUM、COL3A1、SERPINF1、MGP等'],['本轮定位','偏myCAF表达的候选状态','相对更偏iCAF表达，并有矿化相关程序的候选状态']],[92,209,210])
fig(B/'CAF_discovery_stats_20260921/figures/patient_program_effects.png','图9  相对同患者其余候选细胞的程序差异。每点为一位患者，颜色表示三个主统计来源；黑线为中位差。',220)
p('程序使用客户原始面板：myCAF为ACTA2、TAGLN、MYH11、MYL9、TPM2、CNN1；iCAF为IL6、CXCL12、PDGFRA、DPT、HAS1；apCAF为CD74、HLA-DRA、HLA-DRB1。原始面板阳性要求至少80%基因实测、至少2个基因计数≥1，标签允许重叠。')
p('成骨采用客户16基因，包含RUNX2、SP7、ALPL、IBSP、SPP1、BGLAP等；矿化/钙调节采用客户15基因，包含MGP、ENPP1、STC1、ANXA2/5、CASR、TRPM7和OPG/RANKL相关基因。原始客户名称保留；S100A8/S100A9作为髓系鉴别线索另看，不据其阳性直接认定矿化。','note')
#10 stat primary
page();h('9 患者级统计支持不同表达程序')
p('主统计限定为GSE152048、GSE162454和GSE318840：原研究分别为上海第六人民医院11例、广西医科大学6例、河内Vinmec8例的独立招募来源。每项比较使用实际可评估患者数；未知身份及跨来源对应未解决的记录不加入主合并N。')
table(['相对同患者其余候选细胞','患者数/正向','中位差','P','FDR'],[['C7 myCAF','13 / 12','+0.661','0.00122','0.00439'],['C7 成骨','13 / 12','+0.0795','0.000488','0.00220'],['C7 矿化','13 / 13','+0.117','0.000244','0.00146'],['C9 iCAF','15 / 15','+0.489','0.0000610','0.000687'],['C9 矿化','15 / 15','+0.140','0.0000610','0.000687'],['C9 成骨','15 / 11','+0.139','0.0413','0.0790']],[166,78,70,98,99])
p('上表FDR按支持层18群×5程序的90项家族校正。C9成骨为名义显著的探索性支持；C7成骨差异虽然方向一致，但中位效应较小。')
fig(F/'09_C7_C9_paired.png','图10  同一患者C7减C9的直接程序差。myCAF、iCAF及成骨/矿化各有10例满足两群每侧≥10个细胞。',205)
table(['两群直接比较','中位差 C7减C9','P','FDR'],[['myCAF','+0.337','0.00977','0.0244'],['iCAF','-0.388','0.00195','0.00977'],['成骨','-0.00076','0.922','0.922'],['矿化','-0.00622','0.846','0.922']],[162,145,102,102])
p('直接比较支持C7的myCAF更高、C9的iCAF更高；当前不支持将两群定义为成骨或矿化程度显著不同的两个亚型。','note')
#11 within and robustness
page();h('10 myCAF内部保留部分差异与稳健性')
p('仅保留原始面板myCAF阳性细胞，再比较某群与同患者其余myCAF，避免把“较多myCAF细胞”与“myCAF内部程序差异”混为同一问题。')
focus=M[(M.layer=='supported')&M.cluster.isin([7,9])&M.program.isin(['myCAF_customer','iCAF_customer','Osteogenic_customer','Mineral_customer'])].sort_values(['cluster','program'])
names={'myCAF_customer':'myCAF','iCAF_customer':'iCAF','Osteogenic_customer':'成骨','Mineral_customer':'矿化'}
table(['群与程序','患者数','中位差','P','FDR'],[[f'C{r.cluster} {names[r.program]}',str(r.n_patients),f'{r.median_delta:+.3f}',f'{r.P:.5g}',f'{r.FDR:.5g}'] for r in focus.itertuples()],[151,70,95,97,98])
p('C7矿化和C9 iCAF仍达到FDR显著；C7 myCAF增强及C9矿化为名义显著，成骨程序在此范围内均未显著。直接比较两群的myCAF阳性细胞时有8例可评估患者：7/8例C7 myCAF更高、8/8例C9 iCAF更高；该小样本比较按探索性方向呈现。')
sub('敏感性与重要不一致结果')
p('每侧要求从10个提高到20个细胞后，C7 myCAF/成骨/矿化与C9 iCAF/矿化仍达到各自FDR显著。逐一去掉一个主队列，C9 iCAF和矿化仍保持FDR显著；去掉GSE152048后，C7只剩7例，myCAF和矿化方向保留且名义显著，FDR未过。')
p('C7与C9的直接myCAF比较为8/10例同向，精确符号检验P=0.109，而符号秩检验P=0.00977，说明显著性同时利用效应大小信息。C9 iCAF相对C7为10/10例同向，两种检验均支持。')
p('全部三层52群、260项程序及776条描述性候选基因均保留检验记录。基因候选另按层构成FDR家族，不与程序混池；程序、基因和临床主张分别登记，不把某项支持转移给其他结论。')
#12 flow
page();h('11 RNA筛查转向小规模分析型流式')
fig(F/'10_surface_screen.png','图11  FAP+THY1在已知供者/部位记录中的RNA门控表现。目标为各群的原始myCAF阳性细胞；覆盖、门内已知目标比例及CNV高信号混入均为逐记录中位数，不是实测流式纯度。',195)
p('C7在16条有足量目标细胞的记录中，覆盖中位数31.2%、门内已知目标比例3.6%；C9在18条记录中分别为46.5%和11.6%。此前16个表面候选及全部136个单/双组合均已筛查；分母为51个有广义候选样本的全部QC细胞，原始目标门控计数已与来源矩阵核对。')
table(['实验模块','首轮目的与输出'],[['样本质量','记录取材部位、治疗背景、消化条件、活率与有效细胞量；逐样本保留结果'],['基础分析','活死染与单细胞识别；评估CD45/CD31/EPCAM排除标志的分布和保留数量'],['候选表面表达','先看FAP、THY1单标及共表达是否可分辨；报告阳性比例、强度分离及二维分布'],['技术对照','未染、单色及必要FMO对照；记录背景、补偿/解混及门限依据'],['身份衔接','候选阳性群结合独立转录/蛋白标志复核，评估血管壁和恶性间质细胞混入'],['后续决定','信号可分辨且具备身份支持后，再优化组合并评估分选；不清晰者先调整marker']],[92,419])
p('本轮实验定位是小规模分析型流式，而非直接高纯度分选。FAP/THY1双阳性不直接命名为C7或C9；抗体克隆、通道、用量与门限依据仪器、消化条件及现场滴定确定。')
#13 scope and methods
page();h('12 统计方法与结论适用范围')
sub('统计实施')
p('冻结群标签后，每位患者计算群内减群外的平均logCP10k程序差值；每侧至少10个细胞、面板覆盖至少80%。采用双侧精确Wilcoxon配对符号秩检验，差值保留12位小数，并列秩使用整数化秩分布计算，零差剔除。中位差区间由按研究分层的5,000次患者bootstrap获得。')
p('程序主家族为支持层90项；严格层50项、广义层120项作为分层敏感性；myCAF内部比较另按相同家族大小校正。两群直接比较按5类程序校正。不可评估项参与家族时用P=1占位，交付表显示NA。P<0.05但FDR未显著记为名义显著；小于10例的比较保留探索性定位。')
sub('本轮主结论')
p('发现队列识别出C7和C9两个具有不同程序特征的候选状态，并获得患者级探索性统计支持。证据支持C7相对偏myCAF、C9相对偏iCAF；两群都有矿化相关表达，但直接成骨/矿化差异未显著。这一结果可作为小规模分析型流式的实验起点。')
sub('集中说明的解释边界')
p('聚类、候选选择与统计来自同一表达数据，P/FDR未校正选群过程，不能等同预先指定标签的盲验证或新谱系证明。主患者池的独立性依据公开招募来源与项目内映射；全8来源的患者对应仍未全部解决。')
p('表达推断CNV不是DNA层面真值；留出正常细胞用于阈值校准与质量检查，并非外部准确率验证。未满足输入条件者不可解释为CNV阴性。广义层中的CNV分歧、低复杂度、血管壁及环境RNA标记应与程序信号一起查看。')
p('正常肺相关签名曾作为前期背景检查，当前C7/C9尚未获得完整肺参考分类器投射，因此不使用“肺特异”命名。临床A/B1和同步配对仍按资料完整性处理，不以原发标签或同一患者编号补推。')
p('补充组织队列的外部转移/复发关联未获显著支持（C7/C9的FDR均>0.05）。该结果作为延伸证据保留，不作为本报告发现结论或启动小规模分析型流式的必要条件。')
sub('资料与复核入口')
p('数据和代码范围：Phase2至Phase16提供QC/CNV及来源证据；Phase17为三层联合探索；Phase18为逐群证据登记；Phase20为患者级发现统计。本报告图件直接来自这些输出，新增QC损失分解和代表性CNV展示不改变原始筛选。完整矩阵、逐细胞记录、样本阈值与全部P/FDR保留在项目目录。','note')
p('主要来源：GSE152048原研究，doi:10.1038/s41467-020-20059-6；GSE162454原研究，doi:10.3389/fonc.2021.709210；GSE318840原研究，doi:10.1016/j.jbo.2026.100781；ScPCA官方项目SCPCP000017与SCPCP000023及原始metadata。','note')
page();h('补充分析 取消QC阶段的UMI下限')
p('本次按单因素敏感性设计，取消项目QC中的UMI下限，冻结原样本基因数及线粒体阈值。ScPCA仍要求官方processed成员和官方singlet；既有40样本中原先未进入双细胞评估的新增记录单列。正常参考的独立质量要求及三层候选的基因数、谱系和CNV条件保持不变。')
fig(F/'12_no_UMI_sensitivity.png','图12  ScPCA中有新增细胞的7份样本。蓝色为主分析QC细胞，黄色为取消QC UMI下限后通过其余固定条件的官方singlet。其余5份ScPCA数量不变。',210)
table(['范围','仅因UMI可回收','回收记录的后续状态'],[['ScPCA 12份','5,571个','均为官方singlet；其中34个≥300基因，实测4个通过原始分子候选入口'],['既有40份','995个','原先未获得下游双细胞结果；180个≥300基因，其余条件仍需分开评估'],['全项目52份','6,566条记录','其中214个≥300基因；QC数量上限306,049，不能等同最终合格细胞数']],[98,105,308])
p('明确通过其余冻结QC的ScPCA新增5,571个细胞，使全项目已知合格集合从299,483增加至305,054。剩余995条仅是QC回收上限，不计为新增已确认singlet。6,566条回收记录中6,352条低于当前广义候选层300基因要求，因而取消UMI并不等于按相同数量增加CAF。')
p('ScPCA新增原始分子候选为S610的1个、S611的3个，均不能在缺少新CNV证据时称为CAF。S611另外出现2个符合正常内皮参考分子及质量要求的细胞，为此前参考不足提供了可检验的修复路径。其他低复杂度样本未因取消UMI而自动补足正常参考。')
cnvqa=O/'inputs/no_UMI_CNV611_QA.json'
if cnvqa.exists():
 z=json.loads(cnvqa.read_text(encoding='utf-8'));page();h('补充分析 S611的CNV输入与身份重评估')
 p('保留原S611的3,778个QC细胞，并加入9个无UMI门槛回收的官方singlet，共3,787个。原参考池49个加上2个新内皮参考后为51个，按固定规则分为31个训练参考和20个留出校准细胞；在独立目录重新运行inferCNV、CopyKAT和DecontX。原始主分析对象与全队列群标签未覆盖。')
 table(['检查项','结果'],[['双CNV与DecontX','全部完成且退出码为0'],['参考检查','通过' if z['reference_QC_pass'] else '未通过，保留具体原因'],['双方法低信号细胞',str(z['both_CNV_low_cells'])],['严格候选',str(z['strict_candidates'])],['支持候选',str(z['supported_candidates'])],['支持层原始面板myCAF',str(z['supported_raw_myCAF'])],['严格层原始及校正后myCAF',str(z['strict_raw_corrected_myCAF'])]],[220,291])
 p(f"本次重算保留广义探索候选{z['broad_candidates']}个，均保留参考检查未通过标记；其数量不作为已确认CAF数量。")
 acceptance=pd.read_csv(O/'inputs/no_UMI_CNV611_reference_acceptance.tsv',sep='\t')
 table(['方法','可评估留出细胞','低信号比例','参考检查'],[[r['method'],str(r['holdout_evaluable']),f"{r['low_fraction']:.1%}",'通过' if r['pass'] else '未通过'] for r in acceptance.to_dict('records')], [110,130,130,141])
 p('参考检查要求每种方法至少20个留出细胞可评估，且至少90%判为低信号；输入参考数达标与输出检查通过分别记录。若任一方法未通过，本样本不进入严格层或支持层。','note')
 p('S611恢复了双方法CNV运行条件，但CopyKAT仅16个留出参考可评估，未达到固定的20个要求；这16个均为低信号。inferCNV的20个留出参考中95%为低信号。因此本次未通过的原因是CopyKAT可评估参考数量不足，不是留出细胞出现明确CNV高信号。双方法低信号细胞共3,679个，但在参考检查未通过时不自动升级为CAF。')
 p('本页回答“取消QC UMI下限能否恢复该样本的CNV可评估性”。新增候选属于该样本的独立敏感性结果，不直接并入原54,997个支持层细胞，也不重命名原C7/C9。主统计的三个独立招募队列及患者级结果继续按冻结版本报告。')
else:
 p('S611的独立双方法CNV及去污染重算结果在完成后补入；本版不以输入恢复代替最终身份接受。','note')
def footer(c,doc):
 c.saveState();c.setFont('CN',8);c.setFillColor(colors.HexColor('#65727a'));c.drawString(42,24,'骨肉瘤CAF候选状态完整分析报告  2026年9月21日');c.drawRightString(A4[0]-42,24,str(doc.page));c.restoreState()
out=D/'骨肉瘤CAF_完整主线与分析型流式报告_20260921.pdf'
SimpleDocTemplate(str(out),pagesize=A4,leftMargin=42,rightMargin=42,topMargin=35,bottomMargin=42,title='骨肉瘤CAF候选状态完整分析报告',author='').build(story,onFirstPage=footer,onLaterPages=footer)
(D/'报告正文.md').write_text('\n'.join(md),encoding='utf-8')
print('PDF_CREATED',out,'FIGURES',len(figures),'SHA256',hashlib.sha256(out.read_bytes()).hexdigest())
