from pathlib import Path
import json,math,html
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
B=Path(r'C:\Users\31937\CAF_delivery_20261001');I=B/'inputs';out=Path(r'C:\Users\31937\Downloads\骨肉瘤CAF全转录组与空间基质分析_20261002_v2.pdf')
E=I/'phase47_external_candidate_integration_20261001';S=I/'phase46_spatial_functional_20261001';Q=I/'phase48_clinical_programs_delivery_20261001'
REV=Path('C:/Users/31937/CAF_revision_20261002');NS=REV/'inputs/surface';RF=json.loads((REV/'FINAL_REVISION_FACTS.json').read_text(encoding='utf-8'))
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='CNbody',fontName='CN',fontSize=10.2,leading=16.2,spaceAfter=8,wordWrap='CJK'));styles.add(ParagraphStyle(name='CNtitle',fontName='CN',fontSize=23,leading=31,spaceAfter=16,wordWrap='CJK'));styles.add(ParagraphStyle(name='CNhead',fontName='CN',fontSize=16,leading=23,spaceAfter=12,wordWrap='CJK',textColor=colors.HexColor('#234b58')));styles.add(ParagraphStyle(name='CNsmall',fontName='CN',fontSize=8.4,leading=12.2,spaceAfter=5,wordWrap='CJK'));styles.add(ParagraphStyle(name='CNtable',fontName='CN',fontSize=8.7,leading=12.8,wordWrap='CJK'))
story=[]
def p(text,style='CNbody'):story.append(Paragraph(text,styles[style]))
def head(title):p(title,'CNhead')
def page(title):story.append(PageBreak());head(title)
def table(rows,widths=None):
 vals=[[Paragraph(html.escape(str(v)),styles['CNtable']) for v in row] for row in rows];t=Table(vals,colWidths=widths or [480/len(rows[0])]*len(rows[0]),repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9f0f2')),('VALIGN',(0,0),(-1,-1),'TOP'),('FONTNAME',(0,0),(-1,-1),'CN'),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#b4c5cb')),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#dfe6e9'))]));story.append(t);story.append(Spacer(1,10))
def image(path,width=460,maxheight=420):
 from PIL import Image as PIL
 with PIL.open(path) as im:w,h=im.size
 scale=min(width/w,maxheight/h);story.append(Image(str(path),width=w*scale,height=h*scale));story.append(Spacer(1,6))
def fmt(x):
 if pd.isna(x):return 'NA'
 return f'{float(x):.3g}'
facts=json.loads((E/'delivery/REPORT_FACTS.json').read_text());genes=pd.read_csv(E/'tables/main_gene_and_surface_panel_statistics.tsv',sep='\t');sp=pd.read_csv(S/'tables/spatial_patient_statistics.tsv',sep='\t');ext=pd.read_csv(E/'tables/external_patient_statistics.tsv',sep='\t');tf=pd.read_csv(E/'tables/prioritized_TF_regulon_statistics.tsv',sep='\t');source=pd.read_csv(E/'tables/source_projection_statistics.tsv',sep='\t');surface=pd.read_csv(NS/'tables/surface_competing_lineage_statistics.tsv',sep='\t') if (Q/'tables/surface_competing_lineage_statistics.tsv').exists() else None
p('骨肉瘤 CAF 全转录组与空间基质分析','CNtitle');p('F3 F4 功能发现与基质免疫空间证据','CNhead');p('2026年10月2日  |  审计修订结果 v2','CNsmall')
p('本研究围绕客户提出的肿瘤周围胶原与基质保护性环境，优先分析F3和F4，并开放扫描完整可检测转录组。结果将研究重点收敛到F3的选择性基质重塑程序，同时扩展出代谢、干扰素响应、增殖和分泌方向。')
p('<b>主要结论：</b>F3的COL1A1、LOX和MMP14提供患者内一致的基质重塑候选线索；其三基因程序在19名空间患者中与CD8 RNA呈小幅负相关。广义基质主要终点不显著，外部原发患者未整体复现三基因方向，当前证据尚未确认环绕肿瘤的“蛋壳”结构或真实保护作用。')
table([['研究问题','本轮结果'],['F3是否具有基质程序','选择性重塑线索成立；广义胶原程序整体升高未获患者级支持'],['F4是否为胶原屏障群','当前数据不支持；优先关注增殖、分泌及调控程序'],['是否存在免疫空间关联','三基因程序存在小幅CD8 RNA负关联；广义基质主终点阴性'],['能否定义新亚型或直接分选','尚未建立；给出RNA依据明确的候选表面组合与实验验证方案']],[115,365])
p('发现资源包括9个来源、107份组织样本记录。支持层1,224个候选细胞，扩展层8,174个候选细胞；样本记录数和细胞数均不代替独立患者数。完整原始特征并集为63,398，单细胞与空间RNA落实的是全转录组范围。','CNsmall')
page('研究数据与分析方法')
p('保留原始计数与冻结细胞ID，逐来源核对特征ID、HGNC标准名称、别名、全库计数及可测基因覆盖。差异分析使用来源内完整基因范围；共享15,391个基因仅用于既有整合表示。候选分层、CNV及竞争谱系证据沿用已完成的身份审查。')
p('主对比为同患者、同取材部位F3或F4相对其余冻结候选，以及F3对F4、对F2和内部状态对照。原始计数形成患者与部位伪批量，使用对应来源完整文库CP10k归一化，基因效应为log2(CP10k+0.1)差值；部位先聚合，同一来源内每名患者等权。')
p('患者级配对检验采用双侧均值符号翻转；全基因、TF与LR模块在N不超过15时精确枚举，空间与外部模块在N不超过20时精确枚举，超出范围使用登记的Monte Carlo方法。独立临床组采用患者标签置换，配对患者先从独立组移出。候选筛选以预设目标原始P&lt;0.05为名义标准，按基因、程序、通路、TF、临床、通讯、空间及验证各自同类假设家族计算BH FDR。有效患者少于10例的结果按探索性证据呈现。')
p('开放功能分析使用固定2025.1.Hs资源：Hallmark 50、Reactome 1,787、GO BP 7,583和Matrisome 22个定义。完整有效基因排序开展GSEA，原P&lt;0.05的基因按患者平均效应方向开展上调与下调过度代表分析。代表feature在来源内按HGNC对应和固定ID顺序预先选定，再判断可评估性；中位数为零不额外剔除名义基因。GSEA使用1,000次基因集置换；患者平均程序差异另行进行患者级统计。')
p('空间对象包括GSE293065的26份切片、19名患者及53,309个点。程序分数为已测成员log1p(CP10k)均值，CD8 RNA使用CD8A与CD8B，T细胞RNA使用CD3D、CD3E与TRAC。作者IC按原注释聚合为肿瘤、成纤维和T细胞连续有符号分量。')
p('空间关联先在切片内计算秩相关；调整模型将程序和CD8秩分别对肿瘤IC、成纤维IC、完整文库及坐标回归，再计算残差相关。邻域使用第一圈空间邻接，半径为最近正距离中位数的1.05倍。切片依次按事件、部位、患者等权聚合，最终以19名患者检验；未使用spot级P值。')
page('验证与机制分析方法')
p('在验证数据读取前冻结15个功能程序及来自主发现来源的F3、F4签名。GSE237070包含6份样本、4名患者，其中原发患者2名；表达谱按竞争谱系面板注释，主要验证基因及原表面面板从相应注释中留出。参考映射使用患者等权、完整文库归一化的冻结类别均值，并从HVG特征中去除全部冻结程序与签名成员。')
p('内部状态在共享基因表示上选用约2000个HVG和最多50个主成分，经来源与技术Harmony校正后采用15近邻和Leiden分辨率1；分辨率0.5为结构敏感性，随机种子20261001。')
p('患者留出用于检验平均效应方向和单患者贡献；来源投影将已冻结签名映射到各来源的患者效应。F3签名包括74个名义上调和112个下调标准基因；F4签名包括4,320个上调和77个下调标准基因。每个基因沿用实际可评估N，签名覆盖逐来源登记。')
p('TF分析使用CollecTRI有方向、有符号靶关系，构建患者级有符号靶程序响应。配体受体筛选结合角色及膜定位注释，分析CAF候选对T/NK、髓系及其他背景的RNA表达兼容性。来源内患者匹配取材，按预设通讯家族保留完整原P与FDR。')
p('临床比较保留严格骨/肺取材、报告的转移分组和时间尚未确定的历史分组。未知标签保持未知。连续程序以F2、F3、F4的来源患者伪批量值比较；调整模型要求已记录且单一的治疗及技术协变量、满秩设计和足够残差自由度。HC3存在单位杠杆、残差或目标方差数值为零时，保留效应与原审计记录，当前推断P标为NA。')
p('表面候选来自完整基因发现：主来源支持层F3名义上调基因与OmniPath膜受体注释相交。自动得到APP、ATP1A1、BACE1、CA12、FGFR2、LAMP2、MMP14、MRC2、PTK7、SDC1共10个候选，核对竞争谱系及全部45种无序双标组合；候选决策台账完整保留。验证中的映射特征已去除这些发现签名成员。')
p('统计与计数通过独立R核查，包括患者符号翻转、临床标签置换、FDR家族、RNA分数与邻接重建，以及可评估HC3模型的系数、标准误和目标P。完整统计保留零效应、负方向和NA；原始矩阵、脚本、资源版本、哈希及历史失败日志保存在可追溯阶段目录。')
page('内部状态保留结构但未定义新亚型')
p('在冻结类别内，支持层F3的343个细胞得到5个内部状态，扩展层F3的251个细胞得到3个，扩展层F4的5,304个细胞得到17个。分辨率0.5的敏感性分组分别为3、2和14个；这些是计算状态，尚未成为经验证的CAF亚型。')
table([['冻结层与类别','细胞数','内部状态','解释'],['支持层 F3','343','5','来源与患者组成影响内部结构'],['扩展层 F3','251','3','与支持层同编号不等于同一细胞集合'],['扩展层 F4','5,304','17','竞争谱系与CNV冲突保持原标签']],[110,65,70,235])
image(I/'phase42_F3F4_identity_states_20261001/figures/supported_F3_internal_UMAP.png',maxheight=185)
image(I/'phase42_F3F4_identity_states_20261001/figures/expanded_F4_internal_UMAP.png',maxheight=185)
p('内部状态UMAP展示状态、来源与既有身份调用；患者构成另见统计包。少数状态由单一患者或来源主导。分群和全基因差异均保留，候选收敛采用连续程序，避免将来源组成命名为功能亚型。','CNsmall')
page('F3提供选择性基质重塑候选线索')
p('主发现来源支持层F3的COL1A1、LOX和MMP14均名义上调。逐一去除任一患者后，三个基因的平均差异方向仍为正，表明当前方向由多名患者共同贡献。')
rows=[['基因','独立N','平均差异','正向患者','原P','FDR']]
for g in ['COL1A1','LOX','MMP14']:
 r=genes[genes.layer.eq('supported')&genes.symbol.eq(g)].iloc[0];rows.append([g,int(r.N),fmt(r.mean_delta),f'{int(r.positive_patients)}/{int(r.N)}',fmt(r.P),fmt(r.FDR)])
table(rows,[65,50,75,85,95,110]);p('差异单位为log2(CP10k+0.1)变换后的患者平均差值。N=9，证据定位为探索性、名义显著。','CNsmall')
image(I/'phase44_patient_mechanism_20261001/figures/selected_gene_patient_leaveout.png',maxheight=310);p('图1 点为完整患者均值，线段为逐一去除患者后的均值范围，线段不是置信区间。F3与F4均保留来源内生物学N。','CNsmall')
p('三基因合成程序在主发现来源9/9患者为正，平均差异+1.255；来源SCPCP000018的3/3患者方向相同。GSE318840支持层两名患者方向不一致，平均为负。该程序构建于发现结果，发现样本中的合成检验不作为独立验证。')
page('开放发现扩展到代谢与免疫响应')
p('完整排序显示，F3在EMT/基质、糖酵解、干扰素α/γ及部分脂质代谢方向呈模型富集。F4的开放结果更接近增殖、G2M和蛋白分泌。成骨与矿化程序作为独立参考保留，研究范围已扩展到客户原marker之外。')
image(Q/'figures/main_functional_programs.png',maxheight=380);p('图2 固定15程序的患者平均对比差异。F3来自支持层，F4来自扩展层；完整N、原P及FDR见同名统计表。','CNsmall')
p(f'F3的Matrisome Core在修订后的完整排序中NES约{RF["F3_core_matrisome_NES"]:.2f}，患者平均广义程序差异却未达到原P&lt;0.05；广义胶原程序也未整体升高。排名富集与患者均值检验回答不同问题，当前结果支持选择性重塑，不能据排名富集认定所有基质成员一致增加。')
p('干扰素响应与空间CD8 RNA为正关联，糖酵解与CD8 RNA为负关联。免疫响应方向并非单一“免疫抑制”程序，后续优先区分代谢、炎症和选择性基质作用。')
page('F4显示调控程序并保留身份冲突')
p('主来源扩展层F4中，CBFB和TGIF1有符号靶程序在15/15患者上升，均通过对应TF家族FDR校正。F3的HEY2和ZFHX3保留为名义调控线索。')
rows=[['层与类别','程序','N','均值差','原P','FDR','靶数']]
for r in tf.itertuples(index=False):rows.append([r.layer+' '+r.contrast.split('_')[0],r.regulator,int(r.N),fmt(r.mean_delta),fmt(r.P),fmt(r.FDR),int(r.measured_targets)])
table(rows,[95,65,35,70,75,75,65]);image(Q/'figures/TF_reference_network.png',maxheight=260);p('图3 有符号参考靶网络。蓝绿色为正关系，红色为负关系；响应分数来自患者RNA。CBFB仅覆盖3靶、TGIF1覆盖11靶。','CNsmall')
p('F4的RUNX2、ANGPT1和SMAD3存在患者级名义基因上调，但当前签名与胶原屏障群不相符。F4主要来自扩展层，支持层没有独立F4；竞争谱系与CNV身份冲突保留在结论中。')
page('空间主要终点不支持广义基质排斥')
p('19名患者中，广义Matrisome Core与CD8 RNA的调整关联均值为-0.0076，原P=0.371、FDR=0.398，当前数据不支持广义基质程序的CD8排斥关联。未经调整的广义基质与CD8为小幅正关联，调整后趋近于零。')
image(S/'figures/patient_matrix_CD8_effects.png',maxheight=350);p('图4 广义基质主程序的患者级调整秩相关，每点为一名独立患者。','CNsmall')
p('53,309个空间点中53,083个与作者IC匹配，其余点在需要IC的模型中保持缺测。已取得26张低分辨率和26张高分辨率组织图，并核对53,309个RNA坐标；当前未取得独立可验证的肿瘤边界或区域掩膜，“围绕肿瘤边缘的覆盖、连续性及内外CD8定位”这一几何终点不可评估。')
page('选择性重塑与CD8 RNA呈小幅负关联')
p('F3三基因重塑程序在14/19患者中与CD8 RNA负关联；调整后平均ρ为-0.0297，原P=0.00372、FDR=0.00931。邻近空间点的CD8 RNA关联平均ρ为-0.0329，17/19患者为负，原P=0.000565、FDR=0.00169。')
rows=[['程序与终点','平均ρ','负向患者','原P','FDR']]
for name,ep,label in [('F3_SELECTED_REMODELING','CD8_raw','三基因 原始'),('F3_SELECTED_REMODELING','CD8_adjusted','三基因 调整'),('F3_SELECTED_REMODELING','neighbor_CD8_adjusted','三基因 邻域调整'),('NABA_COLLAGENS','CD8_adjusted','胶原程序 调整')]:
 r=sp[sp.program.eq(name)&sp.endpoint.eq(ep)].iloc[0];rows.append([label,fmt(r.mean_effect),f'{int(r.negative_patients)}/{int(r.N)}',fmt(r.P),fmt(r.FDR)])
table(rows,[170,65,85,85,75]);image(Q/'figures/spatial_functional_heatmap.png',maxheight=315);p('图5 全部程序与空间终点的患者平均关联。完整模型与部位敏感性结果均保留。','CNsmall')
p('负关联效应较小。只调整文库及坐标、再分别加入肿瘤或成纤维IC的敏感性分析均保持负方向，原P分别约0.00406、0.00274和0.00449。原始未经调整关联不显著，结果需结合组织成分和归一化解释。')
page('组织图与RNA坐标对齐并显示局部分布')
image(Q/'figures/histology/PT1_tissue_RNA_overlay.png',maxheight=245);image(Q/'figures/histology/MR1_tissue_RNA_overlay.png',maxheight=245);p('图6 预先选择PT1和MR1的组织图叠加。左侧为官方组织图，其余分别叠加三基因重塑RNA、CD8 RNA及作者肿瘤有符号IC；色条为各图连续RNA或IC读出，不代表细胞比例。全部26份叠加图及图像、坐标审核置于统计包。','CNsmall')
p('三基因程序与作者肿瘤IC具有名义正关联，与成纤维IC未形成一致正关联。这表明该空间RNA信号没有被确认仅来源于CAF。观察到的局部共分布与免疫RNA关联是实际数据结果，其细胞来源和肿瘤周边结构仍需组织学定位。')
page('外部验证没有整体确认三基因状态')
p('在GSE237070两个原发患者中，支持层F3参考样细胞的MMP14方向均为正，而COL1A1方向不一致、LOX均为负。F3发现签名也未在两个原发患者中整体同向复现。外部结果保留了关键反例。')
rows=[['基因','N','平均差异','正向患者','原P','FDR']]
for g in ['COL1A1','LOX','MMP14']:
 r=ext[ext.tissue.eq('primary')&ext.layer.eq('supported')&ext.class_id.eq('F3')&ext.kind.eq('gene')&ext.endpoint.eq(g)].iloc[0];rows.append([g,int(r.N),fmt(r.mean_effect),f'{int(r.positive_patients)}/{int(r.N)}',fmt(r.P),fmt(r.FDR)])
table(rows,[80,40,100,90,80,90]);p('外部支持层F3参考样组仅含1和13个细胞，推断N仍为2名患者。N=2的精确符号翻转最小双侧P为0.5；因此以方向、效应与覆盖为主呈现，不能以未通过FDR否定已观察方向，也不能将方向一致称为已验证亚型。','CNsmall')
p('支持与扩展映射分别使用1,260和1,323个留出验证签名后的参考HVG。参考样标签来自表达迁移，外部CNV和双细胞身份未评估。6份样本的4名患者包含原发与癌旁组织，癌旁表达不作为肺转移验证。')
p('来源内患者已经按冻结身份聚合，跨来源完整患者对应仍存在缺口；确认为发现样本重复的外部来源记录已登记并不作为独立复现。各来源结果单列，没有将未知重叠来源合并为独立生物学N。')
page('临床关联尚未形成稳健支持')
p('严格骨/肺取材的191项有效未调整比例检验，以及36项报告或历史分组检验，均没有原P&lt;0.05。连续功能程序完成2,250个注册设计结果，其中105项未调整检验具备有效统计重复，均没有原P&lt;0.05。')
p('连续程序的调整设计因治疗未记录或HC3数学条件不满足而保持不可评估。比例模型经追加数学核查，58个原已报告模型中32个具有有效HC3条件；原计算及完整审计保留。')
p('采用t参考后，F3/F4没有稳健名义比例关联。F1在GSE318840的全部called-cell分母中保留名义信号：N=6，肺相对骨效应约+0.001463，即+0.1463个百分点，原P=0.0303，FDR=1。该次要探索性结果不改变F3/F4主线。')
p('配体受体分析完成2409个可用单基因读出对，443,256项患者级统计，其中130,220项可检验、2849项原P&lt;0.05，无FDR&lt;0.05。575项复合体或未映射读出保持不可评估，候选原P与FDR完整保留。')
p('通讯结果支持后续检测的RNA兼容性候选，当前没有形成实际通讯、T细胞抑制或肿瘤保护的功能证据。CNV-altered/unresolved背景未被称为已确认恶性肿瘤群。')
page('候选检测面板与功能实验方案')
p('优先机制为F3选择性胶原合成、交联与重塑，首选组织检测COL1A1、LOX和MMP14；同时记录糖酵解与干扰素响应，区分基质和炎症环境。F4的CBFB/TGIF1靶程序作为次级调控方向，并在谱系明确的细胞中复核。')
p('按固定规则自动筛得10个名义上调、具有膜注释的候选。原7个之外补齐ATP1A1、BACE1和LAMP2，均进入RNA比较；完整候选决策、竞争谱系与45种双标组合置于修订统计包。膜注释本身不等于适用于活细胞分选的蛋白标志。')
newg=pd.read_csv(NS/'tables/new_surface_discovery_statistics.tsv',sep='\t')
rows=[['候选','N','均值差','原P','FDR']]
for r in newg.sort_values('symbol').itertuples(index=False):rows.append([r.symbol,int(r.N),fmt(r.mean_delta),fmt(r.P),fmt(r.FDR)])
table(rows,[95,45,110,115,115])
p('10个候选均来自主来源9名患者的名义上调，完整全基因家族FDR均未显著。竞争家族扩为10个基因后，MMP14相对壁细胞背景N=10、均值差+3.215、原P=0.00195、FDR=0.0195；相对其余成纤维背景原P=0.00391、FDR=0.0391。')
page('双标组合补全后仍需蛋白验证')
p('主来源支持层F3的318个细胞中，MRC2/MMP14双标RNA共检出137个，即43.1%；MRC2/FGFR2为20个，即6.3%；FGFR2/CA12仅2个，即0.63%，当前不优先采用该稀少组合。这318个细胞来自10名有目标细胞的患者；MRC2/MMP14的患者等权共检出均值为49.8%，中位数47.6%。这些均为描述性RNA检出率。')
p('外部原发F3参考样组的MRC2/MMP14在一名患者1/1细胞共检出，在另一名患者0/13；组合未整体复现。检测优先级据发现和竞争谱系RNA证据制定，蛋白层面仍需验证。')
p('全部45种组合中，APP/MMP14在本发现来源的患者等权RNA共检出率最高，为70.2%，但高检出率不能单独决定谱系特异性或蛋白实验优先级。MRC2/MMP14保留为既有待测组合；当前不将它表述为完整候选比较后已确立的唯一优先方案。')
p('活细胞实验需核实抗体与胞外表位可用性，并记录存活、单细胞及CD45/CD31/EPCAM阴性排除，用壁细胞和恶性成骨对照核查来源。RNA双检出率不替代蛋白门纯度，不能直接将双阳性细胞命名为F3。')
if surface is not None:
 ss=surface[surface.dataset.eq('SCPCP000017')&surface.layer.eq('supported')&surface.class_id.eq('F3')&surface.gene.isin(['MRC2','MMP14','FGFR2'])];ss.to_csv(B/'surface_primary_context_reference_v2.tsv',sep='\t',index=False)
 table([['实验步骤','回答的问题'],['共染与胶原成像','基质程序由哪类细胞产生；是否位于肿瘤周边；CD8位于何处'],['双阳性蛋白检测与谱系对照','候选表面组合能否富集目标程序并排除壁细胞与恶性成骨背景'],['三维基质或组织共培养','基质是否减少T细胞进入；阻断LOX/MMP14等是否改变进入和肿瘤存活']],[155,325])
p('功能实验需使用独立患者来源重复，分别测量基质沉积/交联、T细胞进入及肿瘤存活，才能进一步区分结构作用、免疫作用与真实保护效应。以上为后续实验设计，本轮未执行湿实验。')
page('讨论与证据边界')
p('本轮把候选重点从“所有胶原升高的F3/F4新亚群”收敛到F3的选择性重塑连续程序。三个核心基因的发现方向和患者留出稳定性、以及19名患者中的小幅空间负关联，构成继续进行组织与功能定位的依据。广义基质主终点阴性、来源异质性和外部部分失败同时限定结论。')
p('F3主要发现N=9，名义基因信号不能称为校正后显著；空间患者N=19的部分程序关联通过FDR，但ρ约-0.03，统计证据强度不等于大效应或临床预测能力。空间点数、外部细胞数及技术重复均没有扩大生物学N。')
p('RNA点混合多种细胞，IC为有符号分量而非细胞比例；基质和CD8程序共享总文库分母，调整与邻域分析用于检验结构化关联，但未构成组织内空间随机化因果实验。缺少独立可验证的肿瘤边界标签时，不能确认环形结构、纤维方向、肿瘤内外免疫定位或真实保护。')
p('F4身份冲突、内部状态来源主导、参考样迁移不完全及跨来源患者对应不闭合，意味着目前不能宣称一个可复现CAF新亚型、CAF专属机制或可直接临床使用的分选门。全转录组发现和RNA表面候选分析均未测量DNA全基因组或蛋白活性。')
p('以原P&lt;0.05进入候选筛选，同时保留FDR、效应与患者一致性，使名义线索可继续研究。缺有效重复、真实全零、缺少协变量或模型数学条件不满足时保留NA，不补造推断结果。')
page('补充统计与可复核资源')
p('统计包提供全基因、基因集、通路、患者稳定性、TF、通讯、临床比例与连续程序、空间关联、外部迁移及表面候选完整结果，并附脚本、资源版本、完成标记、退出码、输入引用及SHA256清单。大体积原始计数对象保留于sgs10项目，统计包记录其路径与哈希。')
table([['阶段','结果与当前入口'],['41 全基因覆盖','原始计数对齐、特征与来源覆盖、患者台账'],['42 类别与内部状态','支持/扩展对应、身份与UMAP'],['43 全转录组发现','完整基因与程序统计、全排序通路及ORA'],['44和45 机制','患者稳定性、TF、有条件临床与RNA配体受体'],['46 空间','全26切片图、组织图对齐、患者关联与模型敏感性'],['47 外部与整合','留出签名迁移、来源投影、候选与来源重叠审计'],['48与49 临床与表面','连续程序、当前HC3数学审查表、10个表面候选及45种RNA组合']],[125,355])
p('资源来源：NCBI GEO GSE293065与GSE237070；MSigDB 2025.1.Hs；HGNC完整基因名称集；OmniPath与CollecTRI固定下载。官方入口：','CNsmall')
for url in ['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293065','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237070','https://www.gsea-msigdb.org/gsea/msigdb/','https://www.genenames.org/download/archive/','https://omnipathdb.org/']:
 p(f'<link href="{url}" color="#234b58">{url}</link>','CNsmall')
p('当前通路、上下调汇总和表面候选使用阶段49；当前比例模型仍使用阶段48的D_fraction_HC3_current_statistics.tsv及P_t/FDR_t。修订包补入阶段44的严格部位、报告分组和模型敏感性子目录，CURRENT_RESULT_INDEX.tsv逐项列明当前入口。状态快照、原失败日志和不可评估条目均不替代当前验收记录。','CNsmall')

page('修订改变通路细节并补齐候选与统计交付')
p('本版统一通路代表feature的选择顺序，并按注册的患者平均效应决定名义基因方向。全基因原始效应、P值和患者级生物学N保持原定义；本次修改对应下游输入及解释口径。')
p(f'196个对比与通路库组合共925,316条ORA记录中，{RF["ORA_nominal_status_changes"]}条名义状态及{RF["ORA_FDR_status_changes"]}条FDR状态变化。156个输入变化的GSEA任务重新计算，40个输入相同的任务按哈希和内容一致性复用。')
table([['F4 主来源 ORA','名义条目 原→新','FDR条目 原→新'],['Hallmark','21→20','17→14'],['Reactome','480→403','208→159'],['GO BP','2041→1769','773→596'],['Matrisome','13→13','12→11']],[170,155,155])
p('表中每条为通路与方向的组合，不是去重生物学通路。原始结果和变更明细全部保留。F3三基因的患者级名义差异、CD8 RNA小幅负关联、广义空间主终点阴性，以及外部未整体复现的事实不受本次修改影响。')
p('表面候选补齐为10个及45种无序双标组合；不再用不完整的候选全集确定唯一实验优先级。最终临床修订的42个文件已补入新版统计包，并按当前结果依赖逐项验收。')
p('缺少独立肿瘤边界、真实功能实验、完整跨来源患者对应和蛋白门验证的项目，仍保持原有证据状态；本次计算修订没有将其转为已验证结论。')

def footer(c,d):
 c.setFont('CN',8);c.setFillColor(colors.HexColor('#60737b'));c.drawString(48,27,'骨肉瘤 CAF 全转录组与空间基质分析  2026-10-02');c.drawRightString(A4[0]-48,27,str(d.page))
SimpleDocTemplate(str(out),pagesize=A4,rightMargin=48,leftMargin=48,topMargin=44,bottomMargin=47,title='骨肉瘤CAF全转录组与空间基质分析',author='CAF项目分析').build(story,onFirstPage=footer,onLaterPages=footer)
print('REPORT_CREATED',out)