from pathlib import Path
import html,json,shutil,hashlib,zipfile
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from pypdf import PdfReader
import pypdfium2 as pdfium
B=Path(r'C:\Users\31937\caf_method_audit_20260919');O=B/'results/phase9_method_audit_20260919';D=O/'delivery';D.mkdir(exist_ok=True);Q=B/'pdf_qa';Q.mkdir(exist_ok=True)
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'));styles={k:ParagraphStyle(k,fontName='STSong-Light',fontSize=sz,leading=ld,spaceAfter=7,wordWrap='CJK',keepWithNext=k in ['title','head']) for k,sz,ld in [('title',19,27),('head',13,19),('body',10,15),('cell',8,11)]}
story=[];md=[]
def add(t,k='body'):story.append(Paragraph(html.escape(t),styles[k]));md.append(('# ' if k=='title' else '## ' if k=='head' else '')+t+'\n')
def table(headers,rows,widths):
 md.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)+'\n');t=Table([[Paragraph(html.escape(str(x)),styles['cell']) for x in r] for r in [headers]+rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4edf3')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,-1),.25,colors.HexColor('#cbd5df')),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,8)])
def page():story.append(PageBreak())
add('骨肉瘤 CAF 分群方法与阈值审计','title');add('2026-09-19 | 审计对象：phase9_joint_20260919 | 本轮未修改正式分析标签或参数')
add('联合整合适用，当前阈值尚不能作为已验证的亚型标准','head')
add('结论：保留多批次联合分析与 Harmony 的总体框架；当前 12 群可作为探索性计算分群。首要问题不是强行减少群数，而是各来源身份准入规则不一致、部分排除门槛会明显改变可发现状态，以及患者与技术混杂。参考论文可用于多基因程序对照，暂不适合作为必须分成四群的目标。')
table(['审计项目','判断','本轮证据与处理建议'],[
['联合分析与批次校正','适用','共同测量基因对齐、保留原始 counts 的方案可保留；不以批次完全混合为优化目标。'],
['Harmony theta=2','合理工作基线','已有 theta=1、研究级和样本级敏感性；未证明 theta=2 是最优值，不建议仅为 UMAP 更混合而加大。'],
['Leiden 0.5 / 15近邻','局部稳定','补做 30 个诊断组合；0.3-0.5 区间整体 ARI 0.958-1.000，11-12 群。'],
['HVG 2000 / PC 30','常用起点，尚未充分验证','本轮未重算不同 HVG 数或 PC 数；不能称为最优配置。'],
['身份准入和成骨排除','优先复核','旧新批次规则不同；成骨 0.2% 改为 0.5% 会增加 213 个后续 GEO 候选。'],
['20细胞 / 70%占比','仅操作门槛','可作质量下限和风险提示；不能据此证明统计充分或亚型稳定。']],[102,96,289])
add('软件实现需在 Methods 中补充','head')
add('实际使用 Scanpy 1.11.5、harmonypy 2.0.0。显式 theta=2、最大迭代 30；未指定 lamb 时该版本采用动态 lambda（alpha=0.2），不是固定 lambda=1；Harmony 内部软聚类数为 100，与最终 Leiden 的 12 群不同。主整合在 6 次迭代后收敛。此信息已从安装源码和日志核实。')
page()
add('身份阈值敏感性主要集中于程序排除和最低阳性基因数','head')
add('下表仅针对后续 GEO 的已保留输入，固定其他条件，重建原规则后与 2,065 个主候选逐细胞一致。所有增加数均为“如果改变规则”的诊断结果，未重新认定为 CAF；不能外推为全部原始样本的召回率。')
table(['规则变化','候选数','相对2065','审计判断'],[
['检出基因≥300 / 500','2070 / 2065','+5 / 基线','500 下限可保留；放宽收益很小。'],
['检出基因≥800 / 1000','1896 / 1658','-169 / -407','直接提高会明显删去低复杂度细胞；需按技术审查。'],
['Fib基因≥1 / 2 / 3','2249 / 2065 / 1538','+184 / 0 / -527','≥2 是弱准入证据；数量敏感，不能单独决定身份。'],
['群CNV风险上限10%-50%','均2065','0','在其余严格门槛固定时不敏感，不等于已验证最佳阈值。'],
['竞争/混合程序0.1%','1447','-618','收紧会大量剔除细胞。'],
['竞争/混合程序0.5%','2491','+426','放宽会大量新增；必须区分污染、双细胞、真实重叠状态。'],
['仅成骨排除0.1%','1965','-100','成骨相关状态对准入规则敏感。'],
['仅成骨排除0.5% / 1%','2278 / 2369','+213 / +304','建议保留单独敏感性集合，不能直接全部放回主集。']],[145,92,78,172])
add('0.2% 指程序计数占完整文库计数的百分比，代码比较值为 0.2，不是 20%。成骨冲突还要求至少 3 个该程序基因被检出；本诊断同时改变原始与 DecontX 校正后的成骨阈值。全程序诊断同时改变竞争与混合程序的百分比门槛，破骨专用门槛保持其一半；该效应不能全部归因于成骨程序。')
add('统一准入定义比继续微调聚类分辨率更优先','head')
add('联合流程直接读取三个阶段已接受集合，并未统一重做身份准入。旧 GSE152048 的严格集合构建代码未采用后续 GEO 同样的显式成骨冲突和 DecontX 后准入门槛；最早 188 个细胞也缺少本轮匹配的去环境 RNA 检查。因此当前联合图谱存在选择规则差异，Harmony 不能修复这一差异。')
page()
add('论文四类适合程序对照，不能直接套用分类阈值','head')
add('核对了用户提供论文全文、出版社正文及补充文件 5/6。作者提供 Seurat 分析框架及 CD36、MYH11、C7、FNDC1 的表达特征，但所查方法文字和补充图注未给出足以复现其亚群的完整 QC、分辨率、PC、批次和分类阈值。原单细胞数据声明为按需提供，本项目尚未取得参考对象。[1]')
add('论文的 MYH11 类包含 RGS5、MCAM 等血管壁相关基因；本项目严格规则会排除较强 Mural 程序，两者身份定义不一致。要检验该类是否存在，需另行保留已排除的 Mural/SMC 候选做并列证据核查，不能仅凭严格 CAF 集合的缺失下结论。')
add('CD36 与 MYH11 的单基因阳性数量高度依赖检测门槛','head')
d=pd.read_csv(O/'tables/reference_program_detectability.tsv',sep='\t');g=d.groupby(['program','count_threshold'])[['anchor_detected_cells','anchor_and_ge2_genes_cells']].sum();rows=[]
for name in ['CD36','MYH11','C7','FNDC1']:
 x=g.loc[(name,1)];y=g.loc[(name,2)];rows.append([name,int(x.anchor_detected_cells),int(y.anchor_detected_cells),int(x.anchor_and_ge2_genes_cells),int(y.anchor_and_ge2_genes_cells)])
table(['程序锚定基因','计数≥1','计数≥2','锚定+≥2程序基因：≥1','锚定+≥2程序基因：≥2'],rows,[85,65,65,136,136])
add('分母均为 4,300 个联合候选细胞。后两列要求锚定基因阳性，且包括锚定基因在内至少两个程序基因达到相同计数门槛；不是已验证的亚型分配规则。原始计数≥1 的 CD36 阳性为 84 个（1.95%），≥2 仅 16 个（0.37%），不宜凭单个转录本直接命名 CD36+ CAF 亚群。')
add('四套面板在大多数来源完整测得；FLEX 来源的 C7 六基因面板缺 IGF2，因此只能评估五个成员。参考程序应在明确测量覆盖后按研究内比较，并允许混合和未匹配状态。本轮只评估可检出性，未输出论文四类的强制标签，也未重复其疗效关联结论。')
add('建议的应用方式','head')
add('使用论文多基因程序作为连续注释证据，结合原始/去环境 RNA 表达、竞争谱系、技术覆盖与多患者一致性；如需硬分类，应在独立标注或留出患者数据上校准阈值和拒绝分配规则。当前没有足够依据宣布固定“两个基因阳性”或某个分数就是该亚型。')
page()
add('整体分群较稳定，但 J8/J9 边界与供者主导需要分别判断','head')
add('在冻结主 Harmony 坐标上，检查近邻数 10/15/30、Leiden 分辨率 0.1/0.2/0.3/0.5/0.8、两个随机种子，共 30 组。分辨率 0.3-0.5 区间为 11-12 群、ARI 0.958-1.000；分辨率 0.1 为 6-7 群，0.8 为 13-17 群。没有证据支持为了迎合论文直接指定四群。')
add('局部对应揭示全局 ARI 的盲区：0.3-0.5 区间 J9 最佳 Jaccard 中位数约 0.303、最小 0.265；J8 中位数约 0.723、最小 0.682。J9 更适合暂列细胞周期状态，J8/J9 的分界需专门复核。J10/J11 在邻域参数下很稳定，但分别约 97.2%/92.1% 来自同一已映射供者；数值稳定不能等同于独立生物学复现。')
add('批次校正评价需要同时保留生物信息','head')
add('已有敏感性中，按样本校正与主分析 ARI=0.264、去除 FLEX 后 ARI=0.629、DecontX 可用子集 ARI=0.813。样本校正改变较大不单独证明过校正，但提示“混得更匀”不能作为选参依据。HVG 数、PC 数、文库复杂度、技术测量范围和供者留出仍需分别评价。[2,3]')
add('20 个细胞可作操作下限，降低它不能补足患者重复','head')
table(['每侧最低细胞数','可描述的研究×分群比较','最多配对供者','达到≥3供者的比较'],[['5','27','3','1'],['10','24','3','1'],['20','17','2','0'],['30','13','2','0'],['50','7','1','0']],[113,156,95,123])
add('降到 10 个细胞仅使 GSE152048 的 J9 达到 3 个配对供者，其双侧精确符号翻转检验理论最小 P=2/2³=0.25。达到“≥3”仅满足计算门槛，不意味着具备显著性检验能力。不要为了得到 P 值降低身份或细胞门槛；本轮没有新增推断性检验。')
add('单供者占比 70% 仅作风险标记：门槛设为 60%、70%、80%、90% 时，分别标记 7、5、4、3 群。应报告连续占比、每供者细胞数和跨研究支持，不以是否低于 70% 宣布稳定亚型。')
page()
add('建议顺序：先统一身份定义，再冻结生物分群','head')
table(['优先级','建议','验收依据'],[
['1','统一不同批次的准入证据表','同一细胞保留原标签与统一规则判定；并列严格CAF、成骨冲突候选和Mural/SMC候选，记录剔除原因。'],
['2','校准程序阈值的用途','0.2% 暂留为严格集合规则；新增敏感性集合单独展示。依据正常参考、污染、CNV和谱系证据判断，不以新增数量决定。'],
['3','保留当前整合作为比较基线','在统一集合上再比较 HVG/PC 与供者平衡、留出稳定性；按生物保留和批次处理共同选参。'],
['4','先做连续状态注释','J9 按细胞周期状态复核；成骨、收缩、基质、炎症和论文程序并列展示，不预设四类或十二类必然成立。'],
['5','冻结命名与外部验证口径','仅对有身份证据、跨患者一致且不依赖单一参数的状态命名；构建中使用过的队列不能再当独立验证。']],[47,152,288])
add('CNV 门槛有内部校准依据，尚不是真实身份的金标准','head')
add('inferCNV 信号阈值采用 max(0.05, 留出正常细胞99%分位数, 中位数+5×max(MAD,0.005))；升高还要求至少两条染色体中位偏移≥0.1，或一条≥0.25。CopyKAT 使用其另一路判定，并检查正常参考表现。该做法比单一全局阈值有依据，但同一样本内留出正常细胞用于校准，不能再视为完全独立的阈值验证；双低信号也不证明一定为非恶性。FLEX 测量与弱 CNV 的误判边界仍需保留。')
add('审计范围、验收与尚未验证的事项','head')
add('本轮完成代码/日志核对、30 组冻结嵌入诊断、后续 GEO 身份规则逐项敏感性、成骨单项敏感性、参考程序覆盖和检出性、患者数及供者占比门槛核查。初次审计因 raw_ 字段前缀不匹配退出，修正后基线与原准入标签完全一致，最终主审计及追加诊断退出码均为0；9个来源文件 SHA256 与冻结记录一致。正式 4,300 细胞对象与标签未被重写。')
add('尚未完成：统一全部批次身份规则、HVG/PC 全流程重拟合、去留一供者重拟合、留出亚型分类校准、独立肿瘤身份真值验证。上述事项属于建议的后续实施，不能计入本次已验证结果。')
add('参考来源','head')
add('[1] Zeng et al. Molecular Cancer 25:14. DOI: 10.1186/s12943-025-02516-2；用户提供全文及出版社补充材料5/6。https://pmc.ncbi.nlm.nih.gov/articles/PMC12829093/')
add('[2] Harmony 官方算法/参数文档：https://portals.broadinstitute.org/harmony/reference/RunHarmony.html；具体 Python 默认值以本次存档的 harmonypy 2.0.0 安装源码为准。')
add('[3] Luecken et al. Benchmarking atlas-level data integration in single-cell genomics. https://doi.org/10.1038/s41592-021-01336-8')
def footer(c,d):c.setFont('STSong-Light',8);c.setFillColor(colors.HexColor('#657180'));c.drawString(54,26,'CAF 方法与阈值审计 | 2026-09-19');c.drawRightString(A4[0]-54,26,str(d.page))
out=D/'骨肉瘤CAF_分群方法与阈值审计_20260919.pdf';SimpleDocTemplate(str(out),pagesize=A4,rightMargin=54,leftMargin=54,topMargin=42,bottomMargin=43,title='骨肉瘤CAF分群方法与阈值审计').build(story,onFirstPage=footer,onLaterPages=footer)
(D/'method_threshold_audit.md').write_text('\n'.join(md),encoding='utf-8');pdf=pdfium.PdfDocument(str(out))
for i in range(len(pdf)):pdf[i].render(scale=1.4).to_pil().save(Q/f'page-{i+1}.png')
print('PAGES',len(pdf));text='\n'.join(p.extract_text() for p in PdfReader(out).pages);assert all(x in text for x in ['2065','2278','0.958','0.25','harmonypy'])
print('TEXT_CHECK_PASS')