from pathlib import Path
import json,html,re,hashlib
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4,A3,landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfgen import canvas
from PIL import Image
from pypdf import PdfReader,PdfWriter
import pypdfium2 as pdfium
B=Path(r'C:\Users\31937\caf_joint_20260919');O=B/'results/phase9_joint_20260919';D=O/'delivery';D.mkdir(exist_ok=True);Q=B/'pdf_qa';Q.mkdir(exist_ok=True)
s=pd.read_csv(O/'tables/joint_cluster_review.tsv',sep='\t');m=pd.read_csv(O/'tables/integration_metrics.tsv',sep='\t');r=pd.read_csv(O/'tables/joint_descriptive_top30.tsv',sep='\t');sample=pd.read_csv(O/'tables/sample_manifest.tsv',sep='\t');p=pd.read_csv(O/'tables/pseudobulk_patient_pairs.tsv',sep='\t')
rank=pd.read_csv(O/'tables/joint_descriptive_all_gene_markers.tsv.gz',sep='\t');save=rank[rank.symbol=='S100A4'];save.to_csv(O/'tables/S100A4_joint_cluster_descriptive.tsv',sep='\t',index=False)
texts=[];story=[]
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
styles={k:ParagraphStyle(k,fontName='STSong-Light',fontSize=size,leading=lead,spaceAfter=7,wordWrap='CJK',keepWithNext=(k=='head')) for k,size,lead in [('body',10,15),('title',20,28),('head',13,19),('cell',8,11)]}
def add(t,kind='body'):
 texts.append(('# ' if kind=='title' else '## ' if kind=='head' else '')+t+'\n');story.append(Paragraph(html.escape(t),styles[kind]))
def table(headers,rows,widths=None):
 texts.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)+'\n')
 data=[[Paragraph(html.escape(str(c)),styles['cell']) for c in row] for row in [headers]+rows]
 tb=Table(data,colWidths=widths or [487/len(headers)]*len(headers),repeatRows=1,hAlign='LEFT');tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5edf3')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,-1),.25,colors.HexColor('#c9d2dc')),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([tb,Spacer(1,8)])
def page():story.append(PageBreak())
add('骨肉瘤 CAF 联合整合分析报告','title');add('2026-09-19 | phase9_joint_20260919 | 原计划报告的联合分析增补')
add('第一批与后续批次已纳入同一联合分析','head')
add('本轮将既有 188 个 CAF 候选细胞、GSE152048 的 2,047 个严格 CAF，以及后续 GEO 的 2,065 个主分析 CAF 候选细胞合并，共 4,300 个细胞、5 个研究、17 条样本记录。按研究来源与测序技术实施 Harmony 批次整合，在同一低维空间重新聚类，得到 J0-J11 共 12 个计算分群。')
add('整合后，同批次近邻比例从 97.49% 降至 73.20%，批次轮廓系数从 0.148 降至 0.017，显示批次分离减弱。分群存在患者主导和参数敏感性，当前结论为联合图谱与候选表达状态，未确立 12 种稳定 CAF 亚型。')
table(['来源','纳入细胞','样本记录','用途'],[[ds,str(int(g.cells.sum())),str(len(g)),'联合整合'] for ds,g in sample.groupby('dataset')],[125,75,75,212])
add('方法：共同特征整合，原始计数用于表达比较','head')
add('保留既有身份审查通过的候选细胞。用稳定 Ensembl ID 对齐，选取所有来源共同测得的 15,391 个基因，不将未测得基因填零。以每个细胞完整原始文库计数为分母计算 log1p(CP10k)；按 study_technology 选择 2,000 个高变基因，排除 MT-/RPL/RPS 高变特征，缩放截断值 10，PCA 30 维。')
add('Harmony 主批次变量为研究×技术，共 6 类，theta=2；15 近邻，Leiden 分辨率 0.5，固定随机种子 20260919。UMAP 用于展示；Harmony 只修正低维坐标。marker、检出率及患者级效应均来自原始 RNA 计数，未使用校正坐标代替基因表达。')
add('表达程序分数为各来源完整测量空间内可用成员基因的平均 logCP10k，逐来源覆盖率另列。部位图将 Lung/lung 合并，并将 Femur 归入 bone；原始部位字段保留。')
page()
add('分群表达特征存在差异，五个群由单个来源患者主导','head')
add('下表列出每群前三个描述性排名基因。排名分数为群内外平均 logCP10k 差值乘以正向检出率差值；这是细胞汇总描述，不是患者级差异显著性。供者记录达标数指该群内至少 20 个细胞的已映射记录，不能据此计算跨队列独立样本量。')
rows=[]
for _,z in s.iterrows():
 cl=int(z.cluster);genes=', '.join(r[r.cluster==cl].symbol.head(3));rows.append(['J'+str(cl),int(z.cells),int(z.known_donor_records_ge20),f'{z.max_known_donor_fraction:.1%}',genes])
table(['分群','细胞数','达标供者记录','最大单供者占比','描述性排名基因'],rows,[42,47,73,83,242])
add('J2、J3、J8、J10、J11 的最大已映射单供者占比超过 70%。J0、J1、J6、J7 有至少三条达标供者记录且未超过该阈值，作为多供者候选状态保留。上述阈值用于复核优先级，不是稳定亚型的验证标准；J0 仍有 69.5% 细胞来自同一供者。')
add('程序层面，J9 显示细胞周期相关表达；J8 显示收缩与基质重塑相关表达；J7 的预设成骨程序相对较高。J4/J5 的炎症、抗原呈递与免疫竞争程序，以及 J6/J10 的肌源性竞争程序需结合后续身份复核。程序热图采用群均值的横向标准化，不能将高分直接解释为细胞来源或功能已被证明。')
add('S100A4 在新的联合 J0 描述性 marker 表中排名第 8；完整表另列出各群检出率与效应。这与旧 H0/H2 排名属于不同聚类结果。S100A4 与客户原文 SR1004 不是同一名称，不能用本结果替代 SR1004 的证据。')
page()
add('批次分离减弱，但分群对校正粒度仍敏感','head')
names={'uncorrected':'未整合','joint_main':'主分析：研究×技术','joint_resolution_0.3':'分辨率 0.3','joint_resolution_0.8':'分辨率 0.8','joint_theta1':'theta=1','joint_dataset_only':'仅按研究校正','joint_sample_batch':'按样本校正','no_FLEX':'去除 FLEX','mapped_without_GSE250015':'去除未知供者及 GSE250015','decontx_available':'DecontX 可用子集'}
rows=[]
for _,z in m.iterrows():rows.append([names[z.analysis],int(z.cells),int(z.clusters),f'{z.same_batch_neighbor_fraction:.1%}','-' if pd.isna(z.get('ARI_vs_main')) else f'{z.ARI_vs_main:.3f}'])
table(['分析条件','细胞数','群数','同批次近邻','相对主分析 ARI'],rows,[187,55,45,95,105])
add('研究粒度校正与主分析较一致（ARI=0.948）；按样本校正改变较大（ARI=0.264）。去除 FLEX 后 ARI=0.629，说明测量技术和样本组成会影响分群。ARI 在共有细胞上计算；每个主群的最佳 Jaccard 对应见附图和完整表，不以不同分析中的数字标签相同认定同一群。')
add('患者内效应可报告，当前分群未达到预设推断性检验门槛','head')
add('在每个研究内按已映射患者汇总原始计数；同患者多个样本合并。仅当该患者在目标群和其余 CAF 各有至少 20 个细胞时，计算配对 log2(CPM+1) 差值。保留 17 项“研究×分群”比较，每项仅 1-2 名配对患者，因此本轮只报告中位效应和正负方向人数，不运行推断性检验；P 与 FDR 字段为空。')
add('预设若达到至少 3 名有效配对患者，使用完整符号翻转检验，并在每个研究×分群的共同基因家族内做 BH 校正；本次没有比较达到该条件。不能把空缺 FDR 表述为“FDR 阴性”。')
add('最早批次中 80 个细胞的供者尚未确认，仅参与联合嵌入和展示；GSE250015 的跨队列供者重叠尚未解决，仅作描述。其余记录也未跨研究合并为总体独立 N。')
page()
add('联合分析更新了分群基础，旧结果需保持各自验证口径','head')
add('旧报告采用发现队列加分开验证；本报告改为所有已接受候选细胞联合整合和重新分群。因此，旧 H0/H2 与新 J0-J11 不存在数字一一对应关系。原 GSE21257 对冻结程序的验证仍属于旧程序结果，不能直接归属于新 J 群；进入本次联合构建的 GEO 队列也不再作为新分群的独立外部验证集。')
add('讨论：当前证据支持联合图谱与候选状态','head')
add('研究、技术、患者和组织部位存在混杂；Harmony 降低批次分离不等同于消除全部技术差异，也不能区分全部患者生物差异。特别是 J11 虽全部为肺来源细胞，但 92.1% 来自单个供者，当前数据不支持将其命名为已验证的肺特异 CAF 亚型。')
add('各阶段候选身份筛选规则具有历史差异，本次保留来源标签并统一整合，并非重新证明全部细胞均为真实 CAF。未将参考论文的分群标签强加于数据。后续论文方法对照应先检查基因覆盖、身份竞争证据和多患者一致性，再确定是否适用。细胞数不代替样本量；目前患者级结果均为探索性描述。')
add('完整交付与可复核记录','head')
add('objects/joint_CAF_integrated.h5ad：4,300×15,391 联合对象，保存共同基因的未改动原始 counts、完整文库计数、来源、细胞标签及主 UMAP。完整来源矩阵仍保留在服务器。另附去除 FLEX、供者映射限制、DecontX 敏感性对象。')
add('tables/：样本清单、分群供者支持、完整描述性 marker 与 top30、患者内基因效应、配对患者清单、程序覆盖、混合指标、逐群敏感性对应和 S100A4 专表。resources/、code/、logs/：来源 SHA256、特征映射、冻结程序、可复现脚本及退出记录。')
add('质量检查：整合与复核进程退出码均为 0；9 个来源文件 SHA256 与整合前一致；联合对象回读通过，4,300 个细胞及 12 群数量一致，原始共同基因 counts 未改动。报告图件附后，所有统计表随结果包交付。')
add('方法参考','head')
add('Korsunsky I et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods (2019). DOI: 10.1038/s41592-019-0619-0. https://doi.org/10.1038/s41592-019-0619-0')
add('Scanpy 官方高变基因接口文档：https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html')
body=Q/'body.pdf'
def footer(c,d):c.setFont('STSong-Light',8);c.setFillColor(colors.HexColor('#657180'));c.drawString(54,27,'骨肉瘤 CAF | 联合整合分析 | 2026-09-19');c.drawRightString(A4[0]-54,27,str(d.page))
SimpleDocTemplate(str(body),pagesize=A4,rightMargin=54,leftMargin=54,topMargin=42,bottomMargin=43,title='骨肉瘤 CAF 联合整合分析报告').build(story,onFirstPage=footer,onLaterPages=footer)
figs=[('01_joint_umap.png','整合前后 UMAP：研究、部位及联合分群'),('02_sample_donor.png','联合 UMAP：样本与供者记录'),('03_program_profiles.png','预设表达程序：群均值标准化'),('04_support_sensitivity.png','研究组成与分群敏感性'),('05_marker_panel.png','原始 RNA marker 检出率与表达')]
f=Q/'figures.pdf';c=canvas.Canvas(str(f),pagesize=landscape(A3));w,h=landscape(A3)
for i,(name,title) in enumerate(figs,1):
 c.setFont('STSong-Light',16);c.drawString(35,h-35,'附图 A'+str(i)+' '+title);im=Image.open(O/'figures'/name);iw,ih=im.size;scale=min((w-70)/iw,(h-95)/ih);c.drawImage(str(O/'figures'/name),(w-iw*scale)/2,(h-ih*scale)/2-10,iw*scale,ih*scale);c.showPage();texts.append(f'![{title}](../figures/{name})\n')
c.save();out=D/'骨肉瘤CAF_联合整合分析报告_20260919.pdf';writer=PdfWriter();writer.append(PdfReader(body));writer.append(PdfReader(f));writer.write(out)
(D/'joint_integration_report.md').write_text('\n'.join(texts),encoding='utf-8');doc=pdfium.PdfDocument(str(out))
for i in range(len(doc)):doc[i].render(scale=1.35).to_pil().save(Q/f'page-{i+1}.png')
text='\n'.join(page.extract_text() for page in PdfReader(out).pages);assert all(t in text for t in ['4,300','15,391','S100A4','SR1004','92.1%'])
print('PDF',out,'pages',len(doc));print('body_pages',len(PdfReader(body).pages));print('text_check_passed')