from pathlib import Path
import json,shutil,hashlib,re
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pymupdf
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'full_report_integration';W.mkdir(exist_ok=True);D=B/'review_bundle_final';target=Path('C:/Users/31937/Downloads/骨肺CAF亚群_图文分析与Lung04推荐报告.pdf');backup=W/'original_customer_report.pdf'
if not backup.exists():shutil.copyfile(target,backup)
baseline=pymupdf.open(backup);assert len(baseline)==11
for name,file in [('Song','simsun.ttc'),('Hei','simhei.ttf')]:pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
styles={k:ParagraphStyle(k,fontName='Hei' if k in ['h','title'] else 'Song',fontSize=size,leading=lead,spaceAfter=8,wordWrap='CJK') for k,size,lead in [('body',10.5,17),('small',9,14),('h',15,22),('title',19,27)]}
story=[]
def p(t,k='body'):story.append(Paragraph(t,styles[k]))
def h(t):p(t,'h')
def fig(path,width=490):
 from PIL import Image as PI
 w,z=PI.open(path).size;story.append(Image(str(path),width=width,height=width*z/w))
def table(rows,widths):
 t=Table([[Paragraph(str(v),styles['small']) for v in row] for row in rows],colWidths=widths,repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff3')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#587080')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#d8dee4')),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6)]));story.append(t)
bg=pd.read_csv(D/'tables/normal_IPF_background_comparison.tsv',sep='\t');nodes=['Lung_04','Bone_07','Bone_07_s1','Lung_04_s2'];programs=['mCAF_matrix','collagen_biosynthesis','collagen_crosslinking','fibrillar_I_III_V','collagen_processing_nonCOL'];labs=['Matrix program','Collagen synthesis','Crosslinking','Fibrillar collagen','Non-COL processing'];sel=bg[(bg.query_min_cells==5)&(bg.reference_min_cells==20)&(bg.condition=='Healthy')&bg.reference_type.isin(['Alveolar fibroblasts','Adventitial fibroblasts'])&bg.program.isin(programs)];scale=sel.difference.abs().max()
plt.rcParams.update({'font.size':10,'savefig.dpi':180});f,axs=plt.subplots(1,2,figsize=(10,4.6),sharey=True)
for ax,typ in zip(axs,['Alveolar fibroblasts','Adventitial fibroblasts']):
 a=sel[sel.reference_type==typ].pivot(index='node_id',columns='program',values='difference').reindex(index=nodes,columns=programs).to_numpy();im=ax.imshow(a,vmin=-scale,vmax=scale,cmap='RdBu_r',aspect='auto');ax.set_title(typ);ax.set_yticks(range(4),nodes);ax.set_xticks(range(5),labs,rotation=55,ha='right')
 for i in range(4):
  for j in range(5):ax.text(j,i,f'{a[i,j]:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(a[i,j])>scale*.6 else 'black')
f.subplots_adjust(bottom=.4,right=.89,wspace=.08);f.colorbar(im,cax=f.add_axes([.91,.4,.02,.45]),label='Candidate - healthy score');f.savefig(W/'healthy_background_main.png',bbox_inches='tight');plt.close(f)
notes=json.loads((B/'REPORT_INTERPRETATION.json').read_text(encoding='utf-8'))
h('十、HLCA外部参照：补充候选的细胞身份解释')
p('在原骨肺队列完成亚群筛选和患者级效应评估后，本节进一步回答两个问题：这些候选更接近哪类细胞，以及其基质表达是否超过正常肺组织原有背景。重点分析Lung_04、Bone_07、Bone_07_s1和Lung_04_s2；其余信号群仍按原报告分层保留。')
p('HLCA（Human Lung Cell Atlas，人类肺细胞图谱）整合了多个研究的肺组织单细胞数据。本研究使用其健康核心图谱作为已知细胞身份参照，并从完整图谱提取健康肺及肿瘤邻近肺成纤维细胞进行表达比较。外部参考未并入原6,423个主分析细胞，不改变原分群、UMAP或有效配对患者数。')
fig(B/'figures/identity.png')
p('图12｜HLCA健康核心图谱身份映射。每位患者先计算候选细胞的身份比例，再对患者等权平均；候选至少5细胞。颜色对应成纤维细胞、周细胞、平滑肌和其他谱系。','small')
p('Lung_04的7位可评估患者均来自NCI，患者平均成纤维细胞样比例为94.7%；Lung_04_s2的2位患者为100.0%。Bone_07和Bone_07_s1各4位可评估患者均来自SCPCP000017，其成纤维细胞样比例为82.2%和75.6%，平滑肌样比例为17.8%和24.4%。')
p('保持参考参数冻结，由200轮续训至500轮后，整体身份大类一致率为97.4%；四项候选分别为96.1%、96.9%、91.4%和100.0%。这些结果支持以成纤维细胞样为主的身份解释，并保留骨侧候选的收缩/平滑肌样混合信息。')
story.append(PageBreak());h('十一、正常肺背景：候选的基质程序高于参考水平')
p('肺组织本身含有表达胶原和基质基因的成纤维细胞，因此“检出基质基因”需要结合正常背景解释。本节分别采用肺泡和外膜成纤维细胞作参照，检验主分析发现的基质相关方向。')
fig(W/'healthy_background_main.png')
p('图13｜候选患者功能分数中位数减去健康参考的研究等权平均分数。单位为平均log1p(CP10k)，红色为候选较高，蓝色为参考较高。候选每患者至少5细胞；参考每供者每类至少20细胞。','small')
p('四项候选的基质、纤维型胶原及交联程序均高于健康肺泡和外膜成纤维细胞背景。Lung_04相对健康肺泡背景的分数差分别为+1.704、+1.772和+0.610；其非COL胶原加工程序也偏高，但在正常外膜成纤维细胞中已有一定基础表达，候选的研究等权中位百分位为0.759。')
p('改用肿瘤邻近肺背景后，Lung_04的基质和交联方向仍较高。将外部比较的候选细胞门槛提高至每患者10细胞、参考提高至每供者每类50细胞后，Lung_04仍有7位可评估患者，主要基质方向保持一致。这里无需同患者对照，不能替代原主分析min10的5位有效配对患者。')
table([['候选','原主分析有效配对N（min5）','外部比较N（候选≥5细胞）'],['Lung_04',7,7],['Bone_07',4,4],['Bone_07_s1',3,4],['Lung_04_s2',1,2]],[130,180,180])
p('表6｜两类患者数的口径。原主分析要求同患者目标与对照均满足细胞门槛；外部比较只要求候选自身满足门槛。父子群有包含关系，患者数及细胞数不能相加。','small')
story.append(PageBreak());h('十二、外部参照的方法及与主分析的综合解释')
p('身份映射：复用HLCA v1.1预训练模型及健康核心584,944细胞的参考空间。按来源设置查询批次，以2000基因原始计数适配模型，使用50近邻加权投票；主映射200轮，并进行参考参数冻结下追加300轮的敏感性复核。NCI匹配1,631/2,000个模型特征，缺失特征仅作模型补零。')
p('正常表达背景：从HLCA完整图谱的原始计数提取基质细胞，按研究、供者、条件和细胞类型固定种子抽样，每层最多100细胞。健康肺覆盖140位供者、17项研究，肿瘤邻近肺覆盖33位供者、5项研究；具体类型和细胞门槛下的有效供者数另行统计。参考供者采用study::donor标识。')
p('表达与汇总：联合比较使用378个共同唯一基因，各程序覆盖率至少80%。按全库原始计数归一为CP10k，log1p后对程序基因求均值；先按患者/供者汇总，健康参考再按研究等权。候选至少5细胞、参考每类至少20细胞为主要设置，10/50细胞为敏感性设置。跨数据集比较报告分数差和百分位，不给出跨数据集P值。')
h('对候选推荐的补充')
p('Lung_04继续优先推荐为ECM/胶原生成相关CAF候选。推荐基础仍是原队列7位有效配对患者的基质、胶原交联、纤维型胶原及FNDC1方向一致性；HLCA补充了成纤维细胞样身份和超过正常肺背景的表达证据。Bone_07作为骨侧基质与收缩相关补充候选，其子群及Lung_04_s2保留为细分探索对象。')
p('Lung_04较健康肺泡背景的典型myCAF程序差值仅+0.036，成骨程序为+0.800，矿化程序为+0.072；Bone_07_s1的矿化差值为-0.039。因而正文采用“ECM/胶原生成相关”定位，保留成骨相关转录特征，后续围绕胶原分泌与沉积开展实验验证。')
h('解释与核验')
p('外部图谱提供身份与功能参照，不能计作本项目新增的独立候选患者验证。不同数据集的平台、测序深度和既有过滤影响绝对表达分数；骨侧比较属于跨器官参照。映射投票不是校准后的身份概率，应与原队列标志物及患者证据共同解释。')
p('本次补充分析47项核验通过，涵盖细胞名单、原始计数、程序分数重算、患者汇总、单基因复核及冻结参数检查。原主分析的统计口径与结果保持不变。探索性疾病状态参照及关键敏感性结果列入补充材料S2–S3。')
p('[10] Sikkema L, et al. An integrated cell atlas of the lung in health and disease. Nature Medicine, 2023. DOI: 10.1038/s41591-023-02327-2。HLCA模型：Zenodo 7599104；完整表达图谱：CELLxGENE数据集9f222629-9e39-47d0-b83f-e08d610c7479。','small')
story.append(PageBreak());p('补充材料','title');h('S1｜Lung_04部分标志物的患者配对结果')
p('以下表格沿用原完整版补充表S1，评估来源为NCI，采用原主分析min5配对口径。')
table([['基因','有效N','正向患者','平均效应','原始P','FDR'],['POSTN','7','4/7','-0.163982','0.625000','0.772059'],['ACTA2','7','2/7','-0.023889','0.859375','0.984375'],['TAGLN','7','2/7','+0.029591','0.921875','1.000000'],['LOXL2','7','1/7','-0.093160','0.546875','0.710374']],[72,50,70,100,100,98])
story.append(Spacer(1,20));h('S2–S3的用途：探索性纤维化状态参照')
p('IPF（Idiopathic Pulmonary Fibrosis，特发性肺纤维化）具有肺组织瘢痕化和细胞外基质沉积特征。HLCA包含IPF供者，本补充利用这些成纤维细胞作为纤维化相关状态的外部功能参照，以判断候选的胶原表达是否伴随更广泛的相关程序。该探索性比较不作为CAF筛选或Lung_04推荐的门槛。')
p('原HLCA论文Fig.6h报告了CCL2、COL1A1、CTHRC1、MMP19、SERPINE1和HIF1A等IPF相关基因。S2使用该六基因均分比较候选与参考，S3进一步去除COL1A1并展示单基因证据，判断评分是否主要由胶原表达推动。IPF参考是外部疾病状态，不代表本项目患者患有IPF。')
h('补充比较方法')
p('完整抽样参考共24,983个细胞，包括健康肺、肿瘤邻近肺及IPF基质细胞；IPF覆盖52位供者、4项研究。候选与参考采用正文相同的归一化及供者汇总方法。健康与IPF检验限定同研究、同类型且各至少3位供者，采用双侧Mann-Whitney检验，BH-FDR在同研究和同细胞门槛下的程序/类型对比家族内计算。各候选与外部参考之间只作描述比较。')
p('完整原始P值、FDR、患者汇总和数值核验表保存在服务器Phase122分析目录；以下图表保留影响解释的主要结果。')
story.append(PageBreak());h('S2｜IPF六基因功能参照与候选表现');fig(B/'figures/ipf.png')
for text in notes['ipf']:p(text)
p('补充图S1｜HLCA-IPF六基因均分，单位为平均log1p(CP10k)。点为患者/供者，黑线为中位数；绿色为健康肺泡成纤维细胞，橙色为IPF肺泡成纤维细胞，蓝色为本项目候选。参考每位供者每类至少20细胞，候选每位患者至少5细胞；括号标出各组供者/患者数。该图用于外部功能描述比较。','small')
story.append(PageBreak());h('S3｜单基因及去除COL1A1的敏感性分析');fig(B/'figures/markers.png')
for text in notes['gene_evidence']:p(text)
nc=pd.read_csv(D/'tables/IPF_without_COL1A1_background.tsv',sep='\t');rows=[['候选','去COL1A1后健康百分位','去COL1A1后IPF百分位']]
for n in nodes:
 a=nc[nc.node_id==n].set_index('reference_group');rows.append([n,f"{a.loc['Healthy | Alveolar fibroblasts','median_reference_percentile']:.3f}",f"{a.loc['IPF | Alveolar fibroblasts','median_reference_percentile']:.3f}"])
table(rows,[130,180,180]);p('补充图S2及表S2｜气泡颜色为平均log1p(CP10k)，面积为表达细胞比例。前6列为IPF相关基因，其余为基质和myCAF标志物；先按供者汇总，参考再按研究等权。表中使用肺泡成纤维细胞参照。较高六基因评分主要受COL1A1贡献，当前数据不支持完整典型IPF状态的认定。','small')
add=W/'HLCA_main_and_supplement.pdf';SimpleDocTemplate(str(add),pagesize=(595,842),leftMargin=48,rightMargin=48,topMargin=45,bottomMargin=45).build(story)
added=pymupdf.open(add);assert len(added)==6,len(added)
out=pymupdf.open();out.insert_pdf(baseline)
# Move the existing S1 table into the unified supplement; preserve original main body.
last=out[10];cut=last.search_for('补充表S1')[0].y0-3;last.add_redact_annot(pymupdf.Rect(35,cut,560,805),fill=(1,1,1));last.apply_redactions(images=0,graphics=0)
for i in range(10):assert out[i].get_text()==baseline[i].get_text()
assert re.sub(r'\s+','',''.join(out[10].get_text()))==re.sub(r'\s+','',baseline[10].get_text().split('补充表S1')[0])
out.insert_pdf(added)
for i,page in enumerate(out):page.insert_text((page.rect.width-65,page.rect.height-18),f'{i+1} / {len(out)}',fontsize=8,color=(.4,.45,.5))
out.set_toc([[1,'执行摘要与研究队列',1],[1,'UMAP与标志物图谱',2],[1,'患者构成与支持',4],[1,'患者级功能效应与敏感性',6],[1,'13个信号群与统计方法',8],[1,'文献与综合推荐',10],[1,'HLCA外部身份参照',12],[1,'正常肺背景比较',13],[1,'外部参照方法与综合解释',14],[1,'补充材料',15],[2,'S1 原主分析标志物结果',15],[2,'S2 IPF功能参照',16],[2,'S3 单基因及COL1A1敏感性',17]])
staged=W/'integrated_review.pdf';out.save(staged,garbage=4,deflate=True);out.close();doc=pymupdf.open(staged)
for i in list(range(10,17))+[0]:doc[i].get_pixmap(matrix=pymupdf.Matrix(1.15,1.15)).save(W/f'page_{i+1}.png')
alltext=''.join(p.get_text() for p in doc);assert len(doc)==17;assert '13' in alltext and all(n in alltext for n in nodes);assert 'IPF' not in ''.join(doc[i].get_text() for i in range(11,14));assert all(v in alltext for v in ['-0.163982','0.625000','0.984375','0.710374']);assert all(len(doc[i].get_images())>=len(baseline[i].get_images()) for i in range(10))
qa=dict(status='STAGED_PASS',pages=17,original_pages=11,main_HLCA_pages=[12,13,14],supplement_pages=[15,16,17],original_first10_pages_text_unchanged=True,original_main_conclusion_preserved=True,original_S1_relocated=True,main_HLCA_has_no_IPF=True,original_sha256=hashlib.sha256(backup.read_bytes()).hexdigest(),staged_sha256=hashlib.sha256(staged.read_bytes()).hexdigest(),target=str(target),backup=str(backup))
(W/'INTEGRATION_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(qa,ensure_ascii=False))
