import pathlib,json,csv,re,html,hashlib
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
P=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase109_three_tier_report_20261020')
pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/arphic/uming.ttc'))
pdfmetrics.registerFontFamily('CJK',normal='CJK',bold='CJK',italic='CJK',boldItalic='CJK')
body=ParagraphStyle('body',fontName='CJK',fontSize=10,leading=14.2,wordWrap='CJK',spaceAfter=6,textColor=colors.HexColor('#253543'))
small=ParagraphStyle('small',parent=body,fontSize=8.5,leading=11.6,spaceAfter=4)
title=ParagraphStyle('title',parent=body,fontSize=18,leading=25,spaceAfter=10,textColor=colors.black)
head=ParagraphStyle('head',parent=body,fontSize=12,leading=17,spaceBefore=8,spaceAfter=6,textColor=colors.HexColor('#155C69'))
sub=ParagraphStyle('sub',parent=body,fontSize=11,leading=15.5,spaceBefore=6,spaceAfter=5,textColor=colors.black)
def pa(t,s=body):return Paragraph(t,s)
def read(name):
 with (P/'tables'/name).open(encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
cand=read('tier_candidate_registry.tsv');six=read('final_six_complete_statistics.tsv')
bygene={r['gene']:r for r in six};names=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
def num(v,digits=3):
 if v in ('NA','','nan'):return 'NA'
 return f'{float(v):.{digits}g}'
def eff(v):return f'{float(v):.3f}'
def tb(rows,widths,fs=8.7):
 st=ParagraphStyle('cell',parent=small,fontSize=fs,leading=fs+3,spaceAfter=0)
 t=Table([[pa(str(v),st) for v in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAF2F5')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),.3,colors.HexColor('#D0DAE1')),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
 return t
story=[]
# Page 1
story += [pa('骨肉瘤CAF候选基因三级筛选与生物学解释',title),pa('独立分析报告  2026年10月20日',small),
pa('本报告从全基因差异结果逐级累加证据：原两队列各自名义显著且同向得到14个候选；加入NCI对骨比较及合并总体检验后保留8个；扣除正常组织背景后保留6个。最终6个均为肺部病灶候选细胞中较高表达，其中5个在两项正常背景对比中均通过全基因FDR。'),
pa('三级筛选结果',head)]
rows=[['档次','累加的筛选条件','数量'],
['第一档','原SCPCP000017和SCPCP000023各自P&lt;0.05，肺减骨方向一致','14'],
['第二档','第一档基础上，NCI肺对原骨P&lt;0.05且同向；合并总体比较也P&lt;0.05且同向','8'],
['第三档','第二档基础上，合并和NCI对骨两项扣除正常背景后，仍同向且各自P&lt;0.05','6']]
story += [tb(rows,[62,413,50],9),Spacer(1,7)]
story += [pa('最终保留PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1。PNKD在两项背景对比中均名义显著，但NCI背景对比FDR=0.0715；其余5个两项背景FDR均&lt;0.05。“五个FDR显著”不表示原两个队列的全基因FDR也显著。',small),
pa('14个初始候选的逐档去向',head)]
rows=[['基因','第一档','第二档','第三档']]
for r in cand:
 t2=r['tier2']=='True';t3=r['tier3']=='True'
 last='保留' if t3 else '参考缺失' if r['gene']=='AL096854.1' else '未达P阈值' if t2 else '未进入'
 rows.append([r['gene'],'保留','保留' if t2 else '未达该档',last])
story += [tb(rows,[165,95,115,150],8.6),
pa('第二档保留NEBL、PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1、AL096854.1。第三档中NEBL两项背景P分别为0.0518和0.0921；AL096854.1因正常参考缺失无法评估，不能作为阴性。',small),
pa('样本单位为患者。原队列17为骨12位、肺3位；原队列23为骨4位、肺3位；NCI为肺17位。合并比较为骨16位、肺23位，纳入6,633个细胞。NCI重复文库按已确认患者关系合并。',small)]
# Page 2
story += [PageBreak(),pa('最终六个候选的统计证据',title),pa('效应为肺减骨的log2(CP10k+0.1)均值差；背景对比再减去正常肺减骨差异。P为双侧原始P，q为全有效基因BH-FDR，未在6个内重算。',small),
pa('原两个队列各自名义显著',head)]
rows=[['基因','队列17效应','P17','q17','队列23效应','P23','q23']]
for g in names:
 r=bygene[g];rows.append([g,eff(r['source17_effect']),num(r['P_17']),num(r['FDR_17']),eff(r['source23_effect']),num(r['P_23']),num(r['FDR_23'])])
story += [tb(rows,[90,80,65,50,80,65,95],8.6),
pa('原队列23仅7位患者，按探索性证据解释；两原队列均没有全基因FDR显著基因，第一档采用的是原始P与方向一致性。',small),
pa('合并检验与NCI对骨检验',head)]
rows=[['基因','合并效应','合并P','合并q','NCI效应','NCI P','NCI q']]
for g in names:
 r=bygene[g];rows.append([g,eff(r['mean_lung_minus_bone']),num(r['P']),num(r['FDR']),eff(r['third_effect']),num(r['third_P']),num(r['third_FDR'])])
story += [tb(rows,[90,65,75,70,65,80,80],8.5),pa('上表为患者置换检验结果；它与此前Welch敏感性结果不同，报告不混用两种方法的P值。',small),
pa('扣除正常背景后的效应与显著性',head)]
rows=[['基因','合并效应 [95%CI]','P','q','NCI效应 [95%CI]','P','q']]
for g in names:
 r=bygene[g]
 def ci(side):return eff(r['effect_'+side])+'<br/>['+eff(r['CI_low_'+side])+', '+eff(r['CI_high_'+side])+']'
 rows.append([g,ci('pooled'),num(r['P_pooled']),num(r['FDR_pooled']),ci('NCI'),num(r['P_NCI']),num(r['FDR_NCI'])])
story += [tb(rows,[78,110,54,55,110,58,60],8.1),
pa('PNKD仍符合第三档P&lt;0.05条件，但不属于两项背景FDR均显著的5个。TAF10在未扣背景时仅为名义信号，扣背景后两项均通过全基因FDR；两种证据层级分别保留。',small),
pa('统计方法',head),
pa('每个患者部位汇总原始计数，以完整可用基因UMI为分母归一化。原两队列采用精确患者标签置换；合并及NCI对骨采用65,536次固定种子Monte Carlo置换并加观察值校正，种子为20261012。正常背景对比为“肿瘤肺减骨”减去“正常肺减骨”，用四组Welch-Satterthwaite检验计入参考方差。正常骨5位、肺3位，每位至少20个参考细胞。',small)]
# Pages 3-5
plots=[
('第一档火山图 原两队列',[('tier1_source17','图1  原队列17的全基因结果。14个候选由本图与图2的同向名义显著交集产生，不是本图的全部显著基因。'),('tier1_source23','图2  原队列23只有骨4位与肺3位，精确置换P值呈离散分布。两原队列均无FDR显著点，橙色仍表示有效的名义信号。')]),
('第二档火山图 合并与NCI比较',[('tier2_pooled','图3  原两个队列加NCI的合并总体比较，骨16位对肺23位。在第一档14个内叠加本图及图4的同向P阈值，第二档保留8个。'),('tier2_NCI','图4  NCI肺17位对原两队列合并骨16位。该项提供第三项比较支持，与总体分析共享患者，不是第三个独立骨肺队列。')]),
('第三档火山图 正常背景对比',[('tier3_pooled','图5  合并肺减骨效应扣除正常肺减骨差异，正常骨5位与正常肺3位。第二档8个中，AL096854.1不可评估，NEBL未通过P阈值。'),('tier3_NCI','图6  NCI肺对骨比较同样扣除正常背景。最终6个均同向名义显著，其中PNKD为橙色，其余5个为红色。')])]
legend='红色为全基因FDR&lt;0.05；橙色为仅原始P&lt;0.05；灰色为不显著。深色圈和基因名标出最终6个。横轴为肺减骨效应，纵轴为-log10原始P；虚线为P=0.05。每张图展示该项全部有效检验，图内总数不等同于逐档候选数。'
for heading,pp in plots:
 story += [PageBreak(),pa(heading,title),pa(legend,small)]
 for key,caption in pp:story += [Image(str(P/'figures'/(key+'.png')),width=525,height=267.3),pa(caption,small),Spacer(1,7)]
# Pages 6-7
bios=[
('PNKD 迁移与代谢相关候选',
'证据层级为其他癌种的体内外功能研究，包含亚型相关机制。[1,2]',
'PNKD存在不同蛋白亚型。肿瘤文献中的myofibrillogenesis regulator-1常写为MR-1，与抗原呈递基因MR1不同。非小细胞肺癌研究将MR-1与ITCH介导的NOTCH3胞内段泛素化和溶酶体调控及迁移侵袭联系；2026年肾透明细胞癌研究报告PNKD促进乳酸生成、细胞周期和肿瘤生长，并以PI3K抑制支持PI3K-AKT-mTOR相关解释。[1,2]',
'本项目中PNKD在三个层次均保留，可作为迁移与代谢适应的候选连接点。后续宜先明确表达细胞及亚型，再分别选择迁移或乳酸代谢作为主要终点，配合信号磷酸化、细胞存活及回补实验。共培养研究应区分成纤维细胞端与肿瘤细胞端干预。',
'定位与检测：随亚型而异的胞内蛋白；不应按MR1抗原呈递分子选抗体，也不作为直接活细胞表面标志。'),
('TOR1AIP1 核膜适应与受限迁移候选',
'证据层级为其他癌种的亚型机制研究，存在癌种间作用差异。[3,4]',
'TOR1AIP1编码核内膜相关蛋白LAP1，主要包括LAP1B和LAP1C。黑色素瘤研究显示LAP1C支持核膜重塑、受限迁移及三维胶原中的侵袭；前列腺癌研究则发现TOR1AIP1通过稳定p53抑制增殖与侵袭。两类结果共同说明亚型和细胞背景决定解释，不能概括为普遍促癌作用。[3,4]',
'本项目肺端升高与“不同基质孔隙及机械负荷下的核适应状态”这一研究假设相符。可采用亚型可分辨的RNA和蛋白检测，将核形态、核膜变化及受限迁移相结合。基因总表达量不能代替LAP1C比例或活性。',
'定位与检测：核内膜。核膜定位不等于细胞表面定位，不适合直接用作活细胞表面分选标志。'),
('GRID2 非经典受体表达候选',
'证据层级为基础受体注释和其他癌种的表达预后关联，尚非CAF功能机制。[5,6]',
'GRID2编码δ2型谷氨酸受体亚基GluD2，经典背景为神经系统突触组织。胃癌研究将GRID2纳入自噬相关基因预后模型，并报告整组织表达及结局关联；这不能单独证明GRID2在CAF内调控自噬。GRID2IP是另一基因，其实验机制不能替代GRID2证据。[5,6]',
'本项目的连续同向统计信号支持保留GRID2，但优先问题是表达的细胞归属和定位。先进行转录本特异性、原位RNA及蛋白定位验证，检查是否由少数细胞或其他谱系背景驱动。明确表达和膜定位后，再探索受体信号或自噬通量。',
'定位与检测：经典注释为膜受体；在本项目中仅属潜在膜相关候选，细胞表面可及性和抗体外部表位需另行验证。'),
('TAF10 转录装置与细胞状态候选',
'证据层级为其他癌种的调控和干预实验，细胞类型特异性有限。[7]',
'TAF10属于TATA结合蛋白相关因子，参与通用转录复合物。Z363研究在乳腺、肺及肝癌相关模型中分析TRIP12、TAF10和MYC，包含蛋白互作、泛素化、表达干预及动物实验；TRIP12调控TAF10降解，TAF10影响MYC转录。[7]',
'本项目可将TAF10作为肺端候选细胞转录状态变化的线索。验证时应区分广泛生物合成或增殖背景与成纤维功能，联合评估TAF10、MYC相关输出及细胞周期。干预与回补需配合一般转录和存活对照，避免把普遍转录下降解释为CAF特异功能抑制。',
'定位与检测：细胞核转录复合物；不适合直接活细胞表面分选。'),
('BEST1 离子稳态与细胞体积候选',
'证据层级为其他癌种的细胞功能实验，模型为上皮来源肿瘤细胞。[8]',
'BEST1编码bestrophin 1，与钙激活氯离子通道及离子稳态相关。T84结肠癌研究观察到快速生长细胞中BEST1升高；RNA干扰降低增殖，过表达促进生长，并涉及细胞内钙信号和体积调节。[8]',
'本项目中BEST1在背景扣除后仍保持肺高，可提出离子与细胞体积适应的研究问题。应先确认RNA及蛋白位于目标成纤维候选细胞，再检测细胞体积、离子或钙反应，并结合存活与增殖读数。若研究成骨或矿化，应设置独立基质矿化终点，不能以BEST1表达替代。',
'定位与检测：膜相关离子通道，定位依细胞背景而异；膜注释本身不足以确立可用的活细胞分选抗体。'),
('RGS2-AS1 非编码RNA状态候选',
'证据层级为规范注释及非肿瘤研究中的候选检测，已有阴性验证结果需保留。[9,10]',
'RGS2-AS1对应ENSG00000285280，旧名包括AL390957.1和LINC01032。它是非编码RNA，不能用RGS2蛋白的G蛋白信号机制解释其功能。牙周膜干细胞炎症与牙骨质分化研究将AL390957.1纳入验证，但qPCR未发现条件间显著差异；论文后续机制针对GACAT2与PKM1/2，不适用于RGS2-AS1。[9,10]',
'本项目肺端信号与该非肿瘤阴性结果来自不同细胞和条件，其合理定位为候选状态相关RNA。先统一稳定ID、转录本和探针覆盖，开展链特异性表达、原位及亚细胞定位；确认后再选择干预方式。邻近RGS2是否受影响可以作为研究问题，不能预设其为直接靶标。',
'定位与检测：长非编码RNA，亚细胞定位待验证；不能作为蛋白表面标志，可考虑固定细胞RNA检测。')]
(P/'resources/report_biology_text.json').write_text(json.dumps(bios,ensure_ascii=False,indent=2),encoding='utf-8')
for index in range(2):
 story += [PageBreak(),pa('六个基因的生物学解释 '+('一' if index==0 else '二'),title),
 pa('以下内容整理自2026年10月16日既有文献调查，并与本次统计结论衔接。文献用于提出研究问题和实验依据，不作为筛除候选的额外门槛。',small)]
 for heading,evidence,lit,interpret,loc in bios[index*3:index*3+3]:
  story += [pa(heading,head),pa(evidence,small),pa(lit),pa('本项目解释与验证方向  '+interpret),pa(loc,small),Spacer(1,5)]
# Page 8
story += [PageBreak(),pa('解释边界与参考文献',title),
pa('本报告支持六个基因作为骨肉瘤骨与肺部病灶成纤维表达候选中的探索性研究对象。原数据为扩展候选集合，NCI采用作者CAF标签，均不等同于经实验确证的纯CAF。此前身份筛选敏感性显示支持强度会减弱，因此本报告不宣称对细胞身份定义稳健。'),
pa('NCI为单核测序且仅有肺样本；第三项比较及合并分析共享患者，不能称为三个独立队列验证。正常骨参考来自骨/骨髓MSC相关细胞，正常肺参考为未受累肺组织成纤维细胞，均为外部参考。背景扣除不是严格配对组织或来源批次校正，也不把剩余差异直接定义为肿瘤特异机制。'),
pa('六个基因的统计支持与文献功能证据分别呈现。其他癌种机制不直接转移为骨肉瘤CAF机制；表达、通路或临床关联不作为因果证明。涉及活细胞分选时，应先验证细胞归属、亚型、定位与表位可及性。'),
pa('参考文献',head)]
refs=[
('Zhao W et al. Myofibrillogenesis Regulator-1 Regulates the Ubiquitin Lysosomal Pathway of Notch3 Intracellular Domain Through E3 Ubiquitin-Protein Ligase Itchy Homolog in the Metastasis of Non-Small Cell Lung Cancer. Advanced Science. 2024.','10.1002/advs.202306472','https://doi.org/10.1002/advs.202306472'),
('Zhang X et al. PNKD Mediates Lactate-Driven Tumor Progression in Clear Cell Renal Cell Carcinoma. Human Mutation. 2026.','10.1155/humu/8430177','https://doi.org/10.1155/humu/8430177'),
('LAP1 supports nuclear adaptability during constrained melanoma cell migration and invasion. Nature Cell Biology. 2023.','10.1038/s41556-022-01042-3','https://doi.org/10.1038/s41556-022-01042-3'),
('Li Z et al. TOR1 AIP1 interacts with p53 to enhance cell cycle dysregulation in prostate cancer progression. Molecular and Cellular Biochemistry. 2025.','10.1007/s11010-025-05276-1','https://doi.org/10.1007/s11010-025-05276-1'),
('NCBI Gene. GRID2 glutamate ionotropic receptor delta type subunit 2.','Gene 2895','https://www.ncbi.nlm.nih.gov/gene/2895'),
('Chen L et al. Establishment and verification of prognostic model for gastric cancer based on autophagy-related genes. American Journal of Cancer Research. 2021.','PMCID PMC8085875','https://pmc.ncbi.nlm.nih.gov/articles/PMC8085875/'),
('Xiong Y et al. Small molecule Z363 co-regulates TAF10 and MYC via the E3 ligase TRIP12 to suppress tumour growth. Clinical and Translational Medicine. 2023.','10.1002/ctm2.1153','https://doi.org/10.1002/ctm2.1153'),
('Spitzner M et al. Eag1 and Bestrophin 1 Are Up-regulated in Fast-growing Colonic Cancer Cells. Journal of Biological Chemistry. 2008.','10.1074/jbc.M703758200','https://doi.org/10.1074/jbc.M703758200'),
('NCBI Gene. RGS2-AS1.','Gene 102724954; ENSG00000285280','https://www.ncbi.nlm.nih.gov/gene/102724954'),
('Li X et al. LncRNA GACAT2 binds with protein PKM1/2 to regulate cell mitochondrial function and cementogenesis in an inflammatory environment. Bone Research. 2022.','10.1038/s41413-022-00197-x','https://doi.org/10.1038/s41413-022-00197-x')]
for i,(ref,ident,url) in enumerate(refs,1):story.append(pa(f'[{i}] '+html.escape(ref)+f' <link href="{url}" color="#155C69">'+html.escape(ident)+'</link>',small))
story += [Spacer(1,5),pa('数据核查  本次报告沿用冻结统计结果，不新增或改变检验。三档集合核对为14、8、6；正常背景两项FDR显著为5个。原始P与各自完整检验家族FDR分别保留，火山图展示全有效基因。文献解释保留既往调查的证据层级及关键阴性。',small)]
def footer(c,doc):
 c.saveState();c.setStrokeColor(colors.HexColor('#CFD8DE'));c.line(35,28,A4[0]-35,28)
 c.setFont('CJK',7);c.setFillColor(colors.HexColor('#647681'));c.drawString(35,16,'CAF三级筛选报告  2026-10-20  |  版本 phase109  |  P与全基因FDR分列')
 c.drawRightString(A4[0]-35,16,str(doc.page));c.restoreState()
out=P/'report/CAF_three_tier_six_genes_biology_20261020.pdf'
doc=SimpleDocTemplate(str(out),pagesize=A4,leftMargin=35,rightMargin=35,topMargin=30,bottomMargin=40,title='骨肉瘤CAF候选基因三级筛选与生物学解释',author='CAF项目分析',pageCompression=1)
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print('PDF_COMPLETE',hashlib.sha256(out.read_bytes()).hexdigest(),out.stat().st_size,flush=True)
