from pathlib import Path
import json, math, html
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak, KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PILImage
B=Path(r'C:\Users\31937\CAF_F1F5_round1_20261004'); D=Path(r'C:\Users\31937\Downloads\CAF_F1F5第一轮交付_20261004'); D.mkdir(exist_ok=True)
x=json.loads((B/'report_data.json').read_text(encoding='utf-8'));qa=json.loads((B/'QA_summary.json').read_text())
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\msyh.ttc'));pdfmetrics.registerFont(TTFont('CNBold',r'C:\Windows\Fonts\msyhbd.ttc'))
styles={
 'title':ParagraphStyle('title',fontName='CNBold',fontSize=23,leading=31,textColor=colors.HexColor('#17374E'),spaceAfter=14,wordWrap='CJK'),
 'h':ParagraphStyle('h',fontName='CNBold',fontSize=16,leading=23,textColor=colors.HexColor('#17374E'),spaceAfter=10,wordWrap='CJK'),
 'sub':ParagraphStyle('sub',fontName='CNBold',fontSize=10.5,leading=16,textColor=colors.HexColor('#176B77'),spaceBefore=9,spaceAfter=5,wordWrap='CJK'),
 'body':ParagraphStyle('body',fontName='CN',fontSize=9.3,leading=15.1,spaceAfter=7,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='CN',fontSize=7.4,leading=11.4,spaceAfter=5,wordWrap='CJK'),
 'cell':ParagraphStyle('cell',fontName='CN',fontSize=7.3,leading=10.8,wordWrap='CJK'),
 'head':ParagraphStyle('head',fontName='CNBold',fontSize=7.4,leading=11,textColor=colors.white,wordWrap='CJK')}
story=[];txt=[]
def p(s,sty='body'):
 story.append(Paragraph(s,styles[sty]));txt.append(s.replace('<b>','').replace('</b>',''))
def page(title):
 if story:story.append(PageBreak())
 p(title,'h')
def sub(s):p(s,'sub')
def fmt(v,n=3):
 if v is None or (isinstance(v,float) and not math.isfinite(v)):return 'NA'
 if isinstance(v,str):return v
 return f'{v:.{n}f}'
def prob(v):
 if v is None:return 'NA'
 if v==0:return '0*'
 return f'{v:.2g}' if v<.001 else f'{v:.4f}'
def table(headers,rows,widths,pad=5):
 data=[[Paragraph(html.escape(str(c)),styles['head']) for c in headers]]+[[Paragraph(html.escape(str(c)),styles['cell']) for c in r] for r in rows]
 t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#214F63')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#F0F5F7'),colors.white]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),pad),('BOTTOMPADDING',(0,0),(-1,-1),pad),('LINEBELOW',(0,-1),(-1,-1),.5,colors.HexColor('#B6C8D0'))]));story.append(t);story.append(Spacer(1,6));txt.append('\t'.join(headers));txt.extend('\t'.join(map(str,r)) for r in rows)
def fig(name,width=493,maxheight=550):
 path=B/'figures'/name;w,h=PILImage.open(path).size;scale=min(width/w,maxheight/h);story.append(Image(str(path),width=w*scale,height=h*scale));story.append(Spacer(1,5))
C=['F1','F2','F3','F4','F5'];totals=[695,225,251,5304,1699];NN=[13,12,3,15,13]
summary={
'F1':('F1 显示有限的基质相关候选特征，炎症型命名缺少支持','EXT1、PLXDC2 等出现名义上调；iCAF 程序未升高，基质和细胞周期程序整体较低。'),
'F2':('F2 的收缩与基质程序最清楚，伴随血管壁样特征','myCAF 与基质程序升高，层粘连蛋白互作、锚定纤维及平滑肌收缩 GSEA 富集。'),
'F3':('F3 呈翻译、增殖与能量程序的探索性富集','GSEA 显示 MYC、氧化磷酸化、翻译等方向；主比较仅 3 位患者，基因层面没有 P＜0.05 候选。'),
'F4':('F4 突出成骨相关表达和细胞周期，身份组成异质','RUNX2 等上调，成骨程序和细胞周期富集；复核映射中 Osteo 注释占较大比例。'),
'F5':('F5 以上调的线粒体转录本和偏低的基质程序为特征','名义上调基因主要为线粒体转录本；当前结果不支持将其解释为活跃的基质分泌型。')}
p('CAF 项目｜第一轮全转录组分析','title');p('原始拓展层 F1–F5 系统复核与 GSEA','h');p('交付版本：2026-10-04 · phase58 · 冻结原始 8,174 个细胞','sub')
p('本轮回答 F1–F5 各自具有什么特征。五类使用统一的全可检测转录组框架，涵盖基因、程序、内部状态、样本组成、表面候选及已有临床关联。')
table(['类别','细胞数','主比较患者 N','本轮结论'],[[k,n,nn,summary[k][1]] for k,n,nn in zip(C,totals,NN)],[36,47,65,345])
p('主比较患者数为 SCPCP000017 内同患者、同部位可配对的独立患者数。其他来源单独分析，不以细胞数替代样本量。','small')
sub('交付范围已完成')
p('完成 9 个来源的 45 个五类主比较登记、33 个冻结内部状态的分来源分析、校正前后比较，以及 22 个可评估主比较 × 4 个通路库的 88 份 GSEA。完整 GSEA 共 206,889 条通路结果；所有统计、leading-edge 基因及不可评估记录随附。')
p('独立 R 复核覆盖 25,972 个 P 值和 355 个 FDR 家族。P 值误差为 0，最大 FDR 误差为 2.22e-16。')
p('本轮未重新运行 CNV。沿用原始拓展层身份与已有校正数据；后续“蛋壳状 CAF”研究尚未启动。','small')
page('样本覆盖决定每类结论的证据强度')
p('原始对象保持为 8,174 个细胞、63,398 个 feature 的联合空间。分析按来源使用可检测 feature 并记录样本缺失测量，不把其他来源未测到的基因填成生物学零值。SCPCP000017 的来源空间为 60,319 个 feature。')
fig('01_patient_N.png',maxheight=290)
p('图 1｜每个来源中 Fk 对其余四类的独立可配对患者数。N=0 保留不可评估记录；N=1 只描述方向。','small')
table(['类别','有效基因检验','名义上调 / 下调','FDR＜0.05 上调 / 下调'],[[r['class_id'],r['valid_tests'],f"{r['nominal_up']} / {r['nominal_down']}",f"{r['FDR_up']} / {r['FDR_down']}"] for r in x['gene_summary'] if r['dataset']=='SCPCP000017'],[45,105,155,188])
p('表中基因候选以原始 P＜0.05 提名，完整 FDR 保留。主来源五类基因家族中没有 FDR＜0.05 的基因，本文均称“名义上调/下调”；这不否定其效应方向和校正前后一致性。程序和 GSEA 各有独立检验家族。','small')
program_names={'myCAF_customer':'收缩 / myCAF','iCAF_customer':'炎症 / iCAF','apCAF_customer':'抗原呈递 / apCAF','Osteogenic_customer':'成骨相关','Mineral_customer':'矿化相关','mCAF_matrix':'基质相关','stress_descriptive':'应激描述程序','cycle_descriptive':'细胞周期'}
extra={
'F1':'全局 695 个细胞中，后续盲审映射为 Fibro 的有 555 个（79.9%）。内部原始 C5 状态的 97.9% 细胞来自单个 donor；不能把该状态的聚类分离直接解释为稳定生物学亚型。PLXDC2、PTPRE 仅是转录候选，竞争表达未支持普遍特异性。',
'F2':'原有 6 个内部状态继续使用。主来源 C0/C1 的内部比较各有 6 位配对患者，其余状态为 2–5 位。血管壁竞争中 NOTCH3、MCAM 与 Mural 的差异不显著，PDGFRB 低于 Mural；因此收缩特征不能独自完成 CAF 与血管壁细胞的区分。',
'F3':'3 个内部状态的最大单 donor 贡献分别为 96.6%、95.2%、93.3%。主来源内部比较为 1、2、2 位患者；内部状态以组成和方向描述为主。主比较 N=3 时双侧精确符号翻转检验最小 P=0.25，没有名义差异基因不等于没有表达方向。',
'F4':'17 个内部状态存在明显患者集中，其中 3 个状态完全来自单个 donor，多个状态超过 99%。原始 5,304 个细胞中，后续映射为 Osteo 的有 3,297 个（62.2%），Fibro 为 909 个（17.1%）。这些构成结果与成骨方向一致，但不完成恶性身份或 CAF 身份判定。',
'F5':'4 个原始内部状态延续分析。C1、C12、C6 在主来源的内部比较 N 分别为 10、8、7；C23 无主来源配对。全局盲审映射为 Osteo 694 个、U 548 个、Fibro 65 个，显示该历史标签覆盖异质细胞。'}
for k,n,nn in zip(C,totals,NN):
 page(summary[k][0]);p(summary[k][1]);r=next(z for z in x['QC'] if z['dataset']=='SCPCP000017' and z['class_id']==k)
 p(f"全局 {n:,} 个细胞；主来源观察到 {r['patients']} 位患者、{r['cells']} 个细胞，可配对主比较 N={nn}。主来源中位 UMI={r['median_UMI']:,.0f}，中位检出基因={r['median_detected']:,.0f}。",'small')
 sub('患者级程序效应：原始与校正后')
 rows=[]
 for pn,label in program_names.items():
  a=next(z for z in x['program_comparison'] if z['dataset']=='SCPCP000017' and z['class_id']==k and z['program']=='customer|'+pn)
  rows.append([label,fmt(a['mean_delta_raw']),f"{a['positive_patients_raw']}/{a['N_raw']}",prob(a['P_raw']),prob(a['FDR_raw']),fmt(a['mean_delta_corrected']),prob(a['P_corrected'])])
 table(['程序','Δ 原始','正向/N','P 原始','FDR 原始','Δ 校正','P 校正'],rows,[112,53,48,65,70,70,75])
 p('Δ 为患者级 log2(CP10k+0.1) 差值的程序均值；正值表示该类高于其余四类。校正后的完整 FDR 见 program_correction_comparison.tsv.gz。','small')
 sub('代表性名义上调基因')
 genes=x['class_results'][k]['nominal_up_examples'][:5]
 if genes:table(['基因','Δ','正向/N','原始 P','FDR'],[[z['symbol'],fmt(z['mean_delta']),f"{z['positive_patients']}/{z['N']}",prob(z['P']),prob(z['FDR'])] for z in genes],[120,85,80,104,104])
 else:p('本比较没有原始 P＜0.05 的上调基因；未补造候选。全部基因效应和患者数保留在完整统计表。')
 a=next(z for z in x['gene_correction'] if z['dataset']=='SCPCP000017' and z['class_id']==k)
 if a['raw_up']:p(f"原始 {a['raw_up']} 个名义上调 feature 中，{a['raw_up_corrected_same_direction']} 个校正后同向，{a['raw_up_corrected_P_lt_05']} 个仍满足 P＜0.05。",'small')
 sub('内部组成与身份复核');p(extra[k])
page('患者级程序比较显示 F2 与 F4 的不同方向')
fig('02_program_raw.png',maxheight=280);fig('02_program_corrected.png',maxheight=280)
p('图 2｜SCPCP000017 原始与 DecontX 校正后的患者级程序效应。F2 的收缩和基质上升方向保持；F4 的成骨方向保持；F5 的基质程序下降方向保持。这里展示的是效应，统计结果见各类表及完整程序表。','small')
p('F2：myCAF 校正后 Δ=+2.269，P=0.0029，FDR=0.0092；基质 Δ=+1.615，P=0.0068，FDR=0.0171。F4：成骨程序校正后 Δ=+0.681，P=0.0221，FDR=0.0430。F5：基质校正后 Δ=−1.443，P=0.0046，FDR=0.0128。','small')
page('全排序 GSEA 补充程序结构，保留完整方向')
fig('03_GSEA_all_hallmark.png',width=465,maxheight=624)
p('图 3｜主来源全部 50 个 Hallmark 的 NES。正值表示位于该类相对其他四类更高的排序端。颜色不等于患者级检验显著性；完整模型 P、GSEA FDR q-value、BH FDR 和 leading edge 均在结果包。','small')
page('GSEA 显示五类的不同通路结构')
terms={'F1':['HALLMARK_MYC_TARGETS_V1','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_PROTEIN_SECRETION'],'F2':['HALLMARK_ANGIOGENESIS','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_G2M_CHECKPOINT'],'F3':['HALLMARK_MYC_TARGETS_V1','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_E2F_TARGETS'],'F4':['HALLMARK_MITOTIC_SPINDLE','HALLMARK_G2M_CHECKPOINT','HALLMARK_ANGIOGENESIS'],'F5':['HALLMARK_PROTEIN_SECRETION','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_GLYCOLYSIS']}
rows=[]
for k in C:
 for term in terms[k]:
  a=next(z for z in x['class_results'][k]['hallmark'] if z['Term']==term);rows.append([k,term.replace('HALLMARK_',''),fmt(a['NES']),prob(a['GSEA_model_P']),prob(a['model_BH_FDR'])])
table(['类','代表性 Hallmark','NES','模型 P','模型 BH FDR'],rows,[28,241,60,77,87])
p('这是按各类结果主线选择的说明性展示，非候选筛选阈值。全部 Hallmark 在上一页，全部四库结果在结果包。0* 表示有限置换中未观察到相应尾部事件；原值保留为 0，不代表真实概率为零。','small')
sub('Reactome 的代表方向')
p('F2：LAMININ_INTERACTIONS（NES=2.475）、ANCHORING_FIBRIL_FORMATION（2.343）与 SMOOTH_MUSCLE_CONTRACTION（2.337）均呈正向模型富集，模型 P=0*。F3：翻译延伸、翻译起始及蛋白靶向膜等通路正向富集。F4：核膜分解、纺锤体检查点等周期程序正向富集。')
p('F1/F5 也有补体、免疫或特定代谢名称的局部正向 Reactome 条目。解释这些条目需同时查看实际覆盖和 leading-edge；单个小基因集的富集不直接定义全类功能。F2 的 TGF-β Hallmark 为负向，不能因基质上升而把所有相关通路一律写成激活。')
sub('模型证据与患者证据分别记录')
p('GSEA 使用全量可评估 canonical gene 排名，未先按 P 值筛基因。88 份结果包括 44 份本轮计算和 44 份排序及效应逐项一致后复用的历史计算。F3 的模型富集不扩大独立患者 N，也不替代其患者级统计。')
page('内部状态中存在明显的患者集中')
fig('05_internal_donor_concentration.png',maxheight=542)
p('图 4｜冻结内部状态的最大 donor 贡献比例，汇总全部来源。保留 F2 的 6 个、F3 的 3 个、F4 的 17 个既有状态；F1 的 3 个和 F5 的 4 个状态来自原始拓展层簇。未为本轮叙事调节聚类参数。','small')
p('F1/F2/F3/F4/F5 各状态的全基因效应、1,867 个程序以及患者组成均已交付。高患者集中度的状态先作为来源/患者相关表达结构解释；少于 2 个可配对患者的比较不输出推断性 P 值。')
page('内部比较提供可追溯的状态特征')
rows=[]
for k in C:
 regs=[r for r in x['internal_registry'] if r['dataset']=='SCPCP000017' and r['parent']==k]
 for rr in regs:
  z=[q for q in x['internal_main_source'] if q['contrast']==rr['contrast'] and q['P'] is not None]
  z=sorted(z,key=lambda q:(q['P'],-abs(q['mean_delta'])))
  a=z[0] if z else None
  rows.append([rr['target'].replace('expanded_',''),rr['N'],program_names[a['program'].split('|')[1]] if a else '无有效检验',fmt(a['mean_delta']) if a else 'NA',prob(a['P']) if a else 'NA'])
table(['状态','配对 N','客户程序中最小 P 项（描述）','Δ','P'],rows,[150,43,155,65,80],pad=2)
p('本表每个状态仅显示预设 8 项客户程序中最小 P 项，方便定位，不能作为额外显著性筛选；完整同类假设家族的 FDR 和全部程序结果在 internal_customer_programs.tsv.gz / all_program_statistics.tsv.gz。小样本结果按探索性呈现。','small')
page('表面候选需结合跨谱系竞争表达解释')
p('采用“该类相对其余四类名义上调（P＜0.05）＋已有质膜跨膜注释”提名，得到 111 条来源-类别候选记录、108 个唯一基因。随后回到 107 个样本的 876,648 个细胞，与 8 类背景进行患者级竞争比较，共 888 条统计。')
rows=[]
for k in C:
 a=[z for z in x['surface'] if z['class_id']==k];names=sorted(set(z['symbol'] for z in a));rows.append([k,len(names),'、'.join(names) if k!='F4' else '、'.join(names[:14])+'……（完整 94 个见表）'])
table(['类别','唯一基因数','转录候选'],rows,[38,62,393])
sub('F2 中 ITGB1、MMP14 的相对表达较有区分度')
p('在当前 8 类背景比较中，ITGB1 和 MMP14 均为正向且 P＜0.05。NOTCH3 与 Mural：Δ=−0.044，P=0.9468；MCAM 与 Mural：Δ=−0.560，P=0.3892；PDGFRB 与 Mural：Δ=−1.756，P=0.0007，FDR=0.0031。收缩类候选与血管壁标记的重叠必须保留。')
sub('F1 的两个候选均未形成广泛的特异性证据')
p('PLXDC2 相对其他 Fibro-like：Δ=−0.698，P=0.4487；相对 Myeloid：Δ=−2.250，P=0.0015。PTPRE 相对其他 Fibro-like：Δ=−0.118，P=0.8423。不能将这两个名义上调基因直接交付为特异分选标志。')
sub('F3/F5 没有通过当前提名规则的表面候选')
p('保留真实的空结果。F4 的大量候选需结合其成骨样及其他身份组成解释。此处输出的是 RNA 和注释层面的候选，不包含抗体性能验证、蛋白定位验证或现成分选方案。')
page('已有临床关联保留严格部位与治疗口径')
p('本轮复用已经完成的临床分析及患者清单，并保留来源文件和校验值。五类全转录组特征分析不改变临床分组、不重定义治疗暴露，也不把旧的共享基因临床结果标成新全转录组临床验证。')
sub('SCPCP000023 的记录治疗分组显示捕获比例下降')
table(['对象','No / Yes 独立患者','No 比例','Yes 比例','差值 Yes−No','P / FDR'],[['全部拓展候选','11 / 15','1.4209%','0.0173%','−1.4036 百分点','0.0070 / 0.0358'],['F4','11 / 15','1.3164%','0.0148%','−1.3016 百分点','0.0119 / 0.0358'],['F5','11 / 15','0.0911%','0.0007%','−0.0904 百分点','0.0215 / 0.0429']],[73,80,75,75,95,95])
p('来源：phase33 的 clinical_proportion_statistics.tsv，分母为对应样本完整 called-cell 库。此结果描述记录治疗分组与捕获细胞比例的关联。','small')
sub('严格骨/肺部位复核未得到五类比例的名义差异')
p('严格部位口径下，拓展层主类别没有 P＜0.05 的比例结果。历史 F2 内部 subgroup1 的宽口径骨/肺比较曾为 N=9/3、P=0.0364；严格修订后为 N=9/2、P=0.0727，不再满足名义筛选标准。报告采用严格修订口径。')
sub('调整模型与历史表达分析分别保留')
p('多因素比例模型使用既有 HC3 修订的 t 分布 P_t/FDR_t；旧 normal-P 列仅为历史留存。临床表达与程序的历史结果、严格部位表、临床 roster 及比较可用性清单随包保留；旧分析的共享基因范围不等同于本轮全 feature 范围。')
p('临床比例、表达、程序和调整模型不混入同一个 FDR 校正池。各文件按自身预先定义的检验家族呈现完整效应、原始 P 与 FDR。')
page('质量指标与盲审映射揭示历史标签的组成')
fig('04_QC.png',maxheight=200)
rows=[]
for k in C:
 z=[r for r in x['lineage'] if r['class_id']==k];z=sorted(z,key=lambda r:-r['cells']);total=sum(r['cells'] for r in z)
 rows.append([k,'；'.join(f"{r['review_lineage_code']} {r['cells']}（{r['cells']/total:.1%}）" for r in z[:3])])
table(['类别','后续盲审映射中的前三类（全部来源）'],rows,[45,448])
p('图 5｜主来源 QC 为描述性结果，未据此新增细胞阈值或改写冻结类别。盲审标签仅叠加到原始细胞，用来核查身份构成；并未用新的大范围候选池替换本轮 8,174 个细胞。','small')
sub('校正前后方向具有一致性')
table(['类别','原始名义基因数','校正后同向','校正后仍 P＜0.05'],[[k,a['raw_nominal'],a['raw_nominal_same_corrected_direction'],a['raw_nominal_corrected_P_lt_05']] for k in C for a in x['gene_correction'] if a['dataset']=='SCPCP000017' and a['class_id']==k],[44,142,135,172])
p('校正一致性表明上述方向不是仅由本次所用环境 RNA 校正改变而产生。其含义限于这一敏感性比较，不能代替谱系、恶性身份或功能实验。','small')
page('方法：冻结对象、患者级比较与通路分析')
sub('对象与效应量')
p('从 phase30 原始拓展层冻结细胞和 F1–F5 标签，使用 phase41 的完整基因矩阵及样本 feature 存在矩阵。每个来源内，患者×部位分别汇总目标与参照细胞计数，用各自全库 UMI 转成 CP10k，再计算 log2(CP10k+0.1) 差值。同患者多部位效应等权平均，患者之间等权汇总。某比较涉及样本未测的 feature 标记 NA。程序为固定基因集内有效基因效应均值。')
sub('检验与校正')
p('患者级双侧符号翻转检验：有效 N≤15 时枚举所有符号组合；N＞15 时使用固定种子 20261001 的 19,999 次随机翻转并计入观测。N＜2 或全零效应不作检验。五类主比较在来源内组成全基因家族；内部比较按来源与父类组成家族；程序另设家族。BH 校正对未检验项按 P=1 占位计入冻结家族，再将其结果恢复为 NA。原始目标 P＜0.05 用于名义提名，FDR＜0.05 标记校正后显著。')
sub('全转录组 GSEA 与 ORA')
p('canonical symbol 映射及代表 feature 在评估可检验性之前固定。GSEA 以全部有效、可评估 canonical 基因的平均患者效应降序排序，包含正负方向，不按差异 P 值预筛。使用冻结 Hallmark、Reactome、GO BP、Matrisome 四库；1,000 次置换，种子 20261001，min_size=1、max_size=背景基因数−1。保留 NES、原始模型 P、软件 GSEA FDR q-value、另算同来源/比较/库内 BH FDR，以及 leading-edge。')
p('ORA 独立以原始基因 P＜0.05 的上/下调集合进行，背景使用对应来源可评估 canonical 基因；同一比较、库和方向构成各自的校正家族。完整通路覆盖、缺失项及不可评估记录保留。')
sub('内部状态、校正与竞争表达')
p('沿用冻结状态标签，状态与同父类其余细胞比较。DecontX 敏感性分析复用既有 corrected_counts，不重新估计模型，核验细胞/基因顺序后构建完整 feature 校正对象。表面提名基于 nominal up 和冻结质膜跨膜注释；回到全部 107 样本，按患者×部位与 8 类背景比较。各来源、目标类的候选×背景构成竞争表达家族。')
sub('计算复用与独立复核')
p('相同旧比较在 feature 顺序、患者效应及对比一致时复用原始统计，再按本轮五类家族重算 FDR；GSEA 仅在全排序基因及效应一致时复用。R 独立实现复算 25,972 个选定 P 值（含全部表面竞争 P）及全部新增相关 FDR 家族；检查所有 GSEA 输出基因归属、P 范围与有限 NES。')
page('讨论与交付索引')
sub('本轮建立的主线')
p('五类并非同一个基质程序的强弱排列。F2 的收缩/基质方向最清楚，F4 的成骨/周期方向突出，F5 呈线粒体转录本占优及基质偏低。F1 的炎症型命名缺少支持；F3 提供有价值的程序排序方向，但证据层级受独立患者数限制。')
sub('统一解释边界')
p('所有 Δ 和 NES 均为相对其他类的表达或排序效应，不能直接外推绝对分泌量。多数外部来源样本较少，内部状态也有患者集中；跨来源不合并细胞以制造样本量。F3 的 N=3 和部分内部状态的低 N 限制患者级推断；GSEA 模型 P 不替代生物学重复。小覆盖通路、有限置换零 P 及 overlapping gene sets 需结合完整覆盖表解读。')
p('原始历史 F 标签与后续盲审谱系映射不同，尤其 F4/F5 组成异质。RNA 候选不等于蛋白特异性，记录治疗关联不等于治疗因果效应。本轮没有空间环绕、胶原包被、免疫排斥或保护肿瘤的功能验证，因此没有把任何类命名为“蛋壳状 CAF”。')
sub('结果包主要入口')
table(['目录 / 文件','内容'],[['tables/genes/','各来源主/内部全基因统计，校正统计与逐基因校正比较'],['tables/all_program_statistics.tsv.gz','1,867 个固定程序的完整患者级结果'],['tables/all_GSEA_results.tsv.gz','四库全部 GSEA 结果及 leading-edge；单文件位于 gsea/'],['tables/ora/ 与 pathway_coverage.tsv.gz','完整 ORA 和基因集覆盖'],['tables/surface_*','表面提名、覆盖、患者表达、8 类竞争及汇总'],['tables/reused/','临床历史/严格修订、患者清单及复用结果'],['independent_QA/','独立 R 复算、全部 FDR 核查和 GSEA 输入输出核查'],['code/、resources/、delivery/','脚本、固定定义、输入来源、校验值和验收记录']],[190,303])
p('大型 H5AD 与患者效应 NPZ 保留在 sgs10 的 phase58_F1F5_round1_20261004/objects/，通过交付清单提供路径与 SHA256。结果包包含可直接查阅的完整统计和全部分析代码。第一轮交付后，再另行启动蛋壳状 CAF 的空间与功能定向研究。','small')

def decorate(c,doc):
 c.setStrokeColor(colors.HexColor('#D5E1E6'));c.line(50,38,545,38);c.setFont('CN',7);c.setFillColor(colors.HexColor('#607786'));c.drawString(50,25,'CAF | F1–F5 第一轮全转录组复核 | 2026-10-04');c.drawRightString(545,25,str(doc.page))
path=D/'CAF_F1F5第一轮全转录组与GSEA报告.pdf'
doc=SimpleDocTemplate(str(path),pagesize=(595.28,841.89),rightMargin=50,leftMargin=50,topMargin=45,bottomMargin=50,title='CAF F1–F5 第一轮全转录组与 GSEA 报告',author='CAF project analysis')
doc.build(story,onFirstPage=decorate,onLaterPages=decorate)
(D/'CAF_F1F5报告纯文本.txt').write_text('\n\n'.join(txt),encoding='utf-8');print(path)