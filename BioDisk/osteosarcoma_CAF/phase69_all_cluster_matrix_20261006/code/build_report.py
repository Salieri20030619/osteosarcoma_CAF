
from pathlib import Path
import json,pandas as pd
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
W=Path(r'C:\Users\31937\CAF_all_clusters_20261006');F=W/'figures';D=W/'analysis/tables'
original=W/'report_style_reference.docx'
exec((W/'report_helpers.py').read_text(encoding='utf-8'))
q=pd.read_csv(D/'additional_nonprimary_source_signals.tsv',sep='\t')
allq=pd.read_csv(D/'nominal_upregulated_all_sources_with_sensitivity.tsv',sep='\t')
assert len(allq)==54 and allq.subgroup_id.nunique()==7 and not (allq.FDR_excluded<.05).any()
short=lambda x:x.replace('_original','-').replace('_','-')
for sec in doc.sections:
 for fp in sec.footer.paragraphs:
  fp.text='CAF 全簇胶原与基质证据核查  |  2026-10-06';fp.alignment=WD_ALIGN_PARAGRAPH.CENTER
  for rr in fp.runs:rr.font.size=Pt(8)
p('胶原与基质相关表达\n全簇补充核查报告','Title')
p('骨肉瘤 CAF 项目｜33 个亚群 × 所有已观察来源｜2026-10-06','Subtitle')
h('多个簇具有候选作用，F2-C3 不是唯一重点')
p('全来源核查显示，除 F2-C3 外，F2-C1、F4-C9 和 F1-C4 也有值得跟进的胶原或基质相关上调证据；F4-C5、F4-C8 等簇应保留为低重复数候选。当前数据支持按胶原类型、加工环节及谱系背景区分这些簇。')
p('此前报告以各簇预先确定的主要来源为筛选范围。本次扩展至全部 33 簇、116 个“簇–来源”组合，相关基因的名义上调记录由 45 条扩展为 54 条，涉及 7 个簇；新增 9 条来自 4 个此前未列入基因上调重点的簇。')
table(['簇 / 线索','本轮定位','主要依据'],[
['F2-C3','继续作为纤维型胶原主线','COL3A1、COL5A1/2、PCOLCE、LOXL2；I/III/V 型胶原程序名义上调。'],
['F2-C1','保留 IV 型胶原支路','COL4A1/2、PXDN、LOXL2；IV 型及交联相关程序名义上调。'],
['F4-C9','新增优先候选','COL6A3 在 S17 中 7/7 患者上调；Fibro 限定仍同向。'],
['F1-C4','新增加工相关候选','ADAMTS14 在 S17 名义上调；原始表达、min5、Fibro 均同向。'],
['F4-C5 / F4-C8','低重复数候选','胶原相关程序方向为正，主要来源有效 N 分别为 2 / 1。'],
['F1-C10 / F3-C0 / F2-C5','描述性补充候选','主要来源 I/III/V 型程序及 Fibro 限定方向均为正；有效 N=1。'],
['F2-C0','保留基因级观察','LOX、LUM、COL6A3 等线索；不支持整套合成/沉积过程同步增强。'],
['F4-C1 / F5-C12','有条件保留','存在基因级信号，但受谱系构成或来源差异影响。']
],[1.22,1.65,4.03],font=9.1)
p('表内 S17、S18、S23 分别为 SCPCP000017、SCPCP000018、SCPCP000023。F1/F5 的简写沿用原簇编号，例如 F1-C4 = F1_originalC4。','Caption')
h('建议调整后的结果主线')
p('以 F2-C3 的纤维型胶原表达为核心，同时呈现 F2-C1 的 IV 型胶原支路、F4-C9 的 VI 型胶原线索和 F1-C4 的加工相关线索。F4-C5、F4-C8 放入候选观察层。耐药与抑制机制继续暂缓。')
page()
h('新增 9 条名义上调记录，分布于 4 个簇')
p('筛选统一使用校正表达、剔除模型阳性双联体后的正向效应与原始 P < 0.05；FDR 保留原完整检验家族。以下均为名义显著，未达到 FDR < 0.05。')
rows=[]
for r in q.itertuples():
 rows.append([short(r.subgroup_id),r.dataset.replace('SCPCP0000','S'),r.gene,f'{int(r.positive_patients_excluded)}/{int(r.N_excluded)}',f'{r.mean_effect_excluded:+.3f}',f'{r.P_excluded:.5f}',f'{r.FDR_excluded:.3f}'])
table(['簇','来源','基因','正向/N','效应','原始 P','FDR'],rows,[1.4,.55,1.05,.9,1.,.95,1.05],font=9)
p('效应 = 患者平均 log2[(目标簇 CP10k + 0.1)/(同患者其他候选细胞 CP10k + 0.1)]。新增是本轮纳入展示的既有统计结果，不是重新检验得到的 P 值。','Caption')
fig('02_new_signals_sources','图 1｜新增簇–基因线索的全部可评估来源。* 表示原始 P < 0.05；灰色为无可评估结果。GSE299023 无同患者其他簇参照，未绘入此图。',6.7)
p('F4-C9 的 COL6A3 在 S17、S18 和 GSE318840 分别为 +2.803、+3.544、+1.133，正向患者为 7/7、3/3、3/3；后两者的 P 均为 0.25。F1-C4 的 ADAMTS14 在 S23 为反向。F5-C12 的 SPARC 在 S17/S18 名义上调，但其他来源方向不一致。')
p('ACTB 与 ACTG1 因广义基质组织基因集成员身份进入完整统计清单，不作为胶原沉积的特异支持基因。SDC2 主要对应细胞–基质相互作用；ADAMTS14、COL6A3、SPARC 的功能注释见参考来源 [1–4]。')
page()
h('全簇图谱显示不同类型的表达线索')
p('以下固定使用各簇原主要来源，展示全部 33 簇，防止只展示阳性候选。所有来源的完整结果另存统计表。')
fig('01_all33_programs','图 2｜33 簇的患者级程序效应。红色为正、蓝色为负；有效 N=1 的方向仅作描述。* 为原始 P < 0.05，不代表校正后显著，也不意味着所有程序成员同向变化。',6.6)
page()
h('表达分布补充呈现未达到筛选门槛的簇')
p('高表达、相对上调和已形成基质是不同证据层级。本图保留有表达但患者重复不足、或相对其他候选细胞未上调的簇。')
fig('03_all33_expression','图 3｜全部 33 簇的胶原及基质相关基因点图。各行使用其固定主要来源；点大小为校正计数 ≥ 1 的细胞比例，颜色为平均 log1p(CP10k)。每个来源的完整 31 基因数据保留在表达表；不把跨来源色深视作统计比较。',6.6)
page()
h('F4-C9 与 F1-C4 保留方向，F4-C1 对谱系限定敏感')
fig('04_new_signals_sensitivity','图 4｜新增 9 条记录的描述性敏感性分析。min5 要求同患者、同部位的目标和参照两侧均至少 5 个细胞；Fibro 限定两侧均为成纤维细胞。未新增 P 值。',6.7)
p('F4-C9 / COL6A3：S17 校正效应 +2.803，原始表达 +2.466，Fibro 限定 +1.656（N=7）；min5 为 +0.672，但只剩 N=1。F1-C4 / ADAMTS14：校正 +1.063（N=13），Fibro +0.921（N=8），min5 +1.424（N=4）。这两条均保留上调方向。')
p('F4-C1 / SDC2：Fibro 限定变为 −2.707（N=3），min5 为 −0.201（N=2）。S17 的该簇 53 个保留细胞中仅 3 个为 Fibro、49 个为 Osteo，因此不纳入稳健 CAF 基质候选。')
p('F5-C12 / SPARC：S17 的 Fibro 限定仍为正（+0.412，N=5）；S18 没有 Fibro 细胞，相关信号不能直接归为 CAF 来源。S18 的 76 个细胞中 70 个为未定谱系、6 个为 Osteo。')
h('低重复数候选保留，但不能合并成“已验证沉积簇”')
table(['簇','主要来源 I/III/V 型程序','核查结论'],[
['F4-C5','+1.157；N=2；P=0.50','Fibro 限定 +0.528（N=1）；GSE318840 的 COL3A1/LOXL2 反向。'],
['F4-C8','+1.195；N=1','Fibro 限定 +0.527（N=1）；保留基质组分候选。'],
['F1-C10 / F3-C0','+1.025 / +1.155；均 N=1','Fibro 限定同向，作为描述性胶原相关候选。'],
['F2-C5 / F4-C12','+1.281 / +0.916；均 N=1','F2-C5 的 Fibro 限定同向；F4-C12 原谱系约 99.5% 为 Osteo，单列。']
],[1.35,2.25,3.3],font=9)
page()
h('方法')
p('数据来自当前 Phase66 修复后的细胞证据表，共 8,174 个候选细胞，保留 F1–F5 的 33 个原簇标签；沿用双联体模型标记，剔除 26 个模型阳性细胞后分析 8,148 个细胞。读取既有 658 基因原始/校正矩阵、全库分母与检测覆盖掩码。未检测特征保留缺失，不填零。')
p('基因检验及程序统计继承 Phase66。每一来源内，以同一患者、同一取样部位的其他候选细胞为参照；基因表达先按目标/参照分别聚合为 CP10k，使用 0.1 伪计数计算 log2 表达比，再对部位、患者依次等权平均。程序沿用原评分定义、基因集及 full/drop_dominant_pair 版本。')
p('既有基因检验使用患者级符号翻转：N≤15 时为精确枚举，更大 N 使用既有 19,999 次随机翻转及观测值修正。无有效生物学重复或不可评估时不产生 P。筛选正向效应且原始 P<0.05，FDR 继承各来源、矩阵及相应同类假设家族，不按本次挑选的基因或簇重新缩小校正池。')
p('本次去除“只筛主要来源”的限制，覆盖 116 个簇–来源组合；复算 151,340 条既有基因基线记录，患者数、正向患者数一致，最大效应误差 1.78×10⁻¹⁵。敏感性增加每侧至少 5 细胞与 Fibro 限定的方向描述，不新增推断性检验。所有表保留来源，未跨数据集简单叠加患者数。')
h('讨论')
p('本次结果支持扩大候选簇范围。相对上调检验回答的是“该簇相对同患者其他候选细胞是否更高”，不回答“未上调的簇是否完全不表达或不参与”。N=1–2 的簇即使效应明确，也无法获得本研究筛选门槛下的统计支持，因而应保留方向和表达证据。')
p('F4-C9 的 COL6A3、F1-C4 的 ADAMTS14 是新增优先线索，但仍需考虑来源和谱系。F4-C9 在 S18 以 Osteo/未定谱系为主；各来源同向不等于已完成独立外部验证。F4-C5、F4-C8 等保留为低重复数候选；F3-C1 虽有历史基质相关注释，当前主要来源的胶原程序总体为负，LOX 为 −0.803（N=2），不据名称将其列为交联上调簇。F4-C10 的矿化相关身份也不能替代纤维型胶原装配证据。')
p('完整程序并非在所有候选中同步增强。F2-C3 的广义基质组织程序未上调；F2-C1 去除 COL4A1/2 后的 IV 型剩余程序转负；F2-C0 合成和基质组织程序为负。RNA 数据确立的是表达方向和候选功能关联，不直接测量胞外蛋白沉积、交联程度、酶活性或不同簇对物理基质的贡献比例。')
h('交付与后续重点')
p('本报告为前一版胶原与基质基因报告的全簇补充。完整结果包包含 33 簇逐簇核查表、全部来源统计、患者效应、谱系构成、表达表、敏感性分析和四幅图的 PNG/PDF。后续验证优先并列考虑 F2-C3、F2-C1、F4-C9、F1-C4；耐药与抑制机制继续暂缓。')
h('参考来源')
for t in [
'[1] NCBI Gene: ADAMTS14（胶原前肽加工）https://www.ncbi.nlm.nih.gov/gene/140766',
'[2] NCBI Gene: COL6A3（VI 型胶原及基质组织）https://www.ncbi.nlm.nih.gov/gene/1293',
'[3] NCBI Gene: SDC2（细胞–基质相互作用）https://www.ncbi.nlm.nih.gov/gene/6383',
'[4] NCBI Gene: SPARC（胶原结合）https://www.ncbi.nlm.nih.gov/gene/6678/'
]:p(t,'Caption')
out=Path(r'C:\Users\31937\Downloads\骨肉瘤CAF全簇胶原与基质核查补充报告_20261006.docx')
# Page-break-before on the next heading avoids standalone blank pages after tall figures.
for par in list(doc.paragraphs):
 if par._p.xpath('.//w:br[@w:type="page"]'):
  nxt=par._p.getnext()
  if nxt is not None and nxt.tag==qn('w:p'):
   from docx.text.paragraph import Paragraph
   Paragraph(nxt,doc._body).paragraph_format.page_break_before=True
  par._p.getparent().remove(par._p)
doc.save(out);print(out)
(W/'BUILD_RECORD.json').write_text(json.dumps({'report':str(out),'figures':len(doc.inline_shapes),'expected_pages':6,'all_source_nominal_rows':54},ensure_ascii=False,indent=2),encoding='utf-8')
