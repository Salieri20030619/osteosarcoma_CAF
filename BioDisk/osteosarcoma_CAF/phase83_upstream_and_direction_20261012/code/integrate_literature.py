from pathlib import Path
from datetime import datetime,timezone
import json
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012'
rows=[]
def add(gene,aliases,evidence,context,finding,boundary,judgment,refs,queries):
    rows.append(dict(gene=gene,aliases=aliases,evidence_type=evidence,literature_context=context,literature_finding=finding,transfer_boundary=boundary,combined_interpretation=judgment,references=refs,search_queries=queries))
add('NEBL','LASP2; LIM-nebulette; LNEBL','其他肿瘤功能实验及临床关联；亚型明确','非小细胞肺癌、结直肠癌、肾透明细胞癌；黏附/骨架',
    'NSCLC研究报道LASP2促进迁移侵袭并关联FAK信号；结直肠癌研究观察NEBL与淋巴结转移关联。另有肾透明细胞癌研究显示NEBL过表达抑制迁移侵袭，作用随模型不同。',
    '本次基因级计数不能分辨NEBL不同亚型；LASP2是该基因短亚型，不是LASP1。未取得直接证明骨肉瘤CAF中同机制的研究。',
    '4来源同向、1零差，加上黏附/迁移文献，使其成为优先研究线索；仍保留Fibro限定和配对方向限制，不能称已验证CAF特异marker。',
    [('Lasp2 enhances tumor invasion via facilitating phosphorylation of FAK and predicts poor overall survival of non-small cell lung cancer patients','PMID28667800','https://pubmed.ncbi.nlm.nih.gov/28667800/'),('Nebulette Expression Is Associated with Lymph Node Metastasis in Patients with Colorectal Cancer','PMID30186581; DOI10.15171/mejdd.2018.107','https://pubmed.ncbi.nlm.nih.gov/30186581/'),('Downregulation of NEBL promotes migration and invasion of clear cell renal cell carcinoma by inducing epithelial-mesenchymal transition','DOI10.1016/j.prp.2023.155068','https://www.sciencedirect.com/science/article/abs/pii/S0344033823007690'),('NCBI NEBL gene annotation','Gene10529','https://www.ncbi.nlm.nih.gov/gtr/genes/10529/')],
    ['NEBL osteosarcoma fibroblast metastasis','"NEBL" cancer metastasis','"LASP2" "osteosarcoma"','"LASP2" "lung" migration invasion','"LASP2" "NEBL" gene isoform'])
add('AL096854.1','ENSG00000286540','骨肉瘤计算预后关联；非机制验证','TARGET-OS青少年骨肉瘤bulk RNA',
    '在细胞衰老相关lncRNA预后研究的单因素候选表中出现，HR约0.606、P约0.0478。该文最终六lncRNA模型不包含AL096854.1。',
    '不是已验证独立预后因子，也不是CAF或肺转移机制证据；不能把候选表条目写成最终模型验证。',
    '3来源同向、其余缺测；具已有OS关联背景，可保留为lncRNA探索候选，但不以此推断促进肺转移。',
    [('Cellular Senescence-Related Long Non-coding RNA Signatures Predict Prognosis in Juvenile Osteosarcoma','PMCID PMC11666862; Table1 candidate vs final six-lncRNA model','https://pmc.ncbi.nlm.nih.gov/articles/PMC11666862/')],
    ['"AL096854.1"','"AL096854.1" cancer osteosarcoma','"AL096854.1" "signature" osteosarcoma'])
add('GRID2','glutamate ionotropic receptor delta type subunit2','其他肿瘤计算预后及组织表达证据','胃癌',
    '胃癌自噬相关基因研究把GRID2纳入预后模型，并展示肿瘤组织表达。主要为关联性研究，不是GRID2干预证明自噬机制。',
    'GRID2IP是另一基因，不能引用其结直肠癌实验作为GRID2证据；本轮未取得直接骨肉瘤CAF实验。',
    '3来源同向、2零差，适合保留为部位表达候选；文献仅增加其他癌种关联背景。',
    [('Establishment and verification of prognostic model for gastric cancer based on autophagy-related genes','PMCID PMC8085875','https://pmc.ncbi.nlm.nih.gov/articles/PMC8085875/')],
    ['"GRID2" cancer','GRID2 cancer fibroblast osteosarcoma','"GRID2" "gastric" prognostic'])
add('RGS2-AS1','AL390957.1; LINC01032; ENSG00000285280','规范注释；非肿瘤实验中的检测线索','牙周膜干细胞炎症/牙骨质分化研究；非癌症',
    'NCBI确认RGS2-AS1/LINC01032与ENSG00000285280。旧名称AL390957.1曾进入牙周膜干细胞研究的候选检测，但qPCR未显示炎症与非炎症条件间显著改变；该文主机制是GACAT2。',
    '不能移用GACAT2机制；本轮未检出充分的RGS2-AS1特异骨肉瘤/CAF机制研究。旧文探针与当前亚型对应仍需核对。',
    '3来源同向、1零差、1缺测，且同来源支持层有名义信号；可作lncRNA候选，机制和配对一致性尚不足。',
    [('RGS2-AS1 NCBI Gene','Gene102724954; HGNC49018','https://www.ncbi.nlm.nih.gov/gene/102724954'),('LncRNA GACAT2 binds with protein PKM1/2 to regulate cell mitochondrial function and cementogenesis in an inflammatory environment','PMCID PMC8927299; AL390957.1 negative qPCR result','https://pmc.ncbi.nlm.nih.gov/articles/PMC8927299/')],
    ['"RGS2-AS1"','"RGS2-AS1" cancer','"AL390957.1" cancer','"LINC01032" cancer'])
add('PEF1-AS1','PEF1 and COL16A1 antisense RNA1; ENSG00000235790','本轮仅定位到转录本/基因注释','非编码RNA注释',
    'NCBI存在PEF1-AS1的长非编码RNA参考序列。定向检索未定位到可用于本项目的基因特异肿瘤功能原始研究。',
    '未检出不等于证明文献不存在或已具新颖性；PEF1、COL16A1的研究不能移作PEF1-AS1机制。',
    '3来源同向，但来源18效应仅约+0.036；作为低效应、机制缺口较大的探索线索保留。',
    [('PEF1-AS1 reference transcript variant2','RefSeq NR_184312.1','https://www.ncbi.nlm.nih.gov/nuccore/NR_184312.1')],
    ['"PEF1-AS1" cancer','"PEF1-AS1"'])
add('EBF1','EBF transcription factor1','骨肉瘤细胞功能实验；正常骨髓基质遗传实验','骨肉瘤细胞、鼠骨髓血管周基质',
    '骨肉瘤细胞研究报道EBF1激活FGD5-AS1/miR-124-3p/G3BP2轴。另有Prx1系Ebf1缺失研究显示骨髓基质/骨稳态改变，支持其骨间充质生物学作用。',
    '癌细胞促肿瘤作用与正常骨基质功能不是同一细胞背景；不能据文献消除本数据的反向结果。',
    '可作为骨谱系/机制对照；虽然4来源骨侧较高，来源18及来源17配对方向相反，不作为稳定部位marker的单一首选。',
    [('EBF1-mediated up-regulation of lncRNA FGD5-AS1 facilitates osteosarcoma progression by regulating miR-124-3p/G3BP2 axis as a ceRNA','PMID35761386; DOI10.1186/s13018-022-03181-7','https://pubmed.ncbi.nlm.nih.gov/35761386/'),('Deletion of the transcription factor EBF1 in perivascular stroma disrupts skeletal homeostasis and precipitates premature aging of the marrow microenvironment','PMID39002837; DOI10.1016/j.bone.2024.117198','https://pmc.ncbi.nlm.nih.gov/articles/PMC11410106/')],
    ['"EBF1" osteosarcoma fibroblast','"FGD5-AS1" "G3BP2"','"s13018-022-03181-7" retraction','"EBF1" "Deletion" "perivascular" stroma'])
add('NDUFS7','complex I subunit PSST','骨肉瘤药物研究补充表达线索；其他癌细胞药物响应','骨肉瘤耐阿霉素模型及HEp-2细胞',
    '骨肉瘤线粒体靶向阿霉素研究的补充表达表列出NDUFS7；其他癌细胞药物处理研究也检测其表达变化。',
    '这些结果不等于NDUFS7单基因干预已证明骨肉瘤CAF或肺转移机制，复合体I其他亚基实验也不能直接替代。',
    '有代谢/耐药参考价值；4来源同向、来源18反向，保留为机制探索而非普遍稳定marker。',
    [('Mitochondria-Targeted Doxorubicin: A New Therapeutic Strategy against Doxorubicin-Resistant Osteosarcoma','PMID27466354; DOI10.1158/1535-7163.MCT-16-0048; NDUFS7 in supplement','https://iris.unito.it/retrieve/handle/2318/1610727/259822/Buondonno%2C%20AOP%20Mol%20Cancer%20Ther%20MS%20and%20Supplementary%2C%202016.pdf'),('Modulation of Mitochondrial and Epigenetic Targets by Polyphenols-rich Extract from Araucaria angustifolia in Larynx Carcinoma','PMID30113000','https://pubmed.ncbi.nlm.nih.gov/30113000/')],
    ['"NDUFS7" cancer osteosarcoma','"NDUFS7" cancer proliferation','"NDUFS7" "Targeting Mitochondria" osteosarcoma'])
add('PNKD','paroxysmal nonkinesigenic dyskinesia gene','其他肿瘤基因干预及通路实验','肾透明细胞癌',
    '2026年研究报道PNKD促进乳酸相关代谢、细胞周期和肿瘤生长，涉及PI3K-AKT-mTOR，并用PI3K抑制作功能逆转。',
    '是肾癌模型，不能直接转述为骨肉瘤CAF中已验证的机制。',
    '代谢机制可参考；4来源同向但来源18反向，跨来源marker优先级受限。',
    [('PNKD Mediates Lactate-Driven Tumor Progression in Clear Cell Renal Cell Carcinoma','DOI10.1155/humu/8430177; 2026-03-27','https://onlinelibrary.wiley.com/doi/10.1155/humu/8430177')],
    ['"PNKD" cancer metastasis','"PNKD Mediates Lactate" retraction'])
add('BEST1','bestrophin1','其他肿瘤细胞干预实验','T84结肠癌细胞',
    '原始研究显示BEST1在快速生长T84细胞中增高；RNA干扰降低增殖，过表达可促进生长，涉及细胞内钙信号和体积调节。',
    '不能把结肠癌细胞结果解释为CAF特异性，也不证明钙沉积或矿化功能。',
    '具功能文献，但本数据3来源同向、2反向，且来源23 Fibro限定无差异，作为通用marker的依据有限。',
    [('Eag1 and Bestrophin1 are up-regulated in fast-growing colonic cancer cells','PMID18222922; DOI10.1074/jbc.M703758200','https://pubmed.ncbi.nlm.nih.gov/18222922/')],
    ['"BEST1" cancer fibroblasts','Bestrophin1 cancer proliferation'])
add('CALHM5','FAM26E; C6orf188','蛋白结构原始研究；未取得直接肿瘤机制证据','重组人CALHM5冷冻电镜结构',
    '有直接的人CALHM5冷冻电镜研究，解析其大孔道装配和膜脂相关结构。',
    '结构研究不是骨肉瘤CAF功能或矿化实验；CALHM1/2/6等家族成员的肿瘤免疫作用不能移用于CALHM5。',
    '3来源同向、1反向、1零；两名配对患者同向，但来源23 Fibro限定零效应。可保留结构明确、肿瘤功能待定的探索候选。',
    [('Cryo-EM structures of human calcium homeostasis modulator5','DOI10.1038/s41421-020-00228-z; PMCID PMC7652935','https://pmc.ncbi.nlm.nih.gov/articles/PMC7652935/'),('CALHM5 NCBI Gene','Gene254228','https://www.ncbi.nlm.nih.gov/gene/254228')],
    ['"CALHM5" cancer','"CALHM5" "FAM26E" cancer study','"CALHM5" "Liu" "2020" structure'])
add('ZNF526','zinc finger protein526','乳腺癌遗传关联及细胞存活实验','汉族女性乳腺癌、乳腺癌细胞',
    '外显子组关联研究定位ZNF526编码变异，并有表达与细胞存活相关实验支持。',
    '属于乳腺癌易感性/细胞背景，非骨肉瘤CAF或肺转移直接证据。',
    '文献增加生物学背景，但本数据3来源同向、2反向，不能用文献覆盖方向不一致。',
    [('A Large-Scale, Exome-Wide Association Study of Han Chinese Women Identifies Three Novel Loci Predisposing to Breast Cancer','PMID29572226; Cancer Research2018 78:3087-3097','https://aacrjournals.org/cancerres/article-split/78/11/3087/625133/A-Large-Scale-Exome-Wide-Association-Study-of-Han')],
    ['"ZNF526" cancer','"ZNF526" "rs3810151"'])
add('RTKN','rhotekin; RTKN1 (literature usage)','其他肿瘤细胞敲低实验','肺腺癌A549/SPC-A-1',
    'RTKN敲低研究报道肺癌细胞生长、侵袭和黏附减少，并观察细胞周期及复制相关蛋白变化。',
    '骨肉瘤检索中常出现的是RTKN2论文，属于不同基因；不能算RTKN的骨肉瘤支持。',
    '虽然有肺癌功能文献，但本项目2来源同向、3反向，不作为广泛适用marker优先。',
    [('Inhibition of rhotekin exhibits antitumor effects in lung cancer cells','PMCID PMC4811395','https://pmc.ncbi.nlm.nih.gov/articles/PMC4811395/')],
    ['"RTKN" "osteosarcoma"','"RTKN" cancer knockdown -RTKN2'])
add('TAF10','TATA-box binding protein associated factor10','其他肿瘤细胞和移植瘤机制实验','乳腺、肺、肝癌相关模型；MYC调控',
    'Z363研究通过TRIP12/TAF10/MYC关系解释抑瘤作用，包含细胞干预和小鼠肿瘤实验。',
    '没有把该机制定位到骨肉瘤CAF；广泛转录调控也不等于细胞类型特异marker。',
    '机制参考充分于所用模型，但本数据2来源同向、3反向，不提升为稳健跨来源marker。',
    [('Small molecule Z363 co-regulates TAF10 and MYC via the E3 ligase TRIP12 to suppress tumour growth','DOI10.1002/ctm2.1153; PMCID PMC9839843','https://pmc.ncbi.nlm.nih.gov/articles/PMC9839843/')],
    ['"TAF10" cancer fibroblast','"Small molecule Z363" TAF10 MYC'])
add('TOR1AIP1','LAP1; LAP1B/LAP1C isoforms','其他肿瘤体内外迁移机制；癌种间作用不同','黑色素瘤核力学、前列腺癌',
    '黑色素瘤研究显示短亚型LAP1C支持核膜适应、受限空间迁移和侵袭；前列腺癌研究则报道其稳定p53并抑制增殖侵袭。',
    '亚型和肿瘤背景重要，基因级表达不能直接等同LAP1C功能；非骨肉瘤CAF验证。',
    '具有明确迁移力学文献价值，但本项目2来源同向、3反向，应作为机制对照而非通用marker首选。',
    [('LAP1 supports nuclear adaptability during constrained melanoma cell migration and invasion','PMID36624187; DOI10.1038/s41556-022-01042-3','https://pubmed.ncbi.nlm.nih.gov/36624187/'),('TOR1AIP1 interacts with p53 to enhance cell cycle dysregulation in prostate cancer progression','DOI10.1007/s11010-025-05276-1; PMCID PMC12263815','https://pmc.ncbi.nlm.nih.gov/articles/PMC12263815/')],
    ['"TOR1AIP1" cancer metastasis','"LAP1 supports nuclear adaptability"'])

direction=pd.read_csv(O/'tables/MARKER_DIRECTION_SUMMARY.tsv',sep='\t');assert set(direction.symbol)==set(r['gene'] for r in rows)
frozen=pd.read_csv(B if False else P/'phase82_allCAF_crossdataset_markers_20261012/tables/MARKER_EVIDENCE_REGISTRY.tsv',sep='\t')
table=[];bib=[];queries=[]
for row in rows:
    d=direction.set_index('symbol').loc[row['gene']];s=frozen.set_index('symbol').loc[row['gene']]
    out={k:v for k,v in row.items() if k not in ['references','search_queries']}
    for k in ['expected_direction','same_sources','opposite_sources','zero_sources','same_source_list','opposite_source_list']:out[k]=d[k]
    out['source17_P']=s['SCPCP000017__P'];out['source23_P']=s['SCPCP000023__P'];out['source17_FDR']=s['SCPCP000017__FDR'];out['source23_FDR']=s['SCPCP000023__FDR']
    out['reference_titles']=' | '.join(r[0] for r in row['references']);out['reference_identifiers']=' | '.join(r[1] for r in row['references']);out['reference_urls']=' | '.join(r[2] for r in row['references']);out['searched_at']='2026-10-12';table.append(out)
    for title,ident,url in row['references']:bib.append(dict(gene=row['gene'],title=title,identifier=ident,url=url))
    for q in row['search_queries']:queries.append(dict(gene=row['gene'],query=q,date='2026-10-12'))
pd.DataFrame(table).to_csv(O/'tables/ALL14_LITERATURE_DIRECTION_INTEGRATION.tsv',sep='\t',index=False)
pd.DataFrame(bib).to_csv(O/'resources/LITERATURE_REFERENCES.tsv',sep='\t',index=False)
pd.DataFrame(queries).to_csv(O/'resources/LITERATURE_SEARCH_TERMS.tsv',sep='\t',index=False)
(O/'resources/LITERATURE_CURATED_RECORDS.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
report='''14个CAF候选marker的既往文献与跨来源方向联合评估
检索与核对日期 2026年10月12日

主要判断
NEBL具有较一致的部位表达方向（4来源同向、1零差）和既往黏附/迁移相关研究，可优先作为研究线索；但文献大量涉及LASP2短亚型或其他癌种，当前计数尚未分辨亚型，也没有证明骨肉瘤CAF中的同一机制。
AL096854.1已有骨肉瘤计算关联条目，可作为疾病相关lncRNA线索，但它只在原研究单因素候选表出现，不在最终六基因模型。GRID2主要有其他癌种预后/组织表达关联。RGS2-AS1和PEF1-AS1的直接肿瘤机制证据有限，不能把缺文献写成已证实的新颖性。
EBF1有骨肉瘤细胞及骨髓基质生物学研究，更适合作为骨谱系/机制参照；本项目存在跨来源和配对方向矛盾。RTKN、TAF10、TOR1AIP1虽有其他肿瘤实验研究，但本项目另外三个来源反向，文献不能替代数据一致性。

方法与证据口径
围绕14个主候选，检索正式符号、已核实别名、骨肉瘤/CAF/基质/转移及其他癌种。优先原始论文、PubMed/PMC、出版商页面和原始作者机构存档；基因数据库仅用于名称/亚型核对，不当作功能实验。部分论文以可访问摘要或正文相关段落为依据，记录了来源链接。本次是有范围的定向参考检索，不是穷尽性系统综述；未检出不等于文献不存在。
原始P<0.05的统计候选不因有无文献被增删；数据效应和P/FDR均保持。文献按实验对象与证据类型单独呈现。研究其他肿瘤细胞、正常骨髓基质或分子结构，不能直接称为骨肉瘤CAF机制验证。文献证据与本数据同向/反向均保留，不构造未经校准的总分。
'''
for row in table:
    report+='\n'+row['gene']+'\n'
    report+=f"本数据：{row['same_sources']}来源同向、{row['opposite_sources']}反向、{row['zero_sources']}零差；原两来源P={row['source17_P']:.6g}/{row['source23_P']:.6g}，FDR={row['source17_FDR']:.6g}/{row['source23_FDR']:.6g}。\n"
    for field,label in [('aliases','名称'),('evidence_type','文献类型'),('literature_context','模型/背景'),('literature_finding','已有研究'),('transfer_boundary','证据边界'),('combined_interpretation','综合判断')]:report+=label+'：'+str(row[field])+'\n'
    for ref in [r for r in bib if r['gene']==row['gene']]:report+=f"参考：{ref['title']} [{ref['identifier']}] {ref['url']}\n"
report+='''
必须避免的混淆
RTKN2不是RTKN；GRID2IP不是GRID2；LASP1不是NEBL/LASP2。NEBL的短亚型LASP2与长亚型不能仅凭基因级计数互换。CALHM家族其他成员的功能不能自动赋给CALHM5；PEF1或COL16A1不是PEF1-AS1；GACAT2论文的主要机制不能归给同文候选AL390957.1/RGS2-AS1。
AL096854.1在骨肉瘤文献的单因素候选表有记录，但不是最终六lncRNA模型成员。该文某些条目未达P<0.05仍被作者称独立预后等问题不在此复用；本表只引用可核对的候选条目。

下一步研究取向
将NEBL作为较优先的方向/黏附机制线索；AL096854.1和GRID2作为疾病关联或跨来源候选补充；RGS2-AS1保留探索价值。EBF1用于解释骨部位/间充质背景及反向结果。其余基因按本数据一致性和具体研究问题保留为对照或机制候选。这里不是临床/流式可用性排序，也没有宣布任何marker已完成独立验证。
'''
(O/'delivery/14个CAF候选marker_文献与方向联合评估_20261012.txt').write_text(report,encoding='utf-8')
qa=dict(genes_covered=len(rows),reference_records=len(bib),search_terms=len(queries),same_14_gene_set=True,original_statistics_unchanged=True,evidence_types_separated=True,alias_and_paralog_warnings_recorded=True,completed_at=datetime.now(timezone.utc).isoformat(),review_type='targeted narrative evidence review; not exhaustive systematic review')
(O/'QA/LITERATURE_REVIEW.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False));print(pd.DataFrame(table)[['gene','evidence_type','same_sources','opposite_sources']].to_string(index=False))
