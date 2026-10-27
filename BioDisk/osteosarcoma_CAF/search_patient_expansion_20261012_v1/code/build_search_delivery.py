
from pathlib import Path
import pandas as pd,json,hashlib,zipfile,datetime,gzip
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'search_patient_expansion_20261012_v1'
rows=[
['P1','GSE234187 / GSM7453693','新增公开单细胞候选','1名患者、1份肺转移灶','肺转移','现有九来源未包含此GSE；患者跨库身份尚未有完整对照','GEO矩阵端点HTTP200，字节范围HTTP206；已核对条形码与矩阵维度','原文报告CAF；实际CAF数量及身份待处理','优先下载完整矩阵并作入组前核查，单来源n=1不单独作推断检验','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234187'],
['P1','SCPCP000023','既有原始数据中恢复可评估患者','肺12人中当前Fibro仅2人；另10人有原始对象','骨/肺；原清单有2对','属于现有队列，不是新发现队列；可增加有效N取决于真实CAF能否恢复','相关全细胞对象全部在服务器','9名肺患者无冻结候选细胞；另1名有候选但无Fibro','按样本追查上游QC、谱系与CNV排除原因；不直接放宽阈值','https://scpca.alexslemonade.org/projects/SCPCP000023'],
['P2','OsteoCAR NCI-POB subset','潜在较大新增肺样本来源','26条样本记录，其中19条明确标为lung；独立患者数未核实','肺/其他转移部位','作者标签疑似17个肺患者标识基础，但只是命名推断，不能作N；全atlas与旧数据大量重叠','作者GitHub最新清单可下载；Figshare API本次仍HTTP403，表达对象未取得','尚无当前可核验CAF数','优先取得患者-标本-文库对照及可下载表达对象，再与全部现有来源去重','https://github.com/kidcancerlab/Osteo_atlas/blob/main/allsample_details.csv'],
['P2','HRA013407','新增治疗前后单细胞候选','原文3人6份配对样本；另述5份非配对验证样本','化疗前后肿瘤，未证实新增肺灶','与现有九来源无同一accession；新增患者和验证数据是否完整入库待核','论文给出GSA-Human号，本次门户连接未成功；开放权限与文件完整性未核','原文报告CAF扩张','更适合耐药/治疗反应支线，不能填充骨肺比较中的肺组','https://www.oaepublish.com/articles/cdr.2026.09'],
['P3','Zheng2024 NCC / CodeOcean9535428','未取得矩阵的既有候选资源','作者清单10份样本、9个唯一Patient ID；Patient4两次取材','9份Primary、1份Recurrent；清单未给肺灶','不得将整合的GSE152048/GSE162454再算新增；自有病例待跨库核对','作者GitHub患者清单HTTP200；CodeOcean本次HTTP403','原文含间充质细胞，细胞身份需核','补原发/治疗前后方向，当前不能解决肺端不足','https://github.com/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas'],
['P3','Shen2026 aging study','文献候选','12例，9年轻/3老年','初治；未找到新增肺灶证据','不能因二次文献称public就计入','正文仅称存入GEO，未给出可核验GSE；本次未取得矩阵','待核','优先级低于明确肺样本和可下载对象','https://link.springer.com/article/10.1007/s13402-026-01162-0'],
['Existing','GSE293065','已有空间验证数据','19名患者、26份样本；肺12份来自10人','骨/肺/其他复发转移','Phase35/36已经使用，新增患者计0','已有本地数据与患者清单','空间混合信号，非单细胞CAF','可用于本次功能线的已有外部验证，不作为新增单细胞队列','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293065'],
['Existing','GSE270231 / phs003569.v1.p1','同一队列的双重入口','5名患者、7份肺样本','肺','GSE270231已在现有九来源；dbGaP官方直接链接该GSE、同为5人','已有原始数据','当前2名患者有Fibro审核细胞','不可把phs号当成新增5人；可复查未进入CAF的既有病例','https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003569.v1.p1'],
['Existing','SCPCP000018 / GSE152048 / GSE318840','既有记录性来源','本次未增加患者','各有原骨肺记录','全部已有；继续遵守record-only定位','原始对象存在','既有CAF覆盖/部位证据不足','只记录恢复可行性，不自动纳入主决策或修改原14候选排序','existing_project_manifests'],
['Exclude_newN','GSE319124 / SJOS030605_D1_DYE7028','同患者重复入口及模型','59样本跨癌种；OS患者样本SJOS030605已存在','该患者既有胸壁/肺记录','GSM9511156写明Patient，但SJOS030605已在SCPCP000023，不能新增计数','GEO样本元数据可取得','其余OS条目主要模型','区分患者样本和培养/PDX，不能将59样本计为59名OS患者','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9511156'],
['Exclude','GSE200529 / GSE307240 / GSE270504','PDX或PDX来源细胞','不计临床新增N','模型','即使物种字段Homo sapiens也非直接患者CAF队列','官方设计已核查','不能当作患者原位CAF','不纳入临床患者扩充','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200529'],
['Exclude','GSE198896 / GSE216280','CD45分选或PBMC','不计CAF样本扩充','免疫细胞/血液','不解决组织CAF样本量','GSE198896官方设计为CD45+','CAF被分选排除或不在组织','保留为免疫参考','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198896'],
['Exclude','GSE290255','细胞系LRRC15实验','不计患者N','细胞系','不是患者肿瘤','官方设计已核查','仅机制参考','不作为患者补充','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290255'],
['Exclude','GSE283885 / GSE217792','肺纤维化或骨修复参考','不计OS患者N','非骨肉瘤参考','OsteoCAR引用不代表新增OS病例','官方GSE283885记录已核查','只适合组织参考','不混入骨肉瘤骨肺临床组','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283885'],
['Existing_other_endpoint','GSE237070','已有重建/核查来源','原发2人，癌旁4人；非肺转移设计','骨原发/癌旁及淋巴结','已做Phase37/38来源工作','旧记录及原始路线可复核','不新增肺端','不把淋巴结转移或癌旁当肺转移','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237070']]
columns=['priority','resource','route','patient_sample_count','site_or_endpoint','overlap_status','access_status','CAF_relevance','next_action','primary_reference']
reg=pd.DataFrame(rows,columns=columns);reg.to_csv(O/'tables/EXPANSION_CANDIDATE_REGISTRY.tsv',sep='\t',index=False)
existing=pd.read_csv(O/'tables/existing_patient_expansion_roster.tsv',sep='\t');s23=existing[(existing.dataset=='SCPCP000023')&(existing.site=='Lung')];assert len(s23)==12 and s23.source_counts_exist.all() and s23.Fibro_cells.gt(0).sum()==2
nci=pd.read_csv(O/'tables/OsteoCAR_NCI_POB_candidate_records.tsv',sep='\t');assert len(nci)==26 and nci.location.eq('lung').sum()==19
ncc=pd.read_csv(O/'sources/NCC_patient_metadata.txt',sep='\t');assert len(ncc)==10 and ncc.ID.nunique()==9
probe=json.loads((O/'MATRIX_HEADER_PROBE.json').read_text());assert probe['matrix_shape_and_nonzeros']=='33694 13684 33453706'
with gzip.open(O/'sources/GSE234187_barcodes.gz','rt') as f:barcodes=f.read().splitlines()
assert len(barcodes)==len(set(barcodes))==13684
series=pd.read_csv(O/'tables/GEO_series_screen.tsv',sep='\t');assert len(series)==46
text='''CAF项目补充患者样本可行性复查
日期：2026-10-12（美国东部）；本轮为定向检索、元数据核查和去重，不是新增正式分析。

结论
有补充空间，但目前没有核实到一个可直接并入的、大规模、无重叠的人源肺转移CAF队列。
最现实的两条路线：一是GSE234187的一名肺转移患者作为新增单细胞候选；二是回查S23现有原始数据中未进入Fibro分析的10名肺患者。更大增量线索为OsteoCAR的NCI-POB部分，其19条肺样本记录尚不能算作19名独立新患者。

1. 新增公开肺转移候选：GSE234187
原始论文Expression of PD-1 mitigates phagocytic activities TAM in osteosarcoma（PMID38223729；DOI10.1016/j.heliyon.2023.e23498）明确描述1名患者的1份肺转移灶，并报告CAF群。GEO为1个GSM（GSM7453693），未列入当前九来源清单。
本次取得完整条形码与特征文件，并用HTTP Range核查矩阵头：33,694特征×13,684细胞条目，33,453,706个非零项。完整压缩矩阵端点HTTP200，大小126,169,981字节（约126MB）；范围请求HTTP206。
本次只下载了小型配套文件及矩阵前64KiB片段，未下载完整矩阵，也未完成细胞QC/CAF注释。原文写15,500细胞，与当前公开矩阵13,684条形码不一致，须在正式入组前记录并核对版本；以实际可取得矩阵为处理对象，不填补缺失细胞。
这是1名患者候选，不能把细胞数当样本量，也不能单独建立骨肺两组的推断性验证。跨研究身份、取材/治疗信息和CAF数需完成入组核查。新增“可分析CAF患者”目前仍为0，1是待处理候选人数。
链接：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234187
原文：https://pmc.ncbi.nlm.nih.gov/articles/PMC10784140/

2. 最值得先做的现有数据恢复：SCPCP000023
项目临床清单有12名肺患者，全部相关全细胞对象在服务器；当前只有2人进入Fibro层。因此有10名肺患者可逐样本回查，其中9人未进入冻结候选池，另1人有候选细胞但没有Fibro标签。
这10人是“待核查患者”，不是保证能恢复的10名CAF阳性患者。真实肿瘤、谱系不明确、无可靠正常成纤维细胞或QC不合格的样本不能为凑N直接放回。
骨端同理：15人中当前Fibro2人，另13人有待复查空间。两端患者可能重叠，不能把10与13简单相加为独立新患者。
原临床清单存在2对骨/肺患者SJOS013768、SJOS063833，但当前Fibro比较中没有有效完整配对。若两部位均能获得可靠CAF，才可能增加配对证据。
S17原清单的肺患者已全部进入Fibro层，当前主要不是遗漏肺患者；非配对只有3名是因为另2名按配对设计单列，不应误说还漏了2名。
链接：https://scpca.alexslemonade.org/projects/SCPCP000023

3. 较大但未打通的增量：OsteoCAR NCI-POB
本次重新取得作者GitHub当前样本表：人源样本记录83条；NCI-POB占26条，19条明确标为lung，其余为肝、主动脉旁、胸壁、纵隔等，不能全部放入肺组。
19条肺记录的命名可拆出17个数字标识基础，其中16178/16178b及16183/16183b疑似同患者重复取材。这只是标签推断，不是已核实的患者映射，故独立新增患者数仍记unknown。
整套OsteoCAR不是新队列：旧审计已确认39条人源记录与发现数据同患者/样本，另7条为已知或候选重叠；其余37条待完整crosswalk。GEO旧数据、St Jude、NCH和NCI各部分必须分别去重。
作者宣告的表达对象托管在Figshare10.6084/m9.figshare.31029559。本次服务器API仍HTTP403，未取得新的表达对象。dbGaP条目不能自动等同可公开取得的全部矩阵或独立患者数。
优先需求：患者-部位-取材-文库映射、可读取的人源表达对象/细胞注释、与SCPCP及GEO的对应关系。未发送邮件、未申请受控数据或作任何对外提交。
作者表：https://github.com/kidcancerlab/Osteo_atlas/blob/main/allsample_details.csv
论文：https://doi.org/10.64898/2026.01.13.698472
数据：https://doi.org/10.6084/m9.figshare.31029559

4. 可以扩充其他临床方向，但不直接补肺组
HRA013407：2026年化疗前后单细胞研究，原始发现为3名患者6份配对样本；论文另述5份非配对验证样本。适合耐药、化疗前后CAF研究。论文给出GSA-Human入口，本次门户连接未成功，开放/受控状态与实际文件集合未核实；不能把文献所有样本直接计入已取得数据。GSE317512是其bulk相关队列，不是新增单细胞CAF患者。
论文：https://www.oaepublish.com/articles/cdr.2026.09
NCC/Zheng2024：作者现有患者表10份样本，但仅9个唯一Patient ID，Patient4对应BS4与BS14两次取材；清单9份Primary、1份Recurrent，未列肺转移取材。可以补原发/治疗方向，不能误计为10个独立患者。CodeOcean9535428本次仍HTTP403，表达数据未获得。
元数据：https://github.com/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas
Shen2026年龄队列：论文描述12例初治患者（9年轻、3老年），但正文只泛称GEO存储，未找到可核验GSE；作为低优先级文献候选保留，不纳入新增人数。

5. 去重和排除结果
GSE293065空间数据虽于2026年发表，但项目Phase35/36已使用：19患者26样本，肺12份来自10人。新增队列计0；可作为当前功能问题的已有空间验证，不能当新增单细胞患者。
phs003569.v1.p1官方明确链接GSE270231，且同为5名患者；GSE270231的7份肺样本已经在项目中，不可再加5人。
GSE319124于2026年新增公开，但59样本跨多癌种与模型；其中GSM9511156为SJOS030605_D1患者样本，与已有SCPCP000023患者SJOS030605重叠，不增加独立患者。其余OS条目主要为衍生模型。
GSE200529、GSE307240、GSE270504为PDX/衍生培养体系；GSE290255为细胞系；GSE198896为CD45+分选，GSE216280为PBMC；不能用于扩大人源原位CAF样本量。
GSE283885是肺纤维化/非纤维化参考，GSE217792为骨修复参考；被OsteoCAR引用不表示新增OS病例。
GSE237070是既有原发/癌旁与淋巴结方向，不补肺灶。SCPCP000018、GSE152048、GSE318840仍为既有记录性来源，本次不自动恢复到主分析或更改14候选排序。

6. 执行优先级
第一步：GSE234187完整矩阵的入组前核查，同时按既有样本名单追查S23肺端10名患者和2对骨肺病例。前者是小幅外部增量，后者是同来源内恢复有效N的机会。
第二步：打通OsteoCAR NCI-POB的人源表达对象与患者crosswalk，确认真实独立增量。
第三步：如研究目标转向治疗/耐药，再考虑HRA013407和NCC。它们不能代替骨肺分组的肺病例。
新增来源若只有肺组，不能把“来源A的骨”与“来源B的肺”直接合并检验，造成来源与部位完全混杂。优先保留来源内骨肺比较；肺单臂资源可用于CAF身份、表达、功能方向描述及适用的外部验证。
本次没有修改原细胞池、阈值、标签、统计结果或项目当前状态，未启动新增全量分析。

检索与证据完整性
完成项目原九来源、旧外部资源登记、ScPCA当前页面、GEO定向检索和原论文/作者元数据核对。GEO检索返回125个混合记录，其中46个Series已逐项获取标题与摘要初筛；其余为Sample/其他记录，不能当125队列。此次为定向检索，非穷尽系统综述。
附表含46系列检索结果、候选优先级、现有患者逐人纳入情况与原始对象路径、NCI-POB逐样本待去重清单。公开门户不可达与权限未核实分别记录，不把连接失败解释为数据不存在。
'''
(O/'PATIENT_EXPANSION_REVIEW_ZH.txt').write_text(text,encoding='utf-8-sig')
source_urls={'GSE234187':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234187','paper234187':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10784140/fullTextXML','OsteoCAR':'https://raw.githubusercontent.com/kidcancerlab/Osteo_atlas/main/allsample_details.csv','NCC':'https://raw.githubusercontent.com/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/main/Patient%20metadata','HRA013407_article':'https://www.oaepublish.com/articles/cdr.2026.09'}
record={'completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'search_and_feasibility_verified','GEO_series_screened':46,'new_public_lung_patient_candidate':1,'new_admitted_CAF_patients':0,'S23_existing_lung_patients_to_review':10,'S23_existing_potential_Bone_Lung_pairs':2,'NCI_POB_lung_sample_records':19,'NCI_POB_independent_new_patients':None,'NCC_samples':10,'NCC_unique_patient_IDs':9,'full_GSE234187_matrix_downloaded':False,'partial_header_verified':True,'baseline_changed':False,'source_urls':source_urls}
(O/'SEARCH_ACCEPTANCE.json').write_text(json.dumps(record,ensure_ascii=False,indent=2))
source_manifest=[dict(file=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted((O/'sources').glob('*')) if f.is_file()];(O/'SOURCE_EVIDENCE_MANIFEST.json').write_text(json.dumps(source_manifest,indent=2))
files=[O/'PATIENT_EXPANSION_REVIEW_ZH.txt',O/'SEARCH_ACCEPTANCE.json',O/'SOURCE_EVIDENCE_MANIFEST.json',O/'ACCESS_CHECKS.json',O/'MATRIX_HEADER_PROBE.json']+list((O/'tables').glob('*'))+list((O/'code').glob('*.py'))
manifest=[dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];(O/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
archive=O/'CAF_patient_expansion_search_20261012.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files+[O/'PACKAGE_MANIFEST.json']:z.write(f,str(f.relative_to(O)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'report':str(O/'PATIENT_EXPANSION_REVIEW_ZH.txt'),'registry_rows':len(reg)},ensure_ascii=False))
