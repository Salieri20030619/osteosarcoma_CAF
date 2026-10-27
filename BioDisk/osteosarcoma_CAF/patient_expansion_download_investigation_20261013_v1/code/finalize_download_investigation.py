
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,zipfile,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1';D=O/'delivery';D.mkdir(exist_ok=True)
def read(rel):return pd.read_csv(O/rel,sep='\t',low_memory=False)
def dump(data,path):path.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
required=['GSE234187_download','GSE234187_investigation','GSE_gene_mapping','OsteoCAR_download','OsteoCAR_extraction','primary_alias_resolution','all_human_investigation','S23_investigation','barcode_index']
for name in required:assert (O/'logs'/(name+'.exit_code')).read_text().strip()=='0',name
assert json.loads((O/'resources/S23_MARKER_INVESTIGATION_COMPLETE.json').read_text())['status']=='completed'
raw=[]
for rec in json.loads((O/'resources/OSTEOCAR_DOWNLOAD_COMPLETE.json').read_text())['files']:
 f=Path(rec['path']);assert f.exists() and f.stat().st_size==rec['bytes'] and rec['official_md5_verified'];raw.append(dict(resource='OsteoCAR',file=str(f),bytes=rec['bytes'],md5=rec['md5'],sha256=rec['sha256'],validation='official_MD5_verified'))
for rec in json.loads((O/'resources/GSE234187_DOWNLOAD_COMPLETE.json').read_text())['files']:
 f=Path(rec['file']);assert f.exists() and f.stat().st_size==rec['bytes'] and rec['gzip_integrity']=='PASS';raw.append(dict(resource='GSE234187',file=str(f),bytes=rec['bytes'],md5='',sha256=rec['sha256'],validation='expected_size_and_full_gzip_CRC_verified'))
pd.DataFrame(raw).to_csv(D/'DOWNLOADED_FILES_AND_HASHES.tsv',sep='\t',index=False)
qa_objects=[];frames=[]
for kind in ['patient_prim','patient_mets']:
 m=read(f'tables/OsteoCAR/{kind}_author_metadata.tsv.gz');s=read(f'tables/OsteoCAR/{kind}_cell_scores.tsv.gz');assert m.cell_id.is_unique and s.cell_id.is_unique and set(m.cell_id)==set(s.cell_id);assert len(m)==len(s)
 gs=[c for c in s if c.startswith('gene_log1pCP10k__')];assert len(gs)==14 and np.isfinite(s[gs].to_numpy()).all();assert s.full_library_counts.gt(0).all()
 fm=read(f'tables/OsteoCAR/{kind}_candidate_feature_map.tsv');assert len(fm)==14 and fm.feature_present_in_object.astype(str).str.lower().eq('true').all()
 qa_objects.append(dict(object=kind,cells=len(m),unique_cell_ids=True,score_ids_match=True,canonical_candidate_columns=14,all_candidate_object_features_present=True,zero_available_library_cells=0));frames.append(m)
m=pd.concat(frames);assert m.cell_id.is_unique and len(m)==374041
r=read('tables/OsteoCAR/ALL_HUMAN_SAMPLE_CAF_ROSTER.tsv');assert r.cells.sum()==374041 and r.author_CAF_cells.sum()==7633 and len(r)==61
n=r[r.author_source=='NCI-POB'];assert len(n)==19 and n.author_CAF_cells.sum()==1530 and n.author_CAF_cells.min()==11
cf=m.Ann_Level0.eq('Host')&m.Ann_Level2.eq('CAFs');assert cf.sum()==7633
site_conflict=json.loads((O/'resources/CONFIRMED_SAMPLE_SITE_CONFLICTS.json').read_text());assert len(site_conflict)==1
safe=read('tables/OsteoCAR/SJ_SITE_CONFLICT_EXCLUDED_CAPACITY.tsv');assert len(safe)==2
q=read('tables/GSE234187_cell_investigation.tsv.gz');assert len(q)==13684 and q.cell_id.is_unique;assert q.n_counts.gt(0).all();mapping=read('tables/GSE234187_candidate_stable_gene_mapping.tsv');assert len(mapping)==14 and mapping.assayed.sum()==12
for g in ['AL096854.1','RGS2-AS1']:assert q['gene_log1pCP10k__'+g].isna().all()
blocked=[dict(resource='NCC / CodeOcean9535428',status='not_downloaded',reason='HTTP403 from direct and proxy public routes'),dict(resource='HRA013407',status='not_downloaded',reason='TLS/connection failure from direct and proxy routes; permission status not determined'),dict(resource='Shen2026 aging cohort',status='accession_unresolved',reason='No verified GEO deposit or PubMed-GEO link located')]
pd.DataFrame(blocked).to_csv(D/'UNAVAILABLE_SOURCES.tsv',sep='\t',index=False)
source=read('tables/OsteoCAR/OBJECT_SOURCE_CAF_CAPACITY.tsv')
text='''CAF人源样本下载与调查报告
日期：2026-10-13

结论
可访问的人源表达数据已下载并完成本轮调查，总计5.236 GB：OsteoCAR原发与转移两个患者对象，以及GSE234187完整表达矩阵和配套文件。NCC/CodeOcean、HRA013407仍未取得原始数据；Shen2026未定位到可核验下载号。因此本轮不能表述为“所有候选来源全部下载成功”。
最有价值的扩充候选是OsteoCAR的NCI-POB肺样本：19条样本记录、123,255个细胞、1,530个作者标注CAF，每条记录11—283个CAF。按编号归并疑似同患者记录可形成17组，但正式患者对应表缺失，17与19均不能直接当作已核实的独立新增患者N。
原发对象的25条样本记录均具有已知来源名称与强条形码对应证据，不能算新增25人。S23可评估细胞数受身份注释、CNV和取材部位影响，不能仅把旧候选池中被排除的细胞全部放回。

一、数据与方法
下载对象：patient_prim.qs2（2,029,127,091字节）、patient_mets.qs2（3,080,547,578字节）；GSM7453693_matrix.mtx.gz（126,169,981字节），另有特征与条形码文件。下载的表达矩阵与含counts的Seurat处理对象保存在任务目录的downloads下。
OsteoCAR通过本机代理经仅绑定回环地址的临时SSH转发下载，使用分块断点续传，完成后分别比对作者提供的MD5。GSE234187经直连分块取得，核对文件大小、完整gzip CRC、矩阵维度和条形码。全部获取结束后已关闭本任务临时代理。
复用既有107份原始样本的全细胞对象，建立876,648个原始细胞记录的样本条形码索引。去重调查综合作者样本名称、提供方患者编号、条形码重合及选择基因的计数版本对照。
作者CAF定义按下载对象原注释Ann_Level0=Host且Ann_Level2=CAFs读取。Fibrogenic肿瘤细胞、Fibrogenic_TAMs及其他间充质类别分别保留。
GSE234187沿用项目已冻结的竞争性谱系基因集作描述性表达筛查，保留全部提供方细胞；未新增UMI/基因/线粒体阈值，未将该筛查等同于确认CAF身份。候选基因使用既有稳定Ensembl ID核对别名。
本轮没有进行新的推断性临床差异检验，没有将细胞数作为生物学样本量；原Phase86统计、8174细胞冻结基线、原分类和14候选名单保持不变。

二、OsteoCAR实际对象包含61条样本记录，而非清单里的83条
原发对象：152,770细胞，25条样本记录，4,096个作者CAF。
转移对象：221,271细胞，36条样本记录，3,537个作者CAF。
合计374,041细胞、61条样本记录、7,633个作者CAF。
作者GitHub清单列有83条人源记录，其中22条未出现在本次下载的两个患者对象中。部分属于其他转移部位或未保留的样本；未出现不能解释为没有CAF或没有表达，也不能计入本次有效样本量。
GEO只是存储平台，不是单一生物学队列。作者标为GEO的记录已进一步对照GSE152048与GSE162454，不能把它们当同一来源直接合并。

三、NCI-POB提供19条真实可读取的肺样本及1,530个作者CAF
19条记录全部有作者CAF，最低11、最高283个。与旧107样本的最佳条形码重合比例约0.82%—1.60%，没有发现类似已知重复文库中67%—100%的大面积重合。这支持其不是旧文库的简单重复拷贝，但不证明跨机构或跨时间点患者完全不重叠。
CCDI0057_16178与CCDI0101_16178b、CCDI0058_16183与CCDI0102_16183b按名称提示两组重复患者/重复取材。去掉b后按编号形成17组，此规则只是患者映射假设，未冒充作者明确给出的临床患者表。逐样本记录与17组假设清单均交付。
该批只有肺端。它适合在患者对应关系核实后开展肺转移CAF状态与候选表达的外部验证；不能直接与别的来源的骨端拼接，以“增加N”为由忽略来源和部位的混杂。
本轮未把NCI病例计入原主分析的独立样本量。

四、S23有恢复研究范围的空间，但作者CAF口径下肺端仍小
原始29份S23全细胞数据均在本地。既有全细胞表达复核中，6,314个细胞标记为Fibro，来自15名骨端、11名肺端患者；排除既有双联体标记后，这11名肺患者仍有表达候选，其中两名只有个位数细胞。这一层的谱系置信度为中，恶性身份未定。
旧阶段27候选规则为：Fibroblast表达优先，至少一项CNV无异常，且无相反异常调用，排除既有双联体。29份S23中阶段27候选2,153细胞，后续冻结候选1,862细胞，最终Fibro标签26细胞。该门槛复算一致，未发现实现与规则不符；但规则与RNA-CNV适用性本身需要评估。
在被指定为参考的S23细胞中，共13,975条记录；12,630条落入双方法异常、单方法异常或方法冲突，占90.38%。这不是经DNA真值验证的错误率，但表明不能无条件把这些RNA-CNV调用当作恶性身份定论。
本次作者原发对象中S23有192个CAF，其中190个能按同患者/部位条形码对应到既有复核记录；只有49个被旧表达复核标为Fibro，其余归入其他谱系或未定。对应细胞的旧CNV证据也未落入无异常类别。两种注释存在明显分歧，应保留交叉表，不能静默用任何一方覆盖另一方。
按作者CAF及作者部位记录，S23为骨10人、肺6人有至少1个CAF；按至少5细胞的描述性敏感性口径为骨8、肺3。
但同一标本SJOS001105_M1在本地ScPCA提供方元数据（SCPCS000784，submitter_id完全相同）写为Femur，作者表与对象写为lung；条形码重合98.68%。其Metastasis标签说明病程，不能单独证明取材肺部。该标本不是PDX。
排除这条未解决的部位冲突后，作者CAF覆盖为骨10、肺5人；至少5细胞口径为骨8、肺2。两种口径均没有两部位都具作者CAF的完整配对患者。已有临床清单里2对骨/肺病例不能自动变为2对可分析CAF样本。
因此，“可回查10名肺患者”不等于可以无条件新增10名可靠CAF患者。原有成纤维样候选、作者CAF、CNV状态和临床部位应分别记录，先解决身份与取材问题。

五、GSE234187是一名新来源的肺转移患者，但5,990个表达候选并非纯CAF
完整公开矩阵为33,694特征×13,684细胞、33,453,706非零项，零文库细胞为0。原文报告15,500细胞，与实际公开矩阵不一致，已记录版本/公开范围差异；没有补造缺失细胞。
沿用旧竞争性谱系评分得到5,990个Fibroblast表达候选。它们中IBSP检出76.31%、ALPL59.92%、RUNX2 36.08%，同时也有DCN、LUM、DPT等表达。这显示显著成骨谱系混合，不能把5,990直接称为确认CAF。
该来源能提供一个可读取的单病例验证对象；CAF身份与肿瘤/成纤维边界仍需进一步判定，单病例不能独立作骨肺组间推断检验。
与原107份样本的最大条形码重合约1.11%，未见明显旧文库复用；这不替代跨研究的临床患者去重。
按稳定基因ID核对，14候选中12个存在于该特征参考。CALHM5对应旧名FAM26E，PEF1-AS1对应RP11-73M7.6。AL096854.1和RGS2-AS1的稳定ID未在该公开矩阵特征表中找到，保持不可评估，不能填0或写成阴性。

六、作者对象能定位14候选，但不能直接混用旧矩阵与分数
OsteoCAR两个对象在对象层面均找到全部14个候选feature。PEF1-AS1对应AC114488.2，RGS2-AS1对应AL390957.1；该对应由项目既有稳定Ensembl ID及原始特征映射支持。
这是合并对象层的feature存在性，不等于已核实每个来源、每份样本的基因测量掩码。合并矩阵里的0值不能无条件区分未表达与原来源未测量。
25条原发记录的样本名称与条形码强烈支持既有来源对应；选择的11,000个基因-细胞计数比较中10,384项完全一致，但非零值存在差异，147,090个匹配细胞的作者nCount_RNA与旧来源全库总数也不完全相同。这说明处理版本不同，不能因文库来源对应就覆盖旧原始计数。
作者对象内诊断评分使用其可取得RNA矩阵的库总数；旧流程使用各原始来源的完整库总数。正式合并前需统一计数处理、基因测量状态及统计单元，本轮没有混合这些评分进行P值检验。

七、未完成下载的来源及接续条件
NCC/CodeOcean9535428：直连与代理仍返回HTTP403，未取得表达对象。其现有作者表10份样本、9个患者ID，不得计为10人。
HRA013407：本机、服务器、代理及官方旧入口均未建立有效下载连接；连接错误不能解释为数据不存在或已经需要某类权限。原文3人6份化疗前后配对数据及额外非配对验证仍未获取。
Shen2026：原文缺可核验GEO号，定向GEO与PubMed-GEO链接查询未找到存储对象，未下载。
未发送邮件或申请受控访问。剩余来源需可访问的原始文件/公开入口，或由数据持有人提供的获批入口，方能继续完成获取。

八、后续优先顺序
1. 首先补齐NCI-POB患者-标本-文库对应表，确认17组命名假设与原队列患者是否重叠，明确同患者重复取材。
2. 对S23做身份判定的一致性复核：作者Host/CAFs、旧表达注释、原始标志与CNV证据并列；解决SJOS001105_M1取材部位后再冻结可分析细胞与患者清单。
3. 对GSE234187区分成纤维与纤维/成骨肿瘤状态，保留单病例探索性定位。
4. 再决定新临床比较设计。NCI仅肺端，不能简单并入旧骨肺检验；既有14候选继续保留。所有后续有效检验仍按目标效应原始P<0.05筛选并报告FDR，不能为获得阳性而放宽细胞身份规则。

来源
OsteoCAR数据：https://doi.org/10.6084/m9.figshare.31029559
OsteoCAR作者表：https://github.com/kidcancerlab/Osteo_atlas/blob/main/allsample_details.csv
GSE234187：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234187
GSE234187原文：https://pmc.ncbi.nlm.nih.gov/articles/PMC10784140/
ScPCA S23：https://scpca.alexslemonade.org/projects/SCPCP000023
HRA013407原文：https://www.oaepublish.com/articles/cdr.2026.09

交付范围
下载的表达大文件留存在sgs10的patient_expansion_download_investigation_20261013_v1/downloads；报告包提供下载校验清单、全部逐样本与候选患者清单、原始证据交叉表、基因映射、当前细胞级调查表、脚本及关键验收记录。历史失败尝试保留于服务器logs，当前有效完成记录以FINAL_ACCEPTANCE.json及本报告为准。
'''
text+='\n附表A：实际对象与作者CAF覆盖\n'+source.to_string(index=False)+'\n'
text+='\n附表B：NCI-POB逐样本清单（不是已核实独立N）\n'+n[['author_sample','cells','author_CAF_cells','query_overlap_fraction','unverified_NCI_name_based_patient_key']].to_string(index=False)+'\n'
text+='\n附表C：排除取材冲突后的S23作者CAF患者覆盖\n'+safe.to_string(index=False)+'\n'
(D/'CAF_新增人源数据下载与调查报告_20261013.txt').write_text(text,encoding='utf-8-sig')
qa={'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'available_artifact_integrity':'PASS','all_requested_sources_downloaded':False,'downloaded_expression_objects':3,'downloaded_data_files':len(raw),'downloaded_bytes':sum(z['bytes'] for z in raw),'official_MD5_verified_objects':2,'GSE_full_gzip_and_matrix_checks':'PASS','object_checks':qa_objects,'NCI_source_spelling_normalized':True,'author_CAF_count_reproduced':7633,'NCI_lung_sample_records':19,'NCI_author_CAF_cells':1530,'NCI_patient_groups_from_names_unverified':17,'NCI_independent_new_patient_count':None,'GSE234187_gene_ids_verified':12,'GSE234187_unavailable_candidate_ids_preserved_as_NA':True,'baseline_unchanged':True,'new_inferential_tests':0,'unavailable_sources':blocked}
dump(qa,D/'FINAL_ACCEPTANCE.json')
status={'project':str(P),'task_dir':str(O),'status':'available_sources_downloaded_and_investigated_external_access_pending','updated_at_utc':qa['checked_at_utc'],'all_requested_sources_downloaded':False,'raw_data_dir':str(O/'downloads'),'report':str(D/'CAF_新增人源数据下载与调查报告_20261013.txt'),'baseline_phase':'phase86_functional_bone_lung_20261012','baseline_changed':False,'formal_new_patients_added_to_baseline':0,'NCI_lung_records':19,'NCI_CAF_cells':1530,'NCI_verified_independent_new_patients':None,'unavailable_sources':blocked,'next_action':'Confirm NCI patient crosswalk, resolve S23 anatomical and identity conflicts, then freeze expansion analysis design.'}
dump(status,O/'TASK_STATUS.json');pointer=P/'CAF_PATIENT_EXPANSION_CURRENT.json'
if pointer.exists():shutil.copy2(pointer,O/'resources/previous_patient_expansion_pointer.json')
dump(status,pointer)
files=[f for sub in ['tables','code'] for f in sorted((O/sub).rglob('*')) if f.is_file() and '__pycache__' not in str(f)]
files += [D/'CAF_新增人源数据下载与调查报告_20261013.txt',D/'DOWNLOADED_FILES_AND_HASHES.tsv',D/'UNAVAILABLE_SOURCES.tsv',D/'FINAL_ACCEPTANCE.json',O/'TASK_STATUS.json']
receipts=['GSE234187_DOWNLOAD_COMPLETE.json','OSTEOCAR_DOWNLOAD_COMPLETE.json','GSE234187_INVESTIGATION_COMPLETE.json','OSTEOCAR_PATIENT_INVESTIGATION_COMPLETE.json','MEASUREMENT_AND_COUNT_VERSION_BOUNDARIES.json','CANDIDATE_ALIAS_RESOLUTION.json','CONFIRMED_SAMPLE_SITE_CONFLICTS.json','S23_RETENTION_COMPLETE.json','S23_MARKER_INVESTIGATION_COMPLETE.json','GSE234187_GENE_ID_RESOLUTION_COMPLETE.json','SOURCE_NORMALIZATION_CORRECTION.json','PROXY_ROUTE_CHECK.json','ALTERNATE_PUBLIC_ROUTES.json']
files += [O/'resources'/n for n in receipts if (O/'resources'/n).exists()]
manifest=[dict(path=str(f.relative_to(O)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];dump(manifest,D/'PACKAGE_MANIFEST.json')
archive=D/'CAF_human_data_download_investigation_20261013.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for f in files+[D/'PACKAGE_MANIFEST.json']:z.write(f,str(f.relative_to(O)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
receipt={'archive':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'report':str(D/'CAF_新增人源数据下载与调查报告_20261013.txt'),'data_bytes':qa['downloaded_bytes'],'files':len(files)+1}
dump(receipt,D/'DELIVERY_RECEIPT.json');print(json.dumps(receipt,ensure_ascii=False))
