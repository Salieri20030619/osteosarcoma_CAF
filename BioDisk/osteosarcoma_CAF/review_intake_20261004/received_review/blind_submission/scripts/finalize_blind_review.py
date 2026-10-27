import pathlib,pandas as pd,numpy as np,json,hashlib,datetime,platform,shutil
P=pathlib.Path(__file__).parent;O=P/'review_output';F=O/'blind_submission'
r=pd.read_csv(F/'reviewer_response_template.tsv',sep='\t');d=pd.read_csv(F/'cluster_review_details.tsv',sep='\t');st=pd.read_csv(O/'cluster_structure_review.tsv',sep='\t');q=pd.read_csv(O/'cluster_expression_qc.tsv',sep='\t');src=pd.read_csv(O/'source_coverage_summary.tsv',sep='\t');mix=d[d.lineage_code.str.startswith('Mix')];alt=pd.read_csv(O/'alternative_reference_sensitivity.tsv',sep='\t');g=pd.read_csv(O/'selected_gene_evidence.tsv.gz',sep='\t');idx=pd.read_csv(F/'supporting_gene_index.tsv',sep='\t');check=idx.merge(g,on=['cluster_id','gene','feature_id'],how='left',validate='one_to_one');assert check.detection_fraction.notna().all();assert (check.detection_fraction>0).all();assert ((g.detection_fraction>=0)&(g.detection_fraction<=1)).all();assert (g.positive_patients<=g.N_patients_matched).all();
meta=pd.read_csv(P/'sample_technical_metadata.tsv',sep='\t');raw=pd.read_csv(P/'RAW_COUNT_MANIFEST.tsv',sep='\t');raw['sample_id']=raw.relative_path.str.extract(r'(S\d+)');assert set(raw.sample_id)==set(meta.sample_id);assert raw.set_index('sample_id').cells.sort_index().equals(meta.set_index('sample_id').cells.sort_index());s=pd.read_pickle(O/'cell_structure.pkl');roster=pd.read_csv(P/'anonymous_cell_cluster_roster.tsv.gz',sep='\t');assert s.set_index('cell_id')[['source_id','sample_id','patient_id','r10']].sort_index().equals(roster.set_index('cell_id')[['source_id','sample_id','patient_id','r10']].sort_index())
# Cross-resolution summaries are row counts, not evidence of independently replicated populations.
source_lines=[]
for x in src.itertuples():source_lines.append(f'{x.source_id}: {int(x.cells):,}细胞，{int(x.clusters)}簇；CopyKAT覆盖{x.copykat_coverage:.1%}；未定谱系{int(x.unresolved_cells):,}细胞。')
fb=d[d.lineage_code=='Fibro']; fb_lines=[]
for x in fb.itertuples():fb_lines.append(f'{x.cluster_id}（{x.cells:,}细胞；{x.malignancy_assessment}；CopyKAT覆盖{x.copykat_coverage:.1%}）')
altlines=[]
for sid,a in alt.groupby('sample_id'):
 eligible=a[a.both_reported&(a.immune_endothelial_reference.abs()>=10)&(a.lymphoid_reference.abs()>=10)];opp=(eligible.direction_agrees_both_reported==0).sum();altlines.append(f'{sid}: 双参考均报告且双方|净区域长度|≥10 Mb的{len(eligible)}条染色体记录，方向相反{opp}条。这些记录不是独立生物学样本。')
report=f'''CAF重注释独立盲审报告
完成日期：2026-10-03
证据包：CAF重注释_盲审材料_20261002.zip（phase53）
本报告在读取新注释及揭盲对应表之前完成。项目判读证据仅来自本次盲审包，未读取其他项目记忆、历史注释或项目目录。外部查询仅核查通用谱系marker语义，未用于识别匿名来源或导入既有标签。新注释比较须另建文件，不覆盖本提交。

1. 307个簇全部保留；成纤维样表达不等同于已确认CAF
已填写模板全部307行和6项判读字段，保留原source_id、cluster_id、cells、patients及行序，覆盖876,648个细胞。判为成纤维细胞样15簇、91,012细胞；另有14簇、31,356细胞处于成纤维/成骨样间充质边界。明确保留未定谱系31簇、116,202细胞（{116202/876648:.1%}），混合表达/混合程序{len(mix)}簇、{int(mix.cells.sum()):,}细胞。上述是簇标签覆盖的细胞数，不代表簇内每个细胞都已得到独立验证。
成骨样48簇、237,184细胞；软骨样4簇、18,725细胞；骨/软骨双分化1簇、757细胞。免疫、血管壁、内皮、肌源性及上皮群单列，完整分布见lineage_summary.tsv。CAF正式身份需结合非恶性基质身份及组织上下文；本次未将全部COL1A1阳性细胞归为CAF。

2. D01_K017等基质表达群存在恶性候选证据
D01_R10_K017的COL3A1/LUM/DCN/POSTN支持成纤维样基质表达。CopyKAT覆盖2,666/2,720细胞，同样本表达参照比较可用6患者，幅度比的患者中位数2.93，4/6患者≥2倍；57项样本-簇-染色体方向比较中98.2%一致。该簇仍保留成纤维样表达谱系，但恶性身份列记为CNV支持恶性候选，不直接纳入确定非恶性CAF。
D06_R10_K018、D08_R10_K017、D09_R10_K002等间充质边界簇也需优先核查。其中部分簇主要由单一患者贡献，不能把簇内细胞数或少量散在患者的存在当作跨患者重复。
所有成纤维样候选如下，完整效应、覆盖和备注见fibroblast_boundary_review.tsv：
'''+ '\n'.join(fb_lines)+f'''

3. 状态程序与混合表达已单列
增殖并非一个谱系：D01_K012为成骨样增殖群，D01_K020为T细胞增殖群，D01_K024为髓系增殖群；D03_K017及D07_K042为内皮增殖群。IFN、即时早期应激、线粒体或核糖体主导簇分别保留状态标记。
D01_K039的CD3D检出86.6%、C1QA检出99.0%，由边际检出率可推得至少85.5%的条形码两者同检出。D07_K027的PECAM1与SFTPB对应下界为95.0%。D04_K007和D09_K010也存在广泛髓系/基质同检出。该数学下界确认的是同一条形码中两种转录本至少共同出现的比例，不能区分双细胞、细胞复合体、环境RNA、吞噬/融合或真实生物程序。
B/浆细胞/pDC程序并存者保留混合或未细分。JCHAIN、GZMB、COL1A1、ACTA2、FLT4、MBP等均不作为单基因身份判定。肌源性、肾小管/集合管样及肺上皮方向使用成组证据，未因样本含肿瘤而自动认定恶性。

4. CNV覆盖有明显来源差异，替代参考存在方向敏感性
有效CopyKAT连续指标覆盖568,754/876,648细胞（{568754/876648:.1%}），覆盖83/107样本；其余细胞保留缺失。D05覆盖仅17.4%，许多群每细胞汇总计数不足1,000，不能以其低CNV幅度宣布正常。
'''+ '\n'.join(source_lines)+'''
替代参考单独核对，未与provider幅度混合：
'''+ '\n'.join(altlines)+'''
CNV与表达谱系分栏。两种RNA推断方法存在相关误差，其方向一致不等于独立DNA验证。未把区域分箱数、染色体数或技术样本数当作生物学样本量；替代参考变化及不一致记录完整保留。

5. 输入完整性通过；D08存在可解释但尚未原始复算的分母差异
PACKAGE_FILE_MANIFEST列出的302个文件全部存在，字节数和SHA256均匹配。876,648个cell_id唯一；总roster与9来源cell_structure的细胞、来源、样本、患者、主簇一一一致；107份样本计数与技术元数据、RAW_COUNT_MANIFEST一致。各来源主簇细胞总量、HDF5列序及矩阵feature维度均核对。
8个来源的簇矩阵列和与full_library精确一致。D08全簇矩阵计数合计294,566,790，比full_library少122,334（约0.0415%）；技术说明给出32,732公共feature和37,407并集feature。差额方向与公共feature汇总/完整文库口径一致，但本包没有逐细胞H5AD，本次不能独立确认差额逐基因来源，也不能复验原始计数数组不变这一上游声明。未修改原始输入或补造缺失计数。
107份原始H5AD仅有manifest，合计3,605,632,476字节。本次交付是包内证据的完整盲审；未重跑原始单细胞聚类、双细胞识别、RNA校正或CNV，也未下载服务器上的原始矩阵。

6. 方法与可复核口径
逐来源按模板307簇判读。先阅读各簇前22个描述marker，再检查HDF5中预先指定的多谱系面板、状态面板及相关展示marker，累计导出97,032条feature-簇证据。全部可比feature矩阵参与维度、计数和文库核查；没有声称逐条人工阅读所有基因。相同gene_symbol对应多个feature时保留原feature_id；展示支持基因选取检出率乘以非负患者效应加0.1后最大的feature，不相加重复符号的检出率。
谱系初判由逐簇证据阅读后记录于reviewer_decisions.json；面板排序仅帮助查阅，最终标签不是最大模块分数直接分类。高/中/低指谱系证据的定性把握：高为多个相容特异marker支持，中为方向明确但计数/状态/亚型有混淆，低为身份边界、混合或证据不足；不是概率或经验证准确率。
支持基因字段报告簇检出率、包内患者级平均log2差异及正向患者/N_patients_matched。没有新增P值、FDR或显著marker声称，也未以细胞作为独立重复。本次是描述性身份审查，不具备预设独立检验，因此未套用P<0.05筛选。来源间患者不能保证独立，不汇总为跨来源独立患者队列。
分辨率审查使用r0.5/r1.0/r2.0条形码交叉表，记录主归属比例、反向比例及占主簇≥10%的高分辨率子群数。10%是展示分群结构的描述阈值，不删除细胞；主归属不足80%记边界复核。UMAP只作显示，不用于证明独立谱系或恶性。
簇线粒体比例=MT-符号feature原始计数/簇矩阵总计数；核糖体比例用RPL/RPS符号前缀，包含命中前缀的假基因。均为簇汇总比例，不是细胞中位数。计数/细胞数也不替代逐细胞QC。
CopyKAT描述比较使用同样本表达明确的T/Treg/NK/B/肥大/血管或淋巴内皮群作表达参照。参照需≥30个有效细胞且来自≥2个簇；被比较样本-簇需≥20有效细胞。计算该簇细胞mean_abs_native_signal中位数/同样本参照细胞中位数，同一患者有多样本时再取比值中位数，最后对患者等权取中位数。20/30细胞仅为本次CNV摘要可评估条件，不是原细胞筛选。全部覆盖计数仍保留；参照由表达定义而非经验证二倍体，不能用比值接近1确认正常。
区域CNV按sample_id和cell_id关联主簇。原生i6 state<3为loss、state>3为gain、3为中性；对每区域组、染色体、方向合并区间，长度=end-start+1，以组内细胞数加权到样本-簇。未报告的区域或组保留不可评估，不填作中性。provider保留每条染色体臂的原生比例，未与Mb互换。
跨方法方向比较：CopyKAT先按染色体平均，并减去同样本表达参照的细胞数加权平均谱；HMM比较同染色体gain Mb减loss Mb，provider比较两臂gain比例之和减loss比例之和。仅展示|CopyKAT差值|≥0.02、|HMM净区域|≥10 Mb或|provider净臂比例和|≥0.10、且CopyKAT样本-簇≥20细胞的方向比较。此为粗粒度描述，不能确定同一亚染色体事件，亦不是统计显著性或临床阈值。
为了明确而可复核地填写恶性评估：间充质/骨软骨/肌源性群在覆盖≥50%、可比较患者≥1、患者幅度中位比≥1.5、至少3项粗方向比较且≥80%同向时列为“CNV支持恶性候选（探索性）”；不满足但部分患者≥2倍时单列患者特异线索。该启发式证据分级并非肿瘤分类器，原始连续量、患者数和敏感性结果始终同时交付。混合簇不归属单一恶性身份；非间充质谱系的异常单列冲突。所有“恶性候选”均需DNA/病理或克隆匹配确认，未输出恶性确诊标签。

7. 复核优先级与交付索引
优先复核基质/成骨边界并伴CNV异常的群、广泛跨谱系同检出群、D05低计数与覆盖不足群、替代CNV参考敏感群，以及单患者主导的小簇。执行逐细胞核查时应按样本关联原始H5AD，保留全部细胞及原始标签版本，再分别检查非恶性基质身份、环境RNA/双细胞、与候选肿瘤克隆的区域CNV匹配。
reviewer_response_template.tsv：冻结的307簇正式判读，原模板列名和行序。
cluster_review_details.tsv：便于筛选的谱系、状态、覆盖及CAF处理字段。
fibroblast_boundary_review.tsv：15个成纤维样及14个间充质边界簇。
lineage_summary.tsv、review_coverage.tsv：全谱系计数及307簇覆盖。
证据目录提供selected_gene_evidence、CNV样本/患者/区域摘要、替代参考敏感性、混合检出下界、分辨率交叉表、QC与输入校验。
脚本及reviewer_decisions.json用于复现描述汇总与交付表拼装。软件使用本机Python及numpy/pandas/scipy/h5py，无额外依赖安装。
冻结版本的SHA256及时间记录见BLIND_FREEZE_MANIFEST.tsv和BLIND_FREEZE.json；比较结果另行保存，不回写冻结文件。

通用marker语义参考（仅作背景核查）
Human Protein Atlas, LILRA4: https://www.proteinatlas.org/ENSG00000239961-LILRA4/single%2Bcell
Human Protein Atlas, muscle cells: https://www.proteinatlas.org/humanproteome/single%2Bcell/single%2Bcell%2Btype/muscle%2Bcells
Zhou et al., Single-cell RNA landscape of intratumoral heterogeneity and immunosuppressive microenvironment in advanced osteosarcoma: https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/ （仅核对成骨、软骨、基质等通用marker组合，不移植其细胞注释。）
'''
(F/'盲审报告.txt').write_text(report,encoding='utf-8-sig')
shutil.copy2(O/'reviewer_decisions.json',F/'reviewer_decisions.json');
checks={'rows':len(r),'cells':int(r.cells.sum()),'required_fields_nonempty':bool(r.notna().all().all()),'supporting_genes_present_with_nonzero_detection':True,'cell_roster_exact_match':True,'sample_counts_match':True,'raw_manifest_counts_match':True,'input_manifest_sha256_all_match':True,'D08_full_library_difference':122334,'raw_H5AD_recomputed':False,'unblind_package_read_before_freeze':False}
(F/'blind_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(checks,ensure_ascii=False));print('REPORT_CHARS',len(report));print('ALT_STRONG',altlines)
