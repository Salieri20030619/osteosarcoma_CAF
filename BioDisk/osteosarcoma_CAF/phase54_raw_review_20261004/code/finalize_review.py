from pathlib import Path
import json,hashlib,datetime,zipfile,os
import numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase54_raw_review_20261004';S=P/'phase53_de_novo_reannotation_20261002';V=P/'review_intake_20261004/received_review'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
summary=read(O/'tables/target_cluster_raw_review_summary.tsv').set_index('cluster_id');cnv=read(O/'tables/CNV_target_summary.tsv').set_index('cluster_id');registry=read(O/'tables/review_target_registry.tsv');pairs=read(O/'tables/cluster_pair_codetection.tsv');gene=read(O/'tables/cluster_gene_patient_summary.tsv');dbl=read(O/'tables/doublet_sample_cluster_evidence.tsv');pooled=read(O/'tables/CNV_pooled_osteogenic_similarity.tsv');details=read(V/'blind_submission/cluster_review_details.tsv').set_index('cluster_id')
primary={'D01_R10_K002':('PF4','C1QA'),'D01_R10_K034':('IFITM5','C1QA'),'D04_R10_K007':('CD163','DCN'),'D06_R10_K022':('RUNX2','CD163'),'D07_R10_K027':('PECAM1','SFTPB'),'D09_R10_K010':('C1QA','COL3A1')}
expansion=['D02_R10_K012','D05_R10_K018','D07_R10_K024','D07_R10_K041']
notes={
'D01_R10_K002':'稀疏血小板/髓系共信号；原始PF4/C1QA仅7/91条形码共检出，校正量均>=1仅1/91。70/91缺少有效双细胞分类。保留稀疏混合风险，当前数据不支持将全簇判成双细胞。',
'D01_R10_K034':'原始成骨/巨噬共信号成立；髓系成分对DecontX敏感，IFITM5/C1QA校正量均>=1为0/491。双细胞模型亦有支持。保留混合风险，不能按校正结果自动认定为纯成骨群。',
'D04_R10_K007':'原始髓系/基质共表达及双细胞风险成立；校正后髓系计数大幅降低而DCN/COL3A1主要保留。原Myeloid标签不足以支持全簇纯髓系解释。',
'D06_R10_K022':'成骨/髓系混合风险成立，95.4%细胞来自P057；RUNX2计数主要保留，CD163明显衰减。样本重复及纯谱系解释均受限。',
'D07_R10_K027':'内皮/肺上皮混合风险强；95/100原始共检出PECAM1/SFTPB，校正量均>=1仍49/100，provider标为双细胞99/100。不得作为纯内皮群用于谱系特异解释。',
'D09_R10_K010':'髓系/基质混合风险强；原始C1QA/COL3A1共检出1314/1533，校正量均>=1仍562/1533；已评估细胞1204/1450被模型标为双细胞。未被模型标记者仍有共表达，不能靠单次预测清除全部身份问题。',
'D01_R10_K017':'成纤维样表达经校正仍支持，但多患者CNV异常及与骨源样群的谱形相似性阻止将全簇解释为已确认非恶性CAF；不同患者并非同质。',
'D02_R10_K012':'可保留成纤维样候选方向；DCN/LUM/COL3A1校正后仍有表达，79.3%细胞集中于P039。5名可比较患者CNV幅度中位比0.71，但21125/22662原生CNV方法冲突，不能由低幅度认定正常。',
'D05_R10_K018':'跨患者基质方向得到支持，细胞层面覆盖稀疏；原始文库中位339.5计数，CNV仅覆盖441/6314。保留低计数成纤维样/间充质候选，不扩为确认CAF。',
'D07_R10_K024':'多患者基质方向稳定，可保留成纤维样候选；校正后COL3A1/DCN/LUM均17/19患者相对同样本其余细胞正向，同时RUNX2广泛表达且存在患者特异CNV异常。',
'D07_R10_K041':'基质候选方向较弱；校正后LUM仅7/20患者正向，DCN为13/20，COL3A1为15/20，CNV覆盖47.2%。继续保留间充质未定及成纤维样候选方向并列。',
'D08_R10_K017':'间充质边界仍未解决；上游双细胞预测在已评估者中阳性94.6%，98.3%细胞由单患者贡献。CNV异常不能直接归给纯成纤维细胞或独立肿瘤克隆。'}
over=[]
for r in registry.itertuples(index=False):
 cl=r.cluster_id;q=summary.loc[cl];c=cnv.loc[cl];label=details.loc[cl,'lineage_code']
 if cl in primary:decision='sparse_mixed_signal' if cl=='D01_R10_K002' else 'mixed_expression_risk_supported'
 elif cl=='D01_R10_K017':decision='fibro_expression_supported_CNV_identity_unresolved'
 elif cl in expansion:decision='candidate_direction_retained_limited' if cl in ['D05_R10_K018','D07_R10_K041'] else 'candidate_direction_supported_identity_unresolved'
 elif label=='Fibro':decision='fibro_expression_retained_identity_unresolved'
 elif label=='Mesench':decision='mesenchymal_boundary_unresolved'
 else:decision='mixed_program_unresolved'
 note=notes.get(cl,'保留成纤维样表达候选；本轮已核对原始/校正后表达、患者分布、双细胞及CNV证据，不据此作非恶性CAF确证。' if label=='Fibro' else '原始/校正后表达、患者分布、双细胞和CNV均已汇总；证据仍不足以作确定谱系或恶性身份定名，维持未定/混合状态。')
 over.append({'source_id':r.source_id,'dataset':r.dataset,'cluster_id':cl,'original_cluster':r.original_cluster,'cells':int(q.cells),'patients_observed':int(q.patients),'original_proposed_lineage':r.proposed_lineage,'blind_lineage':r.reviewed_lineage,'phase54_review_decision':decision,'interpretation_ZH':note,'max_patient_fraction':q.max_patient_fraction,'doublet_assessed_any':int(q.doublet_assessed_any),'doublet_positive_any':int(q.doublet_positive_any),'copykat_coverage':c.copykat_coverage,'CNV_comparable_patients':int(c.comparable_patients),'CNV_median_amplitude_ratio':c.median_fold,'formal_CAF_identity_confirmed':False,'formal_annotation_replaced':False})
over=pd.DataFrame(over);assert len(over)==42 and over.cells.sum()==127395;over.to_csv(O/'tables/CLUSTER_REVIEW_DECISIONS.tsv',sep='\t',index=False,encoding='utf-8-sig')
co=read(O/'tables/target_cell_review_overlay.tsv.gz');co=co.merge(over[['cluster_id','phase54_review_decision']],left_on='r10',right_on='cluster_id',validate='many_to_one');assert len(co)==127395 and co.cell_id.is_unique;co.to_csv(O/'tables/cell_review_decisions.tsv.gz',sep='\t',index=False)
# Corrected counts are fractional. Explicitly show >0 sensitivity rather than implying >=1 is an unchanged detection rule.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
selected=[]
for cl,(a,b) in primary.items():selected.append(pairs[(pairs.cluster_id==cl)&(pairs.gene_A==a)&(pairs.gene_B==b)].iloc[0])
q=pd.DataFrame(selected);fig,ax=plt.subplots(figsize=(11,6));yy=np.arange(len(q));ax.barh(yy-.23,q.raw_both_fraction*100,.21,label='Raw: both > 0 (integer counts)',color='#286c8e');ax.barh(yy,q.corrected_both_gt0_fraction*100,.21,label='DecontX: both > 0 (fractional counts)',color='#8bb8cc');ax.barh(yy+.23,q.corrected_both_ge1_fraction*100,.21,label='DecontX: both >= 1 (sensitivity)',color='#d28b37');ax.set_yticks(yy,[f'{r.cluster_id}\n{r.gene_A} + {r.gene_B} ({r.cells} barcodes)' for r in q.itertuples()]);ax.set_xlim(0,102);ax.invert_yaxis();ax.set_xlabel('Percentage of barcodes within cluster');ax.set_title('Raw joint transcription and correction sensitivity');ax.legend(loc='upper center',bbox_to_anchor=(.45,-.13),ncol=1,frameon=False);fig.tight_layout(rect=[0,.10,1,1]);fig.savefig(O/'figures/mixed_codetection.png',dpi=180);fig.savefig(O/'figures/mixed_codetection.pdf');plt.close(fig)
# Complete internal review text with findings, methods, and centralized interpretation limits.
lines=['CAF盲审后原始数据复核报告','阶段：phase54_raw_review_20261004；日期：2026-10-04','本轮为揭盲后的数据复核，不是新的独立盲审。','',
'一、原始数据复核确认混合风险需要分层处理，CAF身份不能由表达标签直接确证',
'已核验107份匿名原始H5AD的SHA256与全部876,648细胞，逐样本核对既有DecontX校正矩阵的细胞/基因顺序、原始文库和校正文库。对15个成纤维样、14个间充质边界及13个混合程序簇共42簇、127,395细胞形成原始与校正后证据和附加判读层。84份CopyKAT细胞×分箱谱覆盖568,754细胞，样本幅度比与盲审结果一致，染色体均值也通过直接重算核对。',
'5个重点簇保留有直接证据的混合表达风险；D01_K002只有稀疏血小板/髓系共信号，不支持把全簇判成双细胞。D01_K017的成纤维样表达和CNV身份冲突同时成立。新增4个簇均保留候选层次，不能整体并入已确认CAF。D08公共feature计数差额已完全解释，关闭这一数据完整性疑点。','',
'二、六个重点簇的原始共检出得到直接测量，校正敏感性与双细胞预测分开呈现',
'下列“校正>=1”表示两个feature的DecontX校正量均至少1，是描述性敏感性口径；它不是细胞过滤或身份判定阈值。完整表同时保留校正>0和>=2。']
for cl,(a,b) in primary.items():
 row=pairs[(pairs.cluster_id==cl)&(pairs.gene_A==a)&(pairs.gene_B==b)].iloc[0];v=summary.loc[cl]
 lines += [f'{cl}（{int(row.cells)}细胞；{int(v.patients)}名来源内患者记录）：{a}/{b}原始共检出{int(row.raw_both)}/{int(row.cells)}（{row.raw_both_fraction:.1%}）；校正>=1为{int(row.corrected_both_ge1)}/{int(row.cells)}（{row.corrected_both_ge1_fraction:.1%}）。现有方法至少一种标为双细胞{int(v.doublet_positive_any)}/{int(v.doublet_assessed_any)}已评估条形码；另{int(v.cells-v.doublet_assessed_any)}未评估。',notes[cl]]
lines += ['DecontX的>0共检出比例与原始值相同，说明大量信号变为非零小数而非被精确置零。不能将>=1口径下的下降写成共表达已消失，也不能据此独立证明环境RNA来源。',
'D04_K007中CD163计数保留12.0%，DCN保留95.1%；D06_K022中CD163保留12.1%，RUNX2保留99.7%；D01_K034中C1QA保留21.4%，IFITM5保留96.5%。这显示髓系成分对模型校正更敏感。D07_K027和D09_K010校正后仍有明显双程序共信号。',
'D09_K010中即使限定到“已评估且未被任何方法标阳性”的246个条形码，C1QA/COL3A1原始共检出仍82.1%、校正>=1仍41.5%。现有双细胞标签不能单独解决全部混合身份。',
'双细胞预测逐方法表保留provider scDblFinder、既有上游预测和phase31补算结果；“任一阳性”只用于描述联合风险，不是独立重复投票。D04_K007的phase31为532/620，上游为469/812；D06_K022的phase31为177/207，上游为121/194；D01_K034的phase31为30/31，上游为73/295。不同方法覆盖数不得互换。','',
'三、D01_K017的成纤维样表达稳健，但非恶性身份不成立为已确认事实',
'2,720细胞中，DecontX后DCN、LUM、COL3A1总量分别保留95.1%、95.9%、92.4%；相对同样本其余细胞的患者级方向三者均6/6正向。RUNX2同样6/6正向。成纤维样基质表达不是仅由可校正背景信号解释。',
'CopyKAT覆盖2,666/2,720，6名可比较患者的幅度中位比2.93，4/6超过2倍。P032占全簇41.0%，其幅度比却为0.805；P066为1.729，其余4名为2.887–3.891。不得把少数高幅度患者结论直接推广到整个簇。',
'不挑选最相似的单簇，使用同样本全部盲审Osteo群按细胞数加权的染色体谱作为比较对象；D01_K017去参考信号的描述性Pearson r按患者为：']
for r in pooled[pooled.cluster_id.eq('D01_R10_K017')].sort_values('patient_id').itertuples(index=False):lines.append(f'{r.patient_id}：r={r.pooled_osteo_chromosome_r:.3f}，目标CopyKAT细胞{r.target_copykat_cells}、Osteo比较细胞{int(r.pooled_osteo_cells)}。')
lines += ['5名患者r为0.798–0.986，P066为-0.058；保留这一矛盾，不能宣布6名患者共享同一克隆。此处相关系数是同一RNA推断方法的描述性谱形相似性，未作染色体/bin级推断检验。',
'原生逐细胞类别中1,332细胞为双方法alteration、1,329为方法冲突，仅5为双方法no alteration。簇标签保留成纤维样表达；非恶性CAF身份与克隆身份仍未解决。其他成纤维样及14个间充质边界簇的对应数据全部保留于CNV_target_summary和逐患者表。','',
'四、新增四个候选应保留不同证据强度，不整体扩充正式CAF集合']
for cl in expansion:
 q=summary.loc[cl];c=cnv.loc[cl];lines += [f'{cl}：{int(q.cells)}细胞、{int(q.patients)}名来源内患者记录；最大患者贡献{q.max_patient_fraction:.1%}；CNV覆盖{c.copykat_coverage:.1%}。',notes[cl]]
 for g in ['COL3A1','DCN','LUM','PDGFRA','RUNX2']:
  a=gene[(gene.cluster_id==cl)&(gene.gene==g)].iloc[0];lines.append(f'  {g}：原始/校正后患者平均log2效应{a.raw_mean_effect:.3f}/{a.corrected_mean_effect:.3f}；校正后正向{int(a.corrected_positive_patients)}/{int(a.patients)}。')
lines += ['D05_K018在r2.0下两大子群分别占45.8%和45.3%，后者保留更多骨源样程序；D07_K024两主要子群占77.6%和10.4%；D02_K012两主要子群占47.1%和51.2%；D07_K041的97.7%落在单一r2.0子群。分群稳定或不稳定均不等于谱系已获确证。',
'14个间充质边界簇全部保留未定。其中D08_K017双细胞阳性比例94.6%、最大患者贡献98.3%，不能把其CNV异常解释为纯成纤维群的独立恶性证据；D06_K018的4名可比较患者CNV幅度中位比2.28，D09_K002只有1名满足比较条件且中位比2.38。完整42簇附加判读见CLUSTER_REVIEW_DECISIONS.tsv。','',
'五、D08计数差额完全来自公共feature之外，原交付簇计数与原始数据一致',
'逐样本对完整H5AD计算全feature文库、32,732公共feature文库及非公共feature计数。非公共feature的非零计数之和为122,334，与盲审发现的差额完全相等。对31个簇×32,732个公共feature逐元素重建，全部与原交付primary_cluster_counts矩阵一致。逐样本差额和每个非公共feature的计数已导出。该疑点关闭，无需补造、删除或重算原输入。','',
'六、方法',
'输入为phase53匿名原始H5AD、细胞/簇对应键、冻结盲审及揭盲比较；原输入共107份、9来源、876,648细胞。原始SHA256、细胞唯一性、原样本细胞顺序、gene/feature顺序及原始/校正文库逐样本核对。复用phase31退出码0且具有完成标记的DecontX结果，使用原生自动广谱聚类设置；同时读取已有provider/项目双细胞字段及其实际评估状态。',
'重点42簇包含15 Fibro、14 Mesench及13混合程序簇；不新增细胞筛选。预先列出八类多基因面板，COL1A1、LOX、MMP14不参与本轮谱系面板。面板feature在来源共同测量集合内按feature ID确定；原盲审成对基因另按冻结feature证据和边际检出率精确定位。不存在把重复符号检出率相加的处理。',
'每样本/簇/基因报告原始计数、校正量、检出率、CP10k，以及相对同样本其余细胞的log2((目标CP10k+0.1)/(其余CP10k+0.1))。同患者多个技术样本的效应先等权平均，再在来源内报告患者方向及效应；未将不同来源的匿名患者视为独立重复。保留全部样本和患者，无额外最低细胞数用于表达方向筛选。',
'全部42簇的r2.0子簇表达摘要及六个重点簇按患者、r2.0、双细胞证据分层的共检出均另表输出。面板均值仅作为描述性结构证据，不按最大分数直接定名。',
'CNV原生幅度复算沿用盲审口径：表达定义的T/Treg/NK/B/内皮/淋巴内皮/肥大群作同样本参考，参考至少30个有效细胞且来自至少2簇，目标至少20个有效细胞；这些是摘要可比较条件，不删除任何输入细胞。按样本原生mean_abs信号中位数之比，再对患者多样本取中位数。',
'84份细胞×分箱HDF5直接重建染色体均值；区域方向比较复用冻结原生区域汇总，结合重新计算的CopyKAT信号，未重跑inferCNV。原逐细胞call冲突单独保留。候选与骨源样群比较同时保留全部单簇比较以及不择优的同样本全部Osteo群加权比较。',
'本轮为描述性身份复核，未设置生物处理效应推断检验，未计算或补造P值/FDR。>=1/>0/>=2为计数敏感性展示，1.5或2倍为已有CNV幅度描述，不是FDR替代或细胞准入标准。','',
'七、讨论与剩余边界',
'本轮确认了真实条形码共检出、校正敏感性、患者分布和RNA-CNV冲突。DecontX没有空液滴背景，校正不会区分双细胞、环境RNA、细胞复合体或真实双程序；既有双细胞预测也不构成病理真值，未评估保持缺失。RNA-CNV与表达相关误差未被独立DNA、配对克隆或病理证据消除，因此不作恶性确诊，也不以“无异常”确认正常。',
'有效患者数与细胞数分开报告；不足10名来源内患者的结果保持探索性定位。即使覆盖较多患者，身份判读仍属于描述性证据。本轮未重新聚类全细胞、重新训练双细胞模型或重跑整套CNV；30个状态主导簇和4个骨软骨亚型取舍不在此次42簇原始深度复核范围内，原比较表继续保留，不宣布全307簇身份已通过。',
'阶段性处置：六个重点簇添加本轮风险结论；14个间充质边界及其余混合程序维持未定；15个成纤维样方向均不升级为已确认CAF。所有结论通过独立附加表保存，原始矩阵、phase53标签、冻结盲审和phase49正式统计包未替换。临床/功能外推没有新增确证。','',
'八、交付索引',
'CLUSTER_REVIEW_DECISIONS.tsv：42簇结论、证据计数和原标签。',
'cell_review_decisions.tsv.gz：127,395重点细胞的附加风险层，非正式注释替换。',
'cluster_pair_codetection.tsv / priority_mixed_stratified_review.tsv：全成对及六重点簇分层复核。',
'cluster_gene_patient_summary.tsv / patient_gene_direction.tsv.gz：表达方向与患者一致性。',
'CNV_target_summary.tsv / CNV_patient_recomputed.tsv / CNV_pooled_osteogenic_similarity.tsv：CNV幅度、冲突和同样本谱形比较。',
'doublet_sample_cluster_evidence.tsv / doublet_provenance.tsv：各方法独立分母及来源。',
'D08_count_reconciliation.tsv / D08_excluded_feature_counts.tsv：已关闭的计数差额。',
'cells/：全部107样本重点细胞的原始/校正后标志基因及面板值。',
'code/、resources/、logs/、QA和SHA256清单：复现与完整性凭据。']
report=O/'delivery/CAF盲审后原始数据复核报告_20261004.txt';report.write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
# Necessary output acceptance checks.
rawqa=json.loads((O/'RAW_REVIEW_QA.json').read_text());assert rawqa['samples']==107 and rawqa['cells']==876648 and rawqa['target_cells']==127395
assert set(over.cluster_id)==set(summary.index)==set(co.r10)
assert co.groupby('r10').size().reindex(over.cluster_id).tolist()==over.cells.tolist()
assert np.all((pairs.raw_both>=0)&(pairs.raw_both<=pairs.raw_A_detected)&(pairs.raw_both<=pairs.raw_B_detected))
assert np.all((pairs.corrected_both_ge2<=pairs.corrected_both_ge1)&(pairs.corrected_both_ge1<=pairs.corrected_both_gt0))
assert read(O/'tables/D08_count_reconciliation.tsv').excluded_feature_counts.sum()==122334
assert read(O/'tables/D08_excluded_feature_counts.tsv').counts.sum()==122334
stages={k:int((O/(k+'.exit_code')).read_text()) for k in ['raw_review','cnv_review','doublet_metadata','synthesis']};assert set(stages.values())=={0}
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':'raw_review_complete_identity_limits_retained','stages':stages,'input_cells':876648,'samples':107,'reviewed_clusters':42,'reviewed_cells':127395,'all_review_cell_ids_unique_and_cluster_counts_match':True,'raw_codetection_bounds_valid':True,'correction_sensitivity_counts_monotonic':True,'D08_missing_count_issue_closed':True,'formal_annotation_replaced':False,'downstream_statistics_recomputed':False,'P_values_generated':False,'full_307_cluster_identity_acceptance_claimed':False,'figure_visual_QA':'CNV plot inspected; mixed plot redrawn to show fractional >0 sensitivity and move legend outside plotting area; final PNG review pending.'}
dump(O/'FINAL_DATA_QA.json',qa)
(O/'PROGRESS.json').write_text(json.dumps({'state':'completed','completed_samples':107,'total_samples':107,'raw_exit_code':0}))
print(json.dumps({'report':str(report),'reviewed_clusters':len(over),'reviewed_cells':len(co),'QA_passed':True},ensure_ascii=True))