import pathlib,pandas as pd,json,hashlib,datetime,shutil,numpy as np
P=pathlib.Path(__file__).parent;O=P/'review_output';C=O/'unblinded_comparison';B=O/'blind_submission';U=P/'unblinding_materials';m=pd.read_csv(C/'cluster_annotation_comparison.tsv',sep='\t');sm=pd.read_csv(C/'comparison_summary.tsv',sep='\t');bd=pd.read_csv(B/'cluster_review_details.tsv',sep='\t');fr=json.loads((B/'BLIND_FREEZE.json').read_text());mix=pd.read_csv(B/'evidence/mixed_program_evidence.tsv',sep='\t');
# Correct a narrative-only sample count while preserving every frozen byte.
original=(B/'盲审报告.txt').read_text(encoding='utf-8-sig');assert '覆盖83/107样本' in original;corrected='勘误版：在揭盲后逐细胞映射复核中，将CopyKAT覆盖样本数由83更正为84；其他判读文字、307簇标签与数值不变。冻结原版及其SHA256继续保留于blind_submission。\n\n'+original.replace('覆盖83/107样本','覆盖84/107样本');(C/'盲审报告_仅样本计数勘误.txt').write_text(corrected,encoding='utf-8-sig')
errata={'field':'CopyKAT连续谱覆盖样本数（盲审报告正文）','frozen_report_value':83,'correct_value':84,'reason':'83为存在所选表达参照细胞的样本数；84为有效CopyKAT连续谱实际样本数。79份达到≥30参照细胞、≥2参照簇的摘要条件。','unchanged':'568754细胞覆盖、逐来源覆盖、307簇全部判读、样本/患者级CNV比较与任何原始输入均未改变。','verified_from':'盲审包copykat_cell_continuous_metrics.tsv.gz的sample_id去重计数，以及冻结evidence/copykat_sample_cluster_comparison.tsv。','timing':'冻结后复核发现，公开勘误；未回写冻结原版'}
(C/'盲审报告勘误.json').write_text(json.dumps(errata,ensure_ascii=False,indent=2),encoding='utf-8')
lines=[f'{r.comparison_category}：{r.clusters}簇（{r.cluster_fraction:.1%}），{r.cells:,}细胞（{r.cell_fraction:.1%}）。' for r in sm.itertuples()]
fba=m[(m.lineage_code=='Fibro')&(m.proposed_lineage!='Fibroblast_like')]
fbtxt='\n'.join(f'{r.cluster_id} = {r.dataset} 原K{r.original_cluster}：{r.cells:,}细胞，谱系把握{r.lineage_confidence}，原标签{r.proposed_lineage}。' for r in fba.itertuples())
mixtxt=[]
for r in m[m.comparison_category=='新增混合表达风险'].itertuples():
 g=mix[mix.cluster_id==r.cluster_id].sort_values('minimum_codetection_fraction',ascending=False);z=g.iloc[0] if len(g) else None;bound=f'{z.gene_A}/{z.gene_B}同检出下界{z.minimum_codetection_fraction:.1%}' if z is not None and z.minimum_codetection_fraction>0 else '边际检出率不能确定同检出比例，需原始细胞复核'
 mixtxt.append(f'{r.cluster_id} = {r.dataset} 原K{r.original_cluster}（{r.cells}细胞）：原{r.proposed_lineage}，盲审{r.reviewed_lineage}；{bound}。')
cn=pd.read_csv(C/'proposed_CNV_cluster_composition.tsv',sep='\t');cnt=cn.groupby('CNV_evidence_current').cells.sum();cntext='\n'.join(f'{k}: {v:,}细胞（{v/876648:.1%}）' for k,v in cnt.items())
report=f'''CAF重注释：冻结盲审与新注释揭盲比对
日期：2026-10-03

1. 主要表达谱系总体相容；候选细化和混合表达构成主要差异
对307个主簇、876,648个细胞完成一一比对。盲审已于{fr['frozen_at_UTC']}冻结，此后才读取新注释包。冻结判读表SHA256为{fr['response_SHA256']}。比对未覆盖盲审表，也未改写新注释。
'''+ '\n'.join(lines)+'''
这里“主要谱系方向一致”是将T/Treg/NK归到T_NK、单核/巨噬/DC归到Myeloid、血管/淋巴内皮归到Endothelial等上位范围后的相容性，不是原字符串相等。41簇“未定范围一致”不算身份已解决；其余差异也不等于任何一方已被证明错误。骨软骨亚型范围差异独立保留，没有挤入完全一致类别。没有已验证金标准，因此不报告准确率、灵敏度、特异度或以相容率宣布全批注释通过。
219簇明确方向一致与41簇未定范围一致合计260/307簇（84.7%），覆盖654,689细胞（74.7%）；这一合计必须同时保留“未定”性质。未发现需要将一个明确主谱系整体改成完全无关谱系的互斥命名，但6个原单谱系簇具有新增混合风险，不能因此宣称全部单簇纯度成立。

2. 原11个成纤维样簇全部得到表达谱系支持，额外4簇仍属候选扩展
新包Fibroblast_like为11簇、45,721细胞，全部落在冻结盲审的Fibro类别。盲审另有4簇、45,291细胞提出成纤维样方向，合计15簇、91,012细胞。新增4簇原来均为Mesenchymal_unresolved：
'''+fbtxt+'''
这4簇的新增量较大，不能仅因方向被细化便自动扩充正式CAF集合。D05_K018的CNV覆盖很低，D07_K041的计数和覆盖有限；D02_K012为患者集中群。保留原未定标签、盲审方向、置信描述及连续证据并行记录更合适。
两方一致的Fibroblast_like也不是已确认非恶性CAF。D01_K017（GSE162454原K0）有明确成纤维样表达，同时存在多患者CNV异常，应在CAF正式使用前复核恶性/克隆身份。原新注释包本身亦明示其Fibroblast_like只是候选，因此这里是强化核查条件，并非指出原包已经将其病理确诊为CAF。
另14个盲审间充质边界簇保留不确定性，不强行并入上述15个成纤维样集合。

3. 六个原单谱系簇建议增加混合风险字段
'''+ '\n'.join(mixtxt)+'''
上述是转录信号层面的混合审查。检出率下界max(0,pA+pB-1)只能推导同一条形码的共检出，不能单独确诊双细胞，也不证明发生转分化。建议将mixed_flag作为后续核查项，在完成原始条形码共表达、环境RNA/双细胞和样本内结构复核前，避免将全簇当作纯细胞群开展谱系特异的推断。
另外6簇与原未定范围相容，但盲审进一步指出混合成分；这类不计为原注释主要谱系错误。D07_K039的B/浆细胞混合程序与原B_Plasma宽标签相容。

4. 30个原未定簇获得候选细化，4簇保留骨软骨亚型取舍
30簇、140,177细胞的主要变化来自完整基因面板对状态主导簇的补充。例如D01_K012、D02_K004、D06_K026、D07_K019、D08_K012/K021及D09_K004在增殖之外保留成骨程序，冻结盲审将其标为成骨样。D01_K028保留IFN状态同时给出成骨样方向。这些细化不依赖原标签回填。
低计数或线粒体主导者如D02_K002、D08_K004、D07_K007仅保留中等把握；D04_K019/K028和D05_K030由“全谱系未定”细化为“间充质边界未定”，不能计作完全解决。
4个原Osteochondrogenic_like簇在盲审中倾向Osteo：D02_K015、D05_K001、D05_K007、D09_K009，共77,269细胞。双方骨软骨间充质上位方向相容；成骨和软骨具体程序的强弱、覆盖和技术背景需继续区分。盲审并未把它们改成成纤维细胞。
D03_K026仅11个细胞，原Immune_unresolved；冻结盲审保留全谱系未定，比原包更保守。没有依据增加这类小簇的推断性统计。

5. CNV身份结论须与表达标签独立比较
原新注释逐细胞CNV类别复算如下：
'''+cntext+'''
原包中的method_conflict占45.0%，不能被本次粗粒度区域方向一致性摘要掩盖。原包是逐细胞原生分类call组合；盲审使用同样本表达参照下的连续幅度及染色体级区域方向，不同分辨率、阈值和覆盖下本就不能将两者直接换算。比对表保留每簇原CNV类别比例与冻结盲审连续证据，供定位冲突。
原包明确将concordant_alteration写为RNA CNV支持、非病理确诊；单方法或冲突均保留身份未定。冻结盲审同样未给病理恶性确证，43簇的“CNV支持恶性候选”是报告中说明的描述性启发分级。其余基因组与表达冲突、参考敏感性及缺失均须保留，不能按某一方的标签数量投票决定真值。
S027、S091的两套替代参考在强区域摘要中仍有方向变化，详见冻结证据表。参考预测也不是独立金标准：原包说明ScPCA曾使用较早版BlueprintEncode；HPCA/Blueprint的生物类别覆盖与命名范围不同。比对表保留原始及原生剪枝标签，不把两个参考预测当作两次独立病理验证。

6. 旧F3/F4对应关系进一步显示候选集合不宜整体沿用
对揭盲包已有历史集合交叉表作精确细胞映射，不读取外部历史项目文件。
supported F3的343细胞：盲审Fibro 214、Osteo 93、Mesench 19、U 14、Chondro 2、Macro 1。214个成纤维样细胞与新包一致，主要为D07_K010中的213个细胞，另1个来自D02_K008。成纤维表达支持并不解除该群CNV覆盖/身份限制。
expanded F4的5,304细胞：盲审Osteo 3,297、Fibro 909、U 593、Mesench 412、其余93。909个Fibro中，164个与新包原Fibroblast_like重叠，745个来自原Mesenchymal_unresolved的候选细化。这745个不是额外验证通过的CAF。两个旧集合均不应整体作为已确认CAF继续解释；完整其他历史类别交叉也已保留。

7. 连接与版本完整性核查通过，另有一处公开勘误
新注释包267个manifest文件的大小及SHA256全部通过；CURATION_LOCKED.json实测SHA256与交付说明一致。主簇键按source_id/cluster_id连接dataset/original_cluster，307簇细胞数与患者数完全相同。876,648细胞通过dataset+原cell_id一一连接匿名cell_id，主簇归属与盲审roster逐细胞相同；逐细胞full_library汇总也与盲审clusters.tsv完全相同。不存在靠行序猜测身份、丢失簇或按标签相似度配对。
盲审报告中“CopyKAT覆盖83/107样本”应为84/107；83为存在所选表达参照细胞的样本数，79份达到参照摘要条件。568,754细胞和逐簇/患者结果一直正确。冻结原版保持不动，单独提供盲审报告勘误.json及仅修正该处样本数的报告副本。这一勘误发生在揭盲后复核，已明确记录，不伪装为冻结前内容。
D08完整文库/公共feature差额仍按冻结盲审保留；新包方法说明也确认使用公共feature。未取得原始H5AD，不声称独立验证逐基因差额。

8. 本次交付的使用顺序
先查看本报告及cluster_annotation_comparison.tsv。fibroblast_candidate_comparison.tsv集中展示15个候选；differences_and_unresolved.tsv保留其余差异和未定情况。cell_annotation_comparison.tsv.gz提供全部876,648匿名细胞的双版本标签连接，不替换原注释。历史集合文件仅为溯源比较。
优先处理6个新增混合风险簇、D01_K017等伴CNV异常的成纤维/间充质群，再讨论4个成纤维样候选扩展及状态主导簇的细分。当前比对支持保留大部分主要表达谱系，并将明确复核事项带入下一版；不支持宣布整批CAF身份已获验证。
本任务止于独立判读、比较与复核建议，未改写原分析配置、未替换正式逐细胞标签，未重算下游正式统计。
'''
(C/'揭盲比对报告.txt').write_text(report,encoding='utf-8-sig')
# Validate newly written comparison, cell mapping, and unchanged frozen submission.
assert len(m)==307 and m.cluster_id.is_unique and m.cells.sum()==876648
assert not m[['dataset','cluster_id','curation_id','proposed_lineage','reviewed_lineage','comparison_category','recommended_action']].isna().any().any()
assert sm.clusters.sum()==307 and sm.cells.sum()==876648
fb=m[m.proposed_lineage=='Fibroblast_like'];assert len(fb)==11 and fb.cells.sum()==45721 and (fb.lineage_code=='Fibro').all()
assert len(fba)==4 and fba.cells.sum()==45291
cc=pd.read_csv(C/'cell_annotation_comparison.tsv.gz',sep='\t',usecols=['cell_id','cluster_id','proposed_lineage','lineage_code']);assert len(cc)==876648 and cc.cell_id.is_unique
assert cc.groupby('cluster_id').size().sort_index().equals(m.set_index('cluster_id').cells.sort_index())
assert (cc.proposed_lineage==cc.cluster_id.map(m.set_index('cluster_id').proposed_lineage)).all();assert (cc.lineage_code==cc.cluster_id.map(m.set_index('cluster_id').lineage_code)).all()
actualcovered=pd.read_csv(B/'evidence/copykat_sample_cluster_comparison.tsv',sep='\t').sample_id.nunique();assert actualcovered==84
fmanifest=pd.read_csv(B/'BLIND_FREEZE_MANIFEST.tsv',sep='\t');fv=[]
for x in fmanifest.itertuples(index=False):fv.append({'path':x.path,'bytes_match':(B/x.path).stat().st_size==x.bytes,'sha256_match':hashlib.file_digest(open(B/x.path,'rb'),'sha256').hexdigest()==x.SHA256})
assert all(x['bytes_match'] and x['sha256_match'] for x in fv);pd.DataFrame(fv).to_csv(C/'frozen_blind_integrity_after_comparison.tsv',sep='\t',index=False)
assert hashlib.file_digest(open(B/'BLIND_FREEZE_MANIFEST.tsv','rb'),'sha256').hexdigest()==fr['freeze_manifest_SHA256']
validation={'status':'complete','completed_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'clusters':307,'cells':876648,'all_cluster_links_one_to_one':True,'all_cell_links_one_to_one':True,'cell_and_cluster_labels_consistent':True,'all_37_frozen_files_unchanged':True,'all_input_checksums_valid':True,'new_package_curation_hash_matches':True,'CNV_actual_samples':84,'narrative_only_erratum_recorded':True,'original_annotation_files_unchanged':True,'no_downstream_statistics_recomputed':True};(C/'COMPARISON_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8')
S=C/'scripts';S.mkdir(exist_ok=True)
for name in ['prepare_comparison.py','compare_cells.py','classify_comparison.py','finalize_comparison.py']:shutil.copy2(P/name,S/name)
(S/'复现说明.txt').write_text('在工作解压根目录运行，需blind_submission冻结原版及unblinding_materials新注释包。依次prepare_comparison.py、compare_cells.py、classify_comparison.py、finalize_comparison.py。脚本使用原精确cell/cluster键，比较映射不更改冻结标签。不要对冻结目录作任何写入。',encoding='utf-8-sig')
print(json.dumps(validation,ensure_ascii=False));print('REPORT_CHARS',len(report))
