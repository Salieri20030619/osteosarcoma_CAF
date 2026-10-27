from pathlib import Path
import json,datetime,hashlib,zipfile,pandas as pd,numpy as np,shutil
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005'
assert json.loads((O/'QA/ACCEPTANCE.json').read_text())['passed']
impact=json.loads((O/'resources/IMPACT_SUMMARY.json').read_text())
w=pd.read_csv(O/'tables/program_before_after_comparison.tsv',sep='\t')
w['direction_flip_above_numeric_tolerance']=w.direction_flip&(w.mean_effect_all.abs()>1e-10)&(w.mean_effect_excluded.abs()>1e-10)
impact['program_direction_flips_above_numeric_tolerance']=int(w.direction_flip_above_numeric_tolerance.sum())
w.to_csv(O/'tables/program_before_after_comparison.tsv',sep='\t',index=False)
(O/'resources/IMPACT_SUMMARY.json').write_text(json.dumps(impact,indent=2))
m=pd.read_csv(O/'tables/desmosome_cell_counts.tsv.gz',sep='\t',low_memory=False,float_precision='round_trip')
raw=m[m.matrix=='raw'];cor=m[m.matrix=='corrected']
summary=pd.read_csv(O/'tables/subgroup_doublet_summary.tsv',sep='\t')
report="""骨肉瘤CAF双联体证据修复及SR1004标志物核查
2026-10-05｜phase66

一、结论

双联体证据传递问题已在当前后续版本中修复：8174个原始细胞均具备来源可追溯的检测状态，5728条漏带记录已补齐。完整上游模型阳性为26个，其中15个此前未进入下游汇总。原始F1–F5、33内部亚群、表达计数及已有置信度保留，另提供排除26个模型阳性的8148细胞敏感性视图。

主线结果没有改变，但辅助来源存在需要降级的结果。26个模型阳性来自GSE318840（15个）、GSE270231（10个）和GSE250015（1个）；按每亚群原始细胞数预先确定的主来源，没有程序效应改变。部分辅助来源的有效患者由2例降为1例，或由1例降为0例，已在新版统计中停止相应推断，不继续沿用旧P值。

SR1004可对应厂商Desmocollin-1抗体克隆，货号NBP3-22640，转录组应检索基因DSC1。它属于桥粒钙黏蛋白家族的Desmocollin支系；SR1004本身不是HGNC基因名。DSC1在本项目候选细胞中并非未覆盖，而是极低检出：8170个有覆盖的细胞中仅3个原始计数阳性，且各只有1个计数；不足以作为稳定CAF marker。

二、Methods

1. 元数据修复
以细胞ID从各样本原始full_called_counts.h5ad中的上游调用字段，以及已成功结束的phase31 scDblFinder结果恢复检测状态；同时保存每个调用者的原始类别、分数、可用阈值状态和文件来源。有效调用者中任一模型阳性即标记为“model_positive_requires_review”，不把模型预测当实验确证。
修复版保留legacy_doublet_assessed_any和legacy_doublet_positive_any；完整源调用独立重建后与修复值逐细胞比对。保守成纤维标记conservative_fibro按原规则重新计算，本次该派生集合实际未发生变化。
统一读取接口current_metadata.py按ID对齐，拒绝重复或未知ID，不把缺失填成False。分类及置信度由既有聚类/谱系证据产生，保留原值；双联体风险作为独立维度记录。

2. 敏感性分析
在固定33亚群上比较保留全部8174细胞与从目标、参考两侧同时排除26个模型阳性细胞。每个来源单独计算，重新检查患者/部位两侧的细胞数；采用原始/全部/至少1细胞、校正/全部/至少1细胞、校正/全部/至少5细胞、校正/Fibro参考/至少1细胞四种条件。
评估phase63的11个完整机制程序、3个胶原结构面板、phase65的3个主导基因对剔除面板，以及658个固定机制基因。以患者/部位伪bulk的全库CP10k加0.1取log2差，部位等权、患者等权；使用与phase63一致的检验函数（N≤15精确符号翻转，否则19999次模拟加观测）。N<2或全零不计算P；N<10按探索性呈现。各来源/矩阵/参考范围/阈值/程序类别/分析视图分别保留原定义的BH家族，基因与程序分开。原始P<0.05作为名义筛选标准，FDR同时报告。
未重跑无输入变化的bulk临床或空间固定RNA背景关联。此次范围为主线程序和固定机制基因敏感性，不宣称全部全转录组GSEA、轨迹或重新聚类结果均已复算。

3. SR1004与表达核查
通过原厂目录及NCBI/HGNC核实SR1004→Desmocollin-1→DSC1，HGNC:3035、Entrez1823、ENSG00000134765、UniProt Q08554。按来源的规范feature映射及样本feature-presence矩阵区分缺测与生物学零值，在完整基因矩阵中读取DSC1及DSC2/3、DSG1–4、DSP、JUP、PKP1–3作描述性背景。
原始计数≥1为主要检出描述；校正矩阵属于连续估计，同时保留>0与严格≥1指标及实际数值，不把接近1的浮点边界解释为表达消失。不为三个低计数细胞补造差异或特异性P值。

三、Results

3.1 5728条漏带记录完成恢复

根因位于phase55构建全量注释时，将phase54“部分目标细胞深度审阅表”左连接到全量细胞；不在该审阅子集的细胞得到缺失QC字段。phase59直接继承这些列，phase62验证的是字段保留一致，而非上游记录覆盖。新版本以完整上游调用重建证据，并增加有效覆盖、顺序对齐、未知/重复ID拦截检查。
8174/8174现有有效检测记录，26个模型阳性；主分析总体仍保留8174原始细胞，8148视图专用于敏感性，不直接把26个模型标记视为应永久删除的真双联体。
"""
a=pd.read_csv(O/'tables/cells_repaired_with_confidence.tsv.gz',sep='\t',low_memory=False)
for cls,t in a.groupby('class_id'):report+=f"{cls}：{len(t)}个细胞，模型阳性{int(t.doublet_positive_any.sum())}个。\n"
report+="\n3.2 F2-C3和F2-C1主来源结果保持不变\n\n"
for state,program,version in [('F2_C3','fibrillar_I_III_V','full'),('F2_C3','fibrillar_I_III_V','drop_dominant_pair'),('F2_C1','collagen_IV','full'),('F2_C1','collagen_crosslinking','full')]:
 r=w[(w.dataset=='SCPCP000017')&(w.subgroup_id==state)&(w.matrix=='corrected')&(w.scope=='all')&(w.minimum_cells==1)&(w.program==program)&(w.version==version)].iloc[0]
 report+=f"{state}，{program}，{version}：N={int(r.N_all)}→{int(r.N_excluded)}；效应{r.mean_effect_all:+.6f}→{r.mean_effect_excluded:+.6f}；P={r.P_excluded:.6f}；FDR={r.FDR_excluded:.6f}。\n"
report+="F2-C3的38个细胞中模型阳性为0。F2-C1模型阳性1/38来自GSE270231，不在上述7患者主来源中。因此主来源统计不变有明确的数据原因，不能据此声称所有来源都不受影响。小样本结果仍为探索性，名义显著不能写为FDR显著。\n\n3.3 辅助来源确有影响，已停止不具备重复的检验\n\n"
report+=f"共比较{impact['program_contrasts']}条程序配置结果，{impact['program_effect_changed']}条效应值发生变化，{impact['program_direction_flips_above_numeric_tolerance']}条在排除数值零后方向改变，{impact['program_lost_P']}条原有检验在排除后不再可评估。这些是多种来源/矩阵/参考条件下的统计行数，不是独立患者数。\n"
report+="其中GSE270231的F1-originalC4、F2-C1和F4-C1在部分主要条件下由2例变为1例；GSE318840的F2-C5也由2例变为1例。GSE250015的F2-C5、GSE270231的F3-C2、GSE318840的F4-C6等部分配置变为0例。对应P值设为缺失，不将失去检验条件写成生物学阴性。\n"
report+=f"程序名义P<0.05入选状态变化为{impact['program_nominal_changes']}条；{impact['gene_contrasts']}条固定机制基因配置的名义入选状态变化为{impact['gene_nominal_changes']}条。说明当前主线候选未被这项排除敏感性推翻，但辅助来源重复支持必须使用新版N和效应。全部受影响行、正负方向与FDR保存在结果表。\n\n3.4 SR1004应映射为DSC1，当前RNA检出极低\n\n"
report+="原厂将SR1004列为Desmocollin-1单克隆抗体克隆，货号NBP3-22640，列出的应用为WB。该克隆号不是表达矩阵中的基因名；不能自动替换成S100A4，也不能和Desmoglein支系的DSG1混用。厂商同一靶点其他克隆的IHC/Flow等验证不应转移到SR1004。\n"
report+="DSC1有8170/8174个细胞具备feature覆盖，另4个记作缺测而非零表达。3个原始阳性分别位于F4-C8、F4-C9和F4-C12，每个细胞只有1个原始计数；F1、F2、F3和F5中已覆盖细胞未检出。\n"
pos=cor[cor.DSC1>0]
report+=f"这3个细胞校正后仍保留约{pos.DSC1.min():.6f}–{pos.DSC1.max():.6f}的连续估计值，总量{pos.DSC1.sum():.6f}；三者均不是本轮模型阳性，排除26个细胞后保持同样的3个检出。严格≥1阈值与接近1的浮点值之间的差别不用于宣称基因表达消失。\n"
report+="因此，现在可以回答“此前未找到”的两个原因：SR1004是试剂克隆号，应按DSC1检索；按正确基因检索后，该基因在当前CAF候选池中也未形成可用的高表达或稳定亚群marker。3个F4细胞不能支持命名新的桥粒相关CAF亚群，更不能证明蛋白阳性或抗体染色特异性。\n\n四、Discussion与后续使用\n\n"
report+="模型阳性不等于已确认双联体，模型阴性也不保证单细胞纯净；现有分类和置信度作为历史可追踪证据保留，双联体风险独立标注。辅助来源样本量下降不会通过重复计算或改阈值补救，后续应采用本阶段更新结果。\n"
report+="DSC1属于胞间桥粒黏附相关分子，与细胞外胶原沉积是不同层面的机制。RNA极低检出不能否定蛋白存在；但当前数据不支持把它升级为CAF分型或胶原屏障主标志。若客户已有相关染色，应核对具体货号、抗体应用、组织部位和阳性区域，再决定是否作为独立蛋白线索。\n"
report+="原有两套真实病理应答的阴性结果、空间来源关联及独立高分辨率病理缺口保持原证据边界。本轮不新增耐药因果或蛋壳几何结论。\n\n五、当前入口与验收\n\n"
report+="修复主表：tables/cells_repaired_with_confidence.tsv.gz；完整双联体证据：tables/current_doublet_evidence.tsv.gz；逐调用者明细：tables/doublet_method_details.tsv.gz；变更日志：tables/repair_change_log.tsv.gz。\n"
report+="敏感性完整表：program_before_after_comparison.tsv、gene_before_after_comparison.tsv.gz、paired_roster.tsv.gz；DSC1及桥粒背景：desmosome_expression_summary.tsv、desmosome_feature_coverage.tsv、DSC1_trace_counts.tsv。\n"
report+="新版current_metadata.py已用于本阶段计算，并验证顺序打乱仍能按ID正确覆盖，未知/重复ID拒绝。独立重建逐调用者并集，重现phase63全部4466行主程序基线与P/FDR，核对每个患者/部位细胞数及受影响来源效应。原始分型与置信度逐列保持一致。历史阶段文件不覆盖，后续从当前入口读取修复版。\n\n来源\n"
report+="原厂产品：https://www.rndsystems.com/jp/products/desmocollin-1-antibody-sr1004_nbp3-22640\n原厂目录：https://www.novusbio.com/primary-antibodies/desmocollin-1\nNCBI DSC1：https://www.ncbi.nlm.nih.gov/gene/1823/\n"
name='CAF双联体证据修复与SR1004核查报告_20261005.txt';(O/'delivery'/name).write_text(report,encoding='utf-8')
pd.DataFrame([{'query':'SR1004','query_type':'antibody_clone','canonical_gene':'DSC1','catalog':'NBP3-22640','source_url':'https://www.rndsystems.com/jp/products/desmocollin-1-antibody-sr1004_nbp3-22640','note':'Project reagent-name lookup only; not an HGNC gene alias.'}]).to_csv(O/'resources/marker_name_registry.tsv',sep='\t',index=False)
# Root manifests are activated only after local delivery verification.
metadata_manifest={'version':O.name,'cells':8174,'assessed':8174,'model_positive':26,'sensitivity_cells':8148,'cell_table':str(O/'tables/cells_repaired_with_confidence.tsv.gz'),'overlay':str(O/'tables/current_doublet_evidence.tsv.gz'),'loader':str(O/'code/current_metadata.py'),'overlay_sha256':hashlib.sha256((O/'tables/current_doublet_evidence.tsv.gz').read_bytes()).hexdigest(),'cell_table_sha256':hashlib.sha256((O/'tables/cells_repaired_with_confidence.tsv.gz').read_bytes()).hexdigest(),'marker_name_registry':str(O/'resources/marker_name_registry.tsv'),'policy':'Preserve8174master;26model positives flagged and8148exclusion sensitivity separately; use strict loader for future8174 analyses.'}
(O/'resources/CURRENT_CELL_EVIDENCE_MANIFEST.json').write_text(json.dumps(metadata_manifest,indent=2))
(O/'QA/EXECUTION_STATUS.json').write_text(json.dumps({'repair':0,'sensitivity':0,'marker':0,'provenance':0,'verify':0,'exit_statuses_observed_from_native_sessions':True},indent=2))
files=[]
for directory in ['code','resources','tables','QA','logs']:
 for p in (O/directory).rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts and p.name!='DELIVERY_MANIFEST.json':files.append(p)
files += list(O.glob('*COMPLETE.json'))+[O/'delivery'/name]
manifest=[{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)]
mp=O/'resources/DELIVERY_MANIFEST.json';mp.write_text(json.dumps(manifest,indent=2))
zp=O/'delivery/CAF_doublet_repair_DSC1_20261005.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files+[mp]:z.write(p,arcname=O.name+'/'+str(p.relative_to(O)))
with zipfile.ZipFile(zp) as z:assert z.testzip() is None;nm=len(z.namelist())
receipt={'zip':str(zp),'bytes':zp.stat().st_size,'sha256':hashlib.sha256(zp.read_bytes()).hexdigest(),'files':nm,'report':name}
(O/'delivery/DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
print(json.dumps(receipt,ensure_ascii=False),flush=True);print('direction flips robust',impact['program_direction_flips_above_numeric_tolerance'])
