from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'geo_identity_20260919';T=O/'tables'
assert (O/'logs/identity.exit_code').read_text().strip()=='0'
assert (O/'TABLES_COMPLETE').exists()
S=pd.read_csv(T/'sample_identity_summary.tsv',sep='\t');A=pd.read_csv(T/'all_cell_identity_evidence.tsv.gz',sep='\t')
F=A[A.lineage_phase3=='FibLike'];state=json.loads((O/'STATUS.json').read_text())
assert len(S)==16 and len(A)==101487 and len(F)==8870 and A.cell_id.is_unique
J=pd.read_csv(T/'decontx_jobs.tsv',sep='\t')
assert (J.exit_code==0).all() and J.complete.all()
I=pd.read_csv(T/'integrity_checks.tsv',sep='\t')
assert I.cell_order_match.all() and I.raw_counts_conserved.all() and I.raw_panel_counts_reproduced.all()
# Verify identity review did not mutate any accepted input that it processed.
for r in pd.read_csv(T/'accepted_input_sha256.tsv',sep='\t').to_dict('records'):
 with open(r['path'],'rb') as f: assert hashlib.file_digest(f,'sha256').hexdigest()==r['sha256'],r['path']
def table(df):
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))
view=S[['sample','site','FibLike','dual_low_FibLike','legacy_rule_candidate','raw_strict_candidate','supported_candidate','positive_count_sensitivity','reference_QC_pass']].copy()
view['sample']=view['sample'].str.split('__').str[-1]
view.columns=['样本','部位','FibLike','双低','基础规则候选','加风险检查后','校正后支持','>0计数敏感性','参考QC']
reason=F.assign(exclusion_reasons=F.exclusion_reasons.fillna('')).exclusion_reasons.str.split(';').explode()
counts=reason[reason!=''].value_counts().rename_axis('原因代码').rename('细胞数').reset_index()
counts.to_csv(T/'overlapping_exclusion_reason_counts.tsv',sep='\t',index=False)
mism=A[A.caf_supported_candidate!=A.caf_supported_positive_sensitivity][['cell_id','sample','caf_supported_candidate','caf_supported_positive_sensitivity','raw_strict_candidate']]
mism.to_csv(T/'corrected_detection_threshold_sensitivity.tsv',sep='\t',index=False)
waterfall=[]
for sm,g in F.groupby('sample'):
 keep=pd.Series(True,index=g.index);reasons=g.exclusion_reasons.fillna('')
 gates=[('FibLike',pd.Series(True,index=g.index)),('dual_low',g.joint.eq('both_low')),('reference_QC',g.reference_QC_pass),('Fib_markers',g.raw_det_Fib>=2),('gene_complexity',g.n_genes>=500),('CTA',~reasons.str.contains('CTA_context_or_expression;',regex=False)),('mixed_lineage',~reasons.str.contains('mixed_lineage_fraction;',regex=False)),('competing_lineage',~reasons.str.contains('competing_lineage_program;',regex=False)),('osteo_program',~g.raw_osteo_conflict),('cluster_both_high',g.cluster_both_high_fraction.fillna(1)<.2),('cluster_any_high',g.cluster_any_high_fraction.fillna(1)<.2),('ambient_stability',g.corrected_lineage_supported)]
 previous=len(g)
 for step,mask in gates:
  keep &= mask
  n=int(keep.sum());waterfall.append({'sample':sm,'step':step,'remaining':n,'removed_at_step':previous-n});previous=n
 assert keep.equals(g.caf_supported_candidate),(sm,'waterfall conservation')
W=pd.DataFrame(waterfall);W.to_csv(T/'identity_gate_waterfall.tsv',sep=chr(9),index=False)

r1=S[S['sample'].str.endswith('GSM9504303')].iloc[0];r2=S[S['sample'].str.endswith('GSM9504304')].iloc[0]
total=int(S.supported_candidate.sum());raw=int(S.raw_strict_candidate.sum());legacy=int(S.legacy_rule_candidate.sum())
sensitivity_only=int(A.count_threshold_sensitive_candidate.sum())
trace_sensitive=int(A.trace_competing_sensitive_primary.sum())
both_rules=int(A.supported_under_both_detection_rules.sum())
lung_sensitive=A[(A['sample']=='GSE318840__GSM9504304')&A.count_threshold_sensitive_candidate]
byds=S.groupby('dataset')[['FibLike','dual_low_FibLike','raw_strict_candidate','supported_candidate']].sum().reset_index()
byds.columns=['队列','FibLike','双低','原始表达支持','校正后支持']
content=f"""# GEO 队列 CAF 联合身份确认报告

日期：2026-09-19。版本：geo_identity_20260919。范围：已验收的 16 份 GEO 样本，保留 101,487 个细胞的完整身份记录。

## 结论：{total:,} 个细胞满足联合证据支持的 CAF 候选定义

在 8,870 个既有 FibLike 细胞中，3,028 个具有双方法低 CNV 信号；加入谱系、参考校准、所在细胞群 CNV 风险与环境 RNA 敏感性检查后，{total:,} 个细胞获得联合证据支持。它们作为肿瘤组织中的成纤维细胞候选进入本轮交付，逐细胞名单与原始表达对象已导出。

复用既有基础门槛并重算新队列谱系证据后获得 {legacy:,} 个候选；加上骨样程序和任一方法群体高信号排除规则后为 {raw:,} 个；去污染后仍满足谱系规则的为 {total:,} 个。原始表达不支持的细胞不会因为校正结果而被提升为候选。另有 {sensitivity_only:,} 个细胞保留为计数检出阈值敏感的候选，单独交付，不将其解释为成纤维信号消失。

两份肺样本中，GSM9504303 最终保留 {int(r1.supported_candidate):,} 个，GSM9504304 最终保留 {int(r2.supported_candidate):,} 个。跨队列患者去重尚未确立，本轮不增加合并独立患者 N。

## Methods

输入为 geo_preprocess_20260918/accepted 中的分析就绪计数对象、geo_cnv_20260918 中的冻结细胞信息，以及 result_review_20260919 中经审计的逐细胞 CNV 判定。按 cell_id 核对 101,487 个细胞的顺序和身份，重现原始成纤维、骨样、周细胞和 CTA 标志基因检出数。原始细胞标签、计数、CNV 阈值和中期报告均保留。

严格候选以既有 FibLike 标签、至少两个成纤维标志、至少 500 个检出基因、inferCNV 与 CopyKAT 双低为基础。每一种已使用正常参考谱系及对应留出组均要求至少 10 个可判定细胞、至少 90% 低信号；逐组无法判定比例另列。按照此前规则，所在 cluster 的双高比例须低于 20%，分母为联合可判定细胞。

本轮在计算前增加两条身份守门规则：骨样程序达到至少 3 个标志且占计数超过 0.2% 时保留为混合/未决；所在 cluster 的任一方法高信号比例达到 20% 时不纳入严格候选，其分母为两方法可进行二分类的细胞。规则文件与原有口径结果均保存，不调整 CNV 判定阈值。

成纤维、骨样、周细胞、免疫、内皮、上皮、成肌及 CTA 基因集来源于既有谱系复核。竞争谱系检出数与计数占比沿用前期口径；破骨细胞使用专属标志的 0.1% 计数比例，其他竞争谱系使用 0.2%。免疫/内皮/上皮程序与间充质程序分别超过 0.2% 的细胞保留混合谱系标记。CTA cluster 或至少两个 CTA 基因均触发排除。完整基因表及各样本覆盖保存在 gene_panels.json 与 gene_panel_coverage.tsv。

对所有含 FibLike 的 {len(J)} 份样本，以全样本细胞矩阵分样本运行 DecontX；输入 z 为已复核的广义谱系标签，seed=20260919，使用默认 delta 与 estimateDelta，不传入空液滴背景。每次最多两个样本。按原始与校正矩阵计算同一套证据，校正后的标志检出主口径为计数≥1，同时保存计数>0的敏感性结果。校正不能恢复原始严格规则未通过的细胞，也不改变 CNV 结果。[DecontX 官方方法说明](https://github.com/campbio/decontX/blob/main/vignettes/decontX.Rmd)

所有分析为描述性身份判定；未进行细胞水平推断性检验、P 值或 FDR。样本来源、技术、患者映射和跨队列独立性状态随细胞保留。

## Results

### 1. 全部细胞得到可追溯的身份归属，原始输入保持一致

16 份样本、101,487 个细胞的 ID 顺序与计数核对通过。15 份样本已有 CNV，覆盖 97,238 个细胞；GSM9504299 的 4,249 个细胞继续标为未评估，其中 381 个 FibLike 未作为 CNV 阴性。非 FibLike 的细胞保留原谱系和 CNV 状态，本轮未扩大候选搜索范围。

{table(view)}

“基础规则候选”用于口径追踪；“加风险检查后”是本轮原始表达下的严格候选；“校正后支持”为交付集合。“>0计数敏感性”列应用相同逻辑但将校正计数的检出界定改为>0。各列是阶段性集合，不能相加。

### 2. 校正计数的检出阈值影响部分候选边界

{len(J)}/{len(J)} 份含 FibLike 的样本 DecontX 运行成功，退出码和完成标记均验证。共有 {raw:,} 个细胞在原始表达下满足本轮严格规则，其中 {total:,} 个在校正后仍满足成纤维谱系支持及竞争/混合谱系排除条件。校正前后每个 marker 的检出比例、均值以及估计污染比例保留在补充表中。

主口径获得 {total:,} 个，计数>0敏感性口径获得 {int(S.positive_count_sensitivity.sum()):,} 个；两种口径发生归属变化的细胞数为 {len(mism):,}。由于检出规则同时用于正向和竞争谱系证据，两者不保证嵌套。两种规则均通过的为 {both_rules:,} 个；主口径未通过但>0口径通过的为 {sensitivity_only:,} 个；主口径通过但>0口径出现额外竞争标志的为 {trace_sensitive:,} 个。后两类均有逐细胞名单，不能将差异一概归因于环境 RNA 去除。

{table(byds)}

### 3. 两份肺样本保留不同的身份结论

GSM9504303 的 1,025 个 FibLike 中，6 个为双低；此前已确认主要候选群存在两方法可见的染色体变化模式。本轮原始严格规则通过 {int(r1.raw_strict_candidate)} 个，校正后支持 {int(r1.supported_candidate)} 个。高信号、分歧及高风险群中的低信号细胞均保留原 CNV 证据，不能通过成纤维标志阳性覆盖异常线索。

GSM9504304 的 562 个 FibLike 中，431 个为双低。本轮原始严格规则通过 {int(r2.raw_strict_candidate)} 个，校正后支持 {int(r2.supported_candidate)} 个。其余双低细胞的复杂度、竞争谱系、参考校准或校正稳定性等原因逐细胞列出，避免将 431 个双低直接解释为身份确认数。

GSM9504304 的 {len(lung_sensitive)} 个检出阈值敏感细胞，其原始成纤维程序计数占比中位数为 {lung_sensitive.raw_pct_Fib.median():.3f}%，校正后为 {lung_sensitive.corrected_pct_Fib.median():.3f}%，估计污染比例中位数为 {lung_sensitive.decontx_contamination.median()*100:.2f}%。方向支持成纤维程序保留；本轮将其单列为敏感性候选，不按“校正后信号消失”排除。

两份肺样本均为 FLEX，保留技术验证开放标记；本轮不将其称为肺转移特异成纤维细胞或已验证临床亚型。

逐样本顺序筛选人数保存在 identity_gate_waterfall.tsv，所有筛选顺序固定，终点与候选名单逐细胞一致。\n\n### 4. 排除依据全部保留，未决细胞未被删除

以下排除原因可重叠，不能相加为总排除数。所有 8,870 个 FibLike 均保留原始判定和本轮证据，其他 92,617 个细胞也在全细胞表中。

{table(counts)}

## Discussion

本轮完成可由现有表达、CNV 与样本来源证据支持的操作性身份判定。满足规则者称为“联合证据支持的 CAF 候选”，不等于通过独立 DNA/病理真值确认非恶性，亦未建立肿瘤特异性。未决并不等于恶性，也不等于无成纤维细胞。

DecontX 无空液滴背景，以已观察谱系估计污染；校正前后稳定属于模型敏感性证据，不是独立实验验证。低深度、同谱系污染以及 FLEX 技术特点仍影响解释。没有为提高 CAF 数量调低检出或 CNV 阈值。严格候选集的缩小不能反向否定此前已经观察到的成纤维表达。

参考检查沿用既有每谱系口径；高无法判定比例仍见 reference_QC.tsv。本轮没有新增 CNV 参考重跑、分段 DNA 验证或修改已有 CNV 判定。因此低 CNV 信号仍不能完全排除低幅异常或其他恶性细胞。

患者与捕获映射按原 manifest 保存。GSE250015 与既有上海队列的重叠未解决；GSE318840、GSE162454 的跨队列独立性未通过本轮核实。不合并患者 N，不进行组间临床或转移推断。

## 交付与验收

- tables/supported_CAF_candidates.tsv：本轮支持的候选名单，带来源、技术和独立性状态。
- tables/sensitivity_only_CAF_candidates.tsv：计数检出阈值敏感候选及其保留的谱系证据。
- tables/candidate_sensitivity_tiers.tsv：主集合、两规则共同支持集合与敏感性集合的精确关系。
- tables/all_cell_identity_evidence.tsv.gz：101,487 个细胞的原始证据、校正证据、CNV、身份判定及全部排除原因。
- objects/*_supported_CAF_candidates.h5ad：逐样本候选对象；原始 counts 保留，decontx_counts 单独存储；坐标为此前初始显示坐标，不宣称新整合。
- tables/sample_identity_summary.tsv、reference_QC.tsv、cluster_CNV_risk.tsv、marker_before_after.tsv、ambient_summary.tsv：完整样本与方法证据。
- tables/sample_provenance.tsv、gene_panel_coverage.tsv、corrected_detection_threshold_sensitivity.tsv：来源、覆盖和敏感性材料。
- decision_rules.json、gene_panels.json、code/、logs/、SHA256SUMS.tsv：规则、代码、日志和校验记录。

![身份判定流程](figures/identity_funnel.png)

![校正后标志表达](figures/supported_markers.png)
"""
if not (O/'figures/supported_markers.png').exists():
 content=content.replace('![校正后标志表达](figures/supported_markers.png)','本轮无满足全部规则的候选，未生成候选标志图。')
(O/'CAF_GEO_identity_review_20260919.md').write_text(content,encoding='utf-8')
state['integrity_QA']='PASS';state['accepted_inputs_hash_recheck']='PASS';state['decontx_successful_samples']=len(J)
state['report']='CAF_GEO_identity_review_20260919.md'
(O/'STATUS.json').write_text(json.dumps(state,indent=2))
print(json.dumps(state),flush=True)
