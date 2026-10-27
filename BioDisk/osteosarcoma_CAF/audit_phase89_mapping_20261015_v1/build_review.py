from pathlib import Path
import json, hashlib, datetime

P=Path(__file__).resolve().parent
A=json.loads((P/'AUDIT_COMPUTATION.json').read_text(encoding='utf-8'))
assert A['all_checks_pass']
stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
zh='''CAF Phase89 固定状态映射审核与进度衔接
审核版本：phase89_fixed_state_mapping_20261015；MAPPING_COMPLETE.json 完成时间 2026-10-15T02:07:32.851096+00:00。

结论：通过本次约定范围的文件完整性、计数、记录判定规则和探索性描述使用核查。没有发现必须修复后才能保留当前探索性结果的实现错误。此通过不涵盖可靠生物学身份确认、独立外部验证或正式临床合并。
本次审核没有改动 Phase89 代码、映射表、图、阈值或旧基线；只生成独立审核派生文件并更新进度元数据。

已确认的进度
1. 原 8,174 个细胞、33 个状态保留；映射模板排除 26 个已有软件双细胞阳性。保留全部 1,268 个审核后的拓展探索候选、全部上游列和风险标记。
2. 1,268 = 状态相容 534 + 仅大类相容 242 + 未定 492；51 个样本/群组合。1,040 个参考校准通过，228 个参考校准失败。11 份保留候选的文库对应 10 个来源内患者/供者标签；原 Phase88 有 12 份文库，EXT_OS_1_L1 没有进入本轮候选。
3. 状态相容分别为 F1_originalC10 178、F1_originalC4 6、F2_C5 26、F4_C10 30、F5_originalC1 1、F5_originalC12 293。全部 6 个 F1_originalC4 与 26 个 F2_C5 候选参考校准失败；其余 502 个状态相容候选参考校准通过。映射相容与 CNV 参考通过是不同证据。
4. 留出来源记录有 283 个 donor/state 表达谱，模板状态与原标签一致 86 个，大类一致 184 个。不是 283 个独立患者。
5. 审核新增描述性诊断：把既定群级接收规则应用到已保存留出预测，100 个获判状态相容的表达谱中 57 个与原状态一致；106 个仅大类相容表达谱中 73 个与原大类一致。阈值和可校准状态取自同一批留出记录，因此这是重复使用校准数据的内部描述性诊断，不是独立测试，不是 534 个新细胞的准确率，也不复现新细胞的完整群级加单细胞联合接收过程。未作推断检验。
6. F1_originalC10 内部接收 4 个、原标签一致 4 个，但样本很少，新候选 178 个仅来自 BS7；不能据此称为稳定外部复制。F5_originalC12 内部接收 28 个、原标签一致 17 个；其参考池 466 个细胞全部原标为探索性，其中身份未定 269、supported_candidate 197，谱系分层为 Other_unresolved 341、Osteo 91、Fibro 33、Mural 1。新映射 293 个不能据此确认新的纯 CAF 免疫亚型。
7. 32 个状态相容单细胞的 margin 低于群级校准阈值。这符合现有协议：阈值用于群平均谱，单细胞需两种方法与群状态一致。未发现协议与实现不符；不能把群级阈值描述成已验证的单细胞置信度。

问题与证据边界
当前探索性范围无无条件阻断项。以下区分升级结论前条件与可选补强；不新增候选筛选硬门槛。

CAF-MAP-01 | P1 | 升级为可靠状态判定前必须解决；当前探索使用已记录边界
位置：code/map_states.py:138-150、156-168；tables/state_transfer_validation.tsv；MAPPING_COMPLETE.json.validation；审核 audit_tables/heldout_predictions_with_mapping_rule.tsv、state_diagnostics.tsv。
证据：校准只取模板预测正确记录的分位数；原摘要提供 recall 和正确计数，未直接给最终接收规则的错配分布。审核得到状态相容原标签一致 57/100、大类相容 73/106；当前单细胞联合接收过程没有独立准确性评估。
影响：state_compatible 可以表示操作性相容，不能升级为可靠身份或已验证新亚群；上述内部一致率不能外推到新细胞。
建议：保留 534/242/492 和全部风险列，解释时附本次诊断及参考状态证据层级。若将标签用于确认亚群或状态分组推断，另行版本化定义评价目标，并分离阈值校准与评价，检验误接收、未知类型和相应群/单细胞使用层级；不按目标细胞数调阈值。
验收：探索版报告可追溯 57/100、73/106、283 的单位和分母，并明确同批校准及参考标签非生物学真值；更强结论须有相应独立评价与预设可接受标准。原候选与标签保留，不自动剔除或改名。

CAF-MAP-02 | P2 | 可选的复现性补强；不阻断当前只读探索审核
位置：code/map_states.py:39-40、46-50、62-68、175；resources/INPUT_HASHES.json。
证据：原始前后哈希只覆盖 4 个对象；查询 analysis_ready.h5ad、特征存在性矩阵及样本索引、PROGRAMS.json 不在该清单。4 个已登记对象本次均重新核对一致。
影响：原始运行时全输入未变的证据范围有限，不能将其扩大成全输入封存证明；没有证据表明漏登对象实际发生了改变。
建议：今后版本补齐全部实际读取输入、代码和环境登记。现在补记的哈希只能称审核时快照，不能追溯宣称运行前后相同。
验收：如采纳，清单含路径、时间、哈希及记录阶段；原始 4 对象清单和旧结果保持不变；明确历史运行时未验证范围。

患者独立性与取材信息衔接（已完成现有资料层核对，未完成跨研究去重）
核对已保存 NCC 作者 Patient metadata、CELLxGENE 作者 metadata 核查及 Phase88 patient_sample_roster.tsv。NCC BS4/BS14 对应同一 Patient4；EXT_OS_1 两文库对应同一供者。保留候选内 8 个 NCC 患者标签、2 个 CELLxGENE 供者标签，不等于已确认新增 10 个独立患者。
NCC 部位仍 not_explicitly_verified；CELLxGENE 仍 limb。Primary/Recurrent、肢体位置或组织学 Fibroblastic 均不自动替代 Bone/Lung 取材证据。EXT_OS_1 只剩 1 个参考失败、仅大类相容细胞，不能据此宣布新增可用 CAF 患者。
本轮无新遗传去重证据、病理/手术部位记录或跨研究对应证据；新增正式患者及骨肺比较患者仍为 0。audit_tables/source_patient_inventory.tsv 按来源/患者标签列出全部 1,268 个候选，不使用细胞数代替推断 N。

沿用原问题编号的复核状态
CAF-EXP-01 | P1 | 未解决：预处理口径与原生 CNV 可比性未在 Phase89 改变；正式同口径合并前提供逐级计数、原生调用与描述性信号独立列、方法版本和基线不变证明。
CAF-EXP-02 | P1 | 部分解决：1,268 个群级审核候选已原样用于独立映射，1,021 个竞争背景未混入；尚未正式纳入。验收仍为 2,289=1,268+1,021、116 群依据和原标识符保留、1 个严格层分歧单独记录；这些已由前轮复核，未重复扩查。
CAF-EXP-03 | P1 | 部分解决：来源内患者/文库对应已核对；跨研究独立性及样本取材仍未解决。验收：用于临床比较的每个样本具可追溯部位证据及患者对应，BS4/BS14、EXT_OS_1 不重复计 N；无依据样本不进入 Bone/Lung 推断。
CAF-EXP-04 | P2 | 未解决但分层保留：1,040 通过、228 失败，另 728 个信号支持候选未进入映射。验收：原调用与全部风险列保留，不把失败/不可评估写成正常；新原生证据须独立版本，不以 1 个未定细胞称新增患者。

实际核查范围与未验证事项
已检查完整 Phase89 脚本、协议、完成标记、全部 1,268 行和 51 群判定、283 条留出预测、校准分位数、邻居来源与供者唯一性、全部上游列保持、登记输入哈希、33 状态摘要、患者对应和图的可读性。分位数回读最大误差约 1.33e-9，符合 float32 与 TSV 精度差异；不是算法错误。图中标签可读，无裁切，图只展示 776 个大类相容细胞。
未重新从表达矩阵计算全部相似度/特征/函数评分，未独立复跑训练，未开展外部验证、原生 CNV 重跑或新推断检验。Phase89 原始进程退出码未在所审目录独立保存；本次依据完成标记、代码与输出核查，不把进程消失作为成功。审核脚本实际退出码另行记录。
下一步：可基于现有候选进行来源分层的探索性解释；若要正式扩充，先解决 CAF-EXP-01/03/04 的方法、部位和独立性证据。若要确认亚群身份或依赖映射标签作推断，先满足 CAF-MAP-01。不得静默改动当前映射或原始基线。
'''
en='''CAF Phase89 review feedback for Claude Code
Please respond in English. This is review feedback, not an authorization to modify artifacts or launch repairs automatically.
Reviewed version: phase89_fixed_state_mapping_20261015, completed 2026-10-15T02:07:32.851096+00:00.

Disposition: PASS for the inspected file integrity, counts, recorded decision rules and descriptive exploratory use. No unconditional implementation blocker was found for retaining the current exploratory outputs. This is not approval of validated biological identity, independent external validation or formal clinical admission. Original code, results, thresholds and baseline were not changed.

Confirmed: 8,174 baseline cells and 33 states retained; 26 existing doublet-positive cells excluded from templates. All 1,268 reviewed candidates and upstream fields/risk flags preserved. Counts are 534 state-compatible, 242 class-only-compatible and 492 unresolved, in 51 sample/cluster groups; reference calibration passes for 1,040 and fails for 228. Retained candidates occupy 11 libraries and 10 source-scoped patient/donor labels; EXT_OS_1_L1 contributes no retained candidate. State-compatible counts: F1_originalC10 178; F1_originalC4 6; F2_C5 26; F4_C10 30; F5_originalC1 1; F5_originalC12 293. All 32 F1_originalC4/F2_C5 cells have failed reference calibration; the other 502 pass.

CAF-MAP-01 | P1 | Required before upgrading to reliable state assignment; boundary documented for current exploratory use
Location: code/map_states.py:138-150,156-168; tables/state_transfer_validation.tsv; MAPPING_COMPLETE.json.validation; audit_tables/heldout_predictions_with_mapping_rule.tsv and state_diagnostics.tsv.
Evidence: thresholds use quantiles of template-correct held-out predictions. The original summary records recall/correct counts but not errors among assignments accepted by the deployed group rule. Applying that rule to the saved held-out profiles yields original-state agreement in 57/100 state-compatible profiles and original-class agreement in 73/106 class-only profiles. There are 283 donor/state profiles, not 283 independent patients. Thresholds and calibratable states use those same records, so this is a reused-calibration internal descriptive diagnostic, not an independent test, not accuracy for the 534 new cells, and not an evaluation of the full new group-plus-single-cell pipeline. No inferential test was performed.
Impact: state_compatible denotes operational compatibility; it does not establish validated identity or a new subtype. Internal agreement must not be extrapolated to new cells.
Recommendation: preserve 534/242/492 and all risk flags, and accompany interpretation with the audit diagnostic and reference-state evidence levels. Before subtype confirmation or inference based on transferred states, separately version the evaluation target, separate calibration from evaluation, and assess false acceptance, unknown types and the relevant group/cell level. Do not tune thresholds to a target cell count.
Acceptance: exploratory reporting must reproduce 57/100, 73/106 and the unit/denominator of 283, disclose reused calibration and non-ground-truth reference labels; stronger claims require corresponding independent evaluation and prespecified acceptance criteria. Preserve candidate IDs and labels; do not automatically exclude or rename them.
Supporting context: F1_originalC10 has internal accepted agreement 4/4 but only four profiles, and all 178 query cells are from BS7. F5_originalC12 has agreement 17/28; its reference consists of 466 exploratory cells (269 identity-unresolved, 197 supported_candidate), with lineage strata Other_unresolved 341, Osteo 91, Fibro 33, Mural 1. Its 293 query matches do not confirm a pure CAF immune subtype. Thirty-two state-compatible query cells fall below the group-level margin cutoff. This follows the written protocol, which thresholds group means and checks cell-level agreement; it is not an implementation mismatch or validated cell-level confidence.

CAF-MAP-02 | P2 | Optional reproducibility improvement; not a blocker for current read-only exploratory review
Location: code/map_states.py:39-40,46-50,62-68,175; resources/INPUT_HASHES.json.
Evidence: the original before/after manifest covers four objects. Query analysis_ready.h5ad files, sample feature-presence matrix/index and PROGRAMS.json are absent. All four registered objects were rehashed and match.
Impact: historical input immutability evidence has limited coverage; it cannot establish that all inputs were frozen. There is no evidence that an omitted input actually changed.
Recommendation: future versions should register all read inputs, code and environment. Hashes collected now are audit-time snapshots, not proof of original-run before/after identity.
Acceptance if adopted: paths, timestamps, hashes and capture stage are explicit; preserve the original four-object manifest and results, and disclose historically unverified scope.

Patient metadata handoff: existing archived NCC author Patient metadata, CELLxGENE author metadata check and Phase88 patient_sample_roster.tsv were reviewed. BS4/BS14 are one Patient4; both EXT_OS_1 libraries are one donor. Eight NCC and two CELLxGENE source-scoped labels are not ten confirmed newly independent patients. NCC site remains not_explicitly_verified; CELLxGENE tissue remains limb. Primary/Recurrent, limb location and Fibroblastic histology are not substitutes for Bone/Lung specimen evidence. EXT_OS_1 has one reference-failed class-only-compatible cell, insufficient to claim a usable new CAF patient. No genetic deduplication, pathology/surgery site evidence or cross-study identity evidence was added. New formal patients and Bone/Lung patients remain zero. The source-scoped inventory is in audit_tables/source_patient_inventory.tsv.

Prior IDs retained:
CAF-EXP-01 | P1 | UNRESOLVED. Preprocessing/native-CNV comparability unchanged. Before comparable formal merging, provide stepwise counts, separate native calls/descriptive signals, method versions and baseline immutability.
CAF-EXP-02 | P1 | PARTIALLY RESOLVED. Exactly 1,268 adjudicated candidates are used in the separate mapping and 1,021 competing-background cells are excluded; no formal admission. Acceptance remains 2,289=1,268+1,021, documented decisions for 116 groups, preserved IDs and a separate record of the one strict-tier discordant cell. These were checked in the prior review and not unnecessarily reaudited.
CAF-EXP-03 | P1 | PARTIALLY RESOLVED. Source-scoped patient/library mappings checked; cross-study independence and specimen sites remain unresolved. Acceptance: every clinical-comparison sample has traceable site and patient evidence; BS4/BS14 and EXT_OS_1 are not double-counted; unsupported sites do not enter Bone/Lung inference.
CAF-EXP-04 | P2 | UNRESOLVED, STRATIFIED. The 1,040 passed and 228 failed reference flags persist; the separate 728 signal-only candidates are not mapped. Acceptance: retain original calls and risks, do not relabel failed/unevaluable as normal; independently version new native evidence and do not claim a new patient from one unresolved cell.

Scope: full mapping code/protocol/completion marker; all 1,268 output rows, 51 group decisions, 283 held-out predictions, calibration quantiles, neighbor source exclusion/donor uniqueness, all upstream fields, registered input hashes, 33-state summaries, patient metadata and figure readability. Quantile round-trip error is at most about 1.33e-9, consistent with float32/TSV precision. The figure is legible and unclipped and displays 776 class-compatible cells only.
Not independently rerun: expression-to-similarity/feature/function computations, model training, external validation, native CNV or inferential testing. An original Phase89 process exit-code record was not found in the reviewed directory. This review uses recorded completion, code and output checks, not process disappearance. The separate audit exit code is recorded independently.
Next: source-stratified exploratory interpretation is permissible within these boundaries. Formal expansion still requires CAF-EXP-01/03/04 evidence; subtype confirmation or inference using transferred labels requires CAF-MAP-01. No silent changes to current mappings or baseline.
'''
for name,s in [('REVIEW_ZH.txt',zh),('CC_FEEDBACK_EN.txt',en)]:
    (P/name).write_text(s,encoding='utf-8')
manifest={str(f.relative_to(P/'input_phase89')).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in (P/'input_phase89').rglob('*') if f.is_file()}
(P/'REVIEWED_PHASE89_HASHES.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
acceptance=dict(status='PASS_EXPLORATORY_REVIEW_ONLY',reviewed_version='phase89_fixed_state_mapping_20261015',reviewed_completion_utc='2026-10-15T02:07:32.851096+00:00',review_utc=stamp,computational_consistency_pass=A['all_checks_pass'],audit_script_exit_code=0,original_phase89_exit_code_record_available=False,new_inferential_tests=False,original_mapping_modified=False,formal_baseline_changed=False,new_formal_patients=0,new_Bone_Lung_patients=0,patient_metadata_review='existing_archived_source_metadata_checked_cross_study_identity_and_site_unresolved',cell_counts=A['query_counts'],reference_pass=1040,reference_fail=228,internal_reused_calibration_diagnostic={'state_agree':57,'state_accepted':100,'class_only_agree':73,'class_only_accepted':106,'profiles':283,'independent_validation':False,'new_cell_accuracy':False},issue_status={'CAF-MAP-01':'documented_for_exploratory_use_open_before_stronger_assignment_claims','CAF-MAP-02':'optional_reproducibility_improvement','CAF-EXP-01':'unresolved','CAF-EXP-02':'partially_resolved','CAF-EXP-03':'partially_resolved','CAF-EXP-04':'unresolved_stratified'},visual_review='labels_legible_no_clipping_776_class_compatible_cells_only')
(P/'ACCEPTANCE.json').write_text(json.dumps(acceptance,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':acceptance['status'],'files':['REVIEW_ZH.txt','CC_FEEDBACK_EN.txt','ACCEPTANCE.json','REVIEWED_PHASE89_HASHES.json']},ensure_ascii=False))
