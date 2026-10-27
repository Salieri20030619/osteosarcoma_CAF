from pathlib import Path
from datetime import datetime,timezone
import json,shutil,hashlib

root=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
out=root/'phase85_functional_objective_audit_20261012'
out.mkdir(parents=True,exist_ok=True)
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
prior=json.loads((root/'CAF_ALL14_RESEARCH_CURRENT.json').read_text(encoding='utf-8'))
assert len(prior['retained_genes'])==14
report='''CAF 研究主线澄清与现有分群诊断（2026-10-12）

真正研究终点
寻找具有 ECM 构建/重塑、免疫调节等功能特征，并与肺病灶、转移或复发相关的 CAF 状态，进一步寻找能够通过抗体流式识别的表面标志组合。功能状态、临床关联和检测可行性分别提供证据。

本次核查结论
现有分析已有内部状态标志和部分患者级名义差异，尚未得到跨患者支持充分、具有目标功能且能用流式区分的 CAF 亚群。不能把这个结果概括成所有 CAF 均无差异，也不能保证重新分群后一定获得目标亚群。
1. 功能与原 F 类并非一一对应。原报告 F1_C10 为 SFRP2/PLA2G2A 基质与分泌相关态，含 THY1、LOX；F2_C3 为 CD248/POSTN 胶原基质相关态，含 VCAN、MMP9。F4 内也有不同基质相关状态。因此目标功能跨越多个旧类，粗类比较存在混合状态的风险。以上为既有描述画像，不代表其与肺转移或耐药的联系已建立。
2. 原整合报告记录：33 个内部簇中有 19 个至少 90% 的细胞来自同一患者记录。继续细分并不必然增加跨患者的生物学信息。
3. 原扩展 F4 含 Fibro-like 909、Osteo-like 3300、其他/未解决 1095 个细胞。它不能直接作为单一纯 CAF 群解释。此处是既有谱系判断，不等于将全部 Osteo-like 细胞重新判为恶性。
4. 已有 S17 扩展 F1 骨/肺比较为 8/3 名独立患者，其中肺侧两人分别仅有 2、1 个 F1 细胞；扩展 F4 为 11/2，其中肺侧一人只有 1 个 F4 细胞。分群后的表达代表性受限，不能用大量细胞替代独立患者数。
5. 既有 14 个基因是合并 CAF 候选细胞后两来源同向的部位差异候选；它们尚非一个特定功能亚群的表面标志集合。

后续分析方案（本轮登记，尚未执行新的重分群或功能评分）
A. 原 F1–F5 和细胞编号保留。使用现有 CNV、谱系、双细胞证据，在较可信的成纤维细胞中分析功能状态，并用完整候选池作敏感性分析。身份不确定细胞单独标识，不凭 ECM 高表达确认 CAF 身份。
B. 在不知道病灶分组结果的前提下预先定义 ECM 合成/装配、ECM 重塑/交联、收缩、炎症/免疫调节等程序；先以连续评分分析其在旧亚群内外的分布。状态可以重叠；若没有稳定离散边界，不强行命名新亚群。不按显著 P 值选择聚类分辨率。
C. 区分两类问题：功能状态相对其他 CAF 的标志差异；同一功能状态在不同病灶的丰度或活性差异。亚群识别标志不必同时是骨—肺差异基因。不要求一个基因同时满足两种不同的用途。
D. 寻找表面标志组合时，不把用于定义状态的同一基因当作独立验证；采用未参与状态定义的表面基因及留出患者/可用独立来源评估，报告覆盖度及非 CAF 背景表达。RNA 差异、表面蛋白表达和抗体流式可用性分别核查。
E. 主统计来源仍为 SCPCP000017 与 SCPCP000023。SCPCP000018、GSE152048、GSE318840 依用户决定仅作记录，不参与候选保留或优先级判断。原报告的三个主要来源中包括 S18，这不等于恢复它的判定资格。若未来需要更改来源角色，须明确记录新的依据和范围。
F. 按独立患者汇总功能分数与状态比例。有效比较按目标原始 P < 0.05 筛名义信号，同时保留效应、样本一致性、FDR 和检验家族；重复不足仅作描述。文献决定研究假设与机制依据，不代替本项目临床检验。
G. 骨—肺解剖部位比较只能回答部位关联。复发、转移史、治疗与耐药需要相应患者标签或结局，缺失则记录未评估，不以部位代替。ECM 转录高表达不直接证明物理包绕、组织增厚、免疫抑制或耐药；组织结构和功能效应需另行验证。
H. 全部 14 个候选继续保留，映射其与功能状态的关系；后续表面标志搜索开放至全基因范围。lncRNA 可用于机制研究，但不能用普通抗蛋白抗体直接作活细胞表面流式标志。

与流式相关的文献参照
Kieffer et al., 2020：在其他肿瘤背景中通过流式验证功能不同的 CAF 簇，支持“功能状态—标志组合”的研究路线，而非本项目验证结果。
https://pubmed.ncbi.nlm.nih.gov/32434947/
Krishnamurty et al., 2022：胰腺癌小鼠模型中 LRRC15 阳性 myCAF 的免疫抑制功能有干预证据，并有流式识别。
https://pmc.ncbi.nlm.nih.gov/articles/PMC9630141/
LRRC15 在骨肉瘤细胞系也有表达，不能直接借用为骨肉瘤 CAF 特异标志。
https://pubmed.ncbi.nlm.nih.gov/32902907/
因此目前尚未形成经过本项目验证的流式 panel；不能仅靠一个 CAF 常用阳性标志或 CD45/CD31 阴性门确认骨肉瘤 CAF。

核查依据
本地已提取的《F1–F5 内部亚群表达、功能与置信度报告》（2026-10-05）、《骨肉瘤 CAF 功能分类一致性分析报告》整合文本，以及 Phase73/74 F1/F3/F4 结果和患者覆盖记录。本轮为既有结果诊断和目标登记，未重新计算统计显著性。
'''
(out/'CAF功能亚群研究_当前问题与后续方案_20261012.txt').write_text(report,encoding='utf-8')
record=dict(phase=out.name,status='objective_clarified_existing_results_audited',recorded_at=stamp,objective='功能相关且可流式识别的 CAF 状态及其临床关联',retained_genes=prior['retained_genes'],all14_retained=True,record_only_sources=prior['selection_policy']['record_only_sources'],primary_statistical_sources=prior['selection_policy']['primary_statistical_sources'],new_state_analysis_status='planned_not_run',new_flow_panel_status='not_validated',report=str(out/'CAF功能亚群研究_当前问题与后续方案_20261012.txt'))
(out/'OBJECTIVE_AND_AUDIT.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
(out/'state_backups').mkdir(exist_ok=True)
for name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_ALL14_RESEARCH_CURRENT.json']:
    p=root/name;shutil.copy2(p,out/'state_backups'/f'{stamp}_{name}')
    obj=json.loads(p.read_text(encoding='utf-8'))
    obj['current_research_objective']=record
    obj['next_analysis_plan']='跨旧亚群功能状态分析及表面标志搜索；尚未执行'
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(root/'CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(status=record['status'],report=report,record=record),ensure_ascii=False))
