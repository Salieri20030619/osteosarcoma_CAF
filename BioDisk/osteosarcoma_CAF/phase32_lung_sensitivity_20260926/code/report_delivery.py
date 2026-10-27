from pathlib import Path
import json,datetime
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase32_lung_sensitivity_20260926');S=O/'sensitivity_fixed_features'
def read(f):return pd.read_csv(f,sep='\t')
d=read(S/'tables/sensitivity_summary.tsv')
lines=['# 肺参考映射、分群敏感性与抗体方案复核','', '本轮保留支持层 1,224、扩展层 8,174 及 F2 内部 90/225 个候选细胞，未新增或删减客户阈值。','', '## 分群敏感性已完成，内部小群尚不支持稳定亚型命名','', '20/20 个固定原始特征的敏感性模型完成并通过独立核验。原始高变基因集合、细胞 ID、随机种子和分群参数固定；比较去污染、周期调整、应激调整及联合调整。原始分类未覆盖。','', '| 范围 | 去污染 | 周期调整 | 应激调整 | 应激＋周期 | 去污染＋应激＋周期 |','|---|---:|---:|---:|---:|---:|']
for layer,name in [('supported','支持层'),('expanded','扩展层'),('supported_myCAF','支持层 F2 内部'),('expanded_myCAF','扩展层 F2 内部')]:
 z=d[d.layer.eq(layer)].set_index('scenario');keys=['decontx','raw_adjust_cycle','raw_adjust_stress','raw_adjust_stress_cycle','decontx_adjust_stress_cycle']
 lines.append('| '+name+' | '+' | '.join(f'{z.loc[k,"ARI"]:.3f}' for k in keys)+' |')
lines+=['','表内为与原分群比较的 ARI，不设自拟“通过”阈值。扩展层整体对应关系较接近，F2 内部小群边界则明显依赖处理方式。','', '扩展层 F2 在最大对应新群中的保留比例为 80.4%–93.3%；该对应新群内原 F2 占比为 65.6%–87.9%。这些数值描述对应关系，不等于流式纯度或生物学亚型验证。','', '现阶段主线保留候选大类和连续功能程序。支持层原 5 个、扩展层原 6 个 myCAF 内部计算群继续标为探索性状态，不将其逐一命名为已确认的新亚型。','', '首轮重估高变基因遇到单细胞批次方差无定义的问题，诊断记录保留；正式对照统一使用原始特征集合，没有按批次细胞数删除候选。','', '## 肺参考映射','']
if (O/'delivery/HLCA_ACCEPTANCE_QA.json').exists() and (O/'HLCA_JOINT_EVIDENCE_COMPLETE').exists():
 qa=json.loads((O/'delivery/HLCA_ACCEPTANCE_QA.json').read_text());calls=read(O/'tables/HLCA_cell_predictions.tsv.gz')
 lines += [f'HLCA 映射覆盖全部 {len(calls):,} 个候选细胞，参考核心细胞 {qa["reference_core_cells"]:,} 个。已完成基因对应、模型投射、近邻权重和逐细胞概率完整性核验。','', '保留未经置信度过滤的标签及连续不确定度；以下为参考表型相似性，不是正常来源追踪、CAF 确证或肿瘤特异性。','', '主要研究 SCPCP000017 的患者级参考权重如下。每位患者先在类别内取平均，再对患者取平均；支持层与扩展层重叠，不作为独立验证。','', '| 层/类别 | 患者 N | AdvF | AlvF | Peribronchial | Myofibroblast | Pericyte | Smooth muscle |','|---|---:|---:|---:|---:|---:|---:|---:|']
 summary=read(O/'tables/HLCA_class_study_summary.tsv')
 for (layer,cls),g in summary[summary.dataset.eq('SCPCP000017')].groupby(['layer','class_id']):
  row=g.set_index('phenotype');values=[row.loc[p,'mean_patient_weight'] if p in row.index else 0. for p in ['Adventitial fibroblasts','Alveolar fibroblasts','Peribronchial fibroblasts','Myofibroblasts','Pericytes','Smooth muscle']]
  lines.append('| '+layer+' '+cls+' | '+str(int(g.donor_records.iloc[0]))+' | '+' | '.join(f'{v:.3f}' for v in values)+' |')
 lines += ['', '其余参考类别的权重保留在完整统计表中，未强制把所有候选映射为成纤维细胞。骨肺配对的细胞数、未经置信度过滤的标签比例和连续权重均已另存；分母为零时比例缺失，不当作零。']
else:lines += ['模型或联合验收仍未完成，此处不填入预测比例。']
lines+=['','## Gong 参考面板','', '已核实并计算原文人源 PTGS2 共表达面板：PTGS2、IL6、CCL2、CXCL1、PTX3、TNFAIP6，覆盖 8,174 个候选细胞，保留原始和去污染后连续评分。此六基因面板不是鼠源完整差异签名，也没有设“表达两个基因”等阳性门槛。来源：[Gong 等，Immunity 2022，Figure 3H 对应结果](https://pmc.ncbi.nlm.nih.gov/articles/PMC9830653/)。','', '## 抗体方案已具备分析验证设计，分选特异性尚未成立','', '已核实 PDGFRB 28D4、MCAM P1H12、THY1 5E10、FAP 427819 的人体流式产品；阴性排除和 CSPG4 壁细胞重叠观察通道、条件性双管配色及实验对照见独立抗体方案。','', '关键问题是全细胞背景的 RNA 模拟双阳性门特异性弱：PDGFRB＋MCAM 模拟门中 F2 占比中位数约 0.31%，CNV 异常/冲突比例在支持/扩展口径中约 88.4%/94.0%。这不是实际流式纯度预测，但不能忽略。保留至少两个阳性表面标志的客户要求，同时通过实际蛋白分析和恶性/壁细胞身份复核评估组合。','', 'FAP 不强制作为全部 F2 的阳性入门条件，CSPG4 也不自行设为一票排除标志。客户未给定的阳性强度、覆盖率和纯度阈值不自拟。','', '## 交付范围','', '本轮结果来自既有冻结数据范围；未导入的外部人源队列不宣称已映射。替代 CNV 参考重跑仍是全项目独立未结项，不能用本次应激/周期分群敏感性替代。流式滴定、消化表位及分选后身份/纯度需要真实实验数据，当前交付不称为已验证分选方案。','', '计算方法来源：[scvi-tools HLCA 官方映射教程](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/hub/query_hlca_knn.html)。']
(O/'delivery/肺映射_敏感性_抗体综合复核_20260926.md').write_text('\n'.join(lines)+'\n')
print('INTEGRATED_REVIEW_WRITTEN')
