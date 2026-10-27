from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
O=Path('C:/codex_test/CAF_discovery_stats_20260921');B=Path('C:/codex_test')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t')
R=read('discovery_program_tests.tsv');M=read('within_myCAF_program_tests.tsv');D=read('C7_vs_C9_direct_tests.tsv');G=read('discovery_marker_tests.tsv');E=read('primary_patient_program_differences.tsv')
assert len(R)==260 and (R[R.layer=='supported'].shape[0]==90)
assert (R.loc[R.P.notna(),'P']>=R.loc[R.P.notna(),'minimum_attainable_P']-1e-12).all()
assert E.groupby(['layer','cluster','program','donor_record']).size().max()==1
assert set(E.dataset)=={'GSE152048','GSE162454','GSE318840'}
def fmt(v):return 'NA' if pd.isna(v) else f'{v:.6g}'
def table(df,cols,headers):
 rows=['|'+'|'.join(headers)+'|','|'+'|'.join(['---']*len(cols))+'|']
 for r in df.to_dict('records'):rows.append('|'+'|'.join(fmt(r[c]) if isinstance(r[c],float) else str(r[c]) for c in cols)+'|')
 return '\n'.join(rows)
f=R[(R.layer=='supported')&R.cluster.isin([7,9])].sort_values(['cluster','program']);m=M[(M.layer=='supported')&M.cluster.isin([7,9])].sort_values(['cluster','program'])
lines=['# 发现队列支持C7与C9具有不同的表达程序','',
'本轮主结论是发现队列中候选群的特征表达存在差异，不以外部结局验证显著作为必要条件。患者级分析支持C7更偏myCAF程序、C9更偏iCAF程序；两群均有矿化相关表达，但当前直接比较不支持把它们称为成骨或矿化程度明显不同的两个亚型。','',
'## C7与C9的直接患者内比较支持不同程序','',
'在10例两群均有足量细胞的患者中，C7相对C9的myCAF程序中位差为+0.337 logCP10k，P=0.00977，FDR=0.0244；iCAF程序中位差为−0.388，10/10例均为C9更高，P=0.00195，FDR=0.00977。成骨与矿化程序的两群直接差异未显著（P分别0.922、0.846）。',
table(D[D.population=='all_candidate_cells'],['program','n_patients','C7_higher_patients','C9_higher_patients','median_C7_minus_C9','P','FDR'],['程序','患者数','C7较高','C9较高','C7减C9中位差','P','FDR']),'',
'上述方向与“存在不同表达状态”一致。它不等同于已发现独立新谱系、已证实与转移相关，或已获得可用流式分选群。','',
'## 相对其余候选细胞的程序差异得到统计支持','',
'C7：myCAF程序在13例患者中12例升高，中位差+0.661，P=0.00122，FDR=0.00439；成骨程序差异较小（中位+0.0795）但方向一致，P=0.000488，FDR=0.00220；矿化程序13/13例升高，中位+0.117，P=0.000244，FDR=0.00146。',
'C9：iCAF和矿化程序均在15/15例升高，P均为0.0000610，FDR均为0.000687；myCAF也较群外升高。成骨程序为名义显著（P=0.0413、FDR=0.0790），记为探索性支持，不写成校正后显著。',
table(f,['cluster','program','n_patients','positive_patients','median_delta','CI_low','CI_high','P','FDR'],['群','程序','患者数','正向患者','中位差','中位差CI下限','CI上限','P','FDR']),'',
'差值均为同患者的群内平均程序分数减群外平均程序分数；置信区间为按研究分层的患者bootstrap中位差95%区间。表中apCAF为负向差异，不能将其P值显著写成apCAF富集。中位差区间与Wilcoxon检验的统计对象不同，可能出现区间跨零而秩检验名义显著。','',
'## 限定myCAF内部后仍保留可推进的差异','',
'在仅保留客户原始面板阳性myCAF细胞后，C7相对其余myCAF的矿化程序达到校正后显著（n=12，P=0.00342，FDR=0.0385）；C7的myCAF增强为名义显著（P=0.0342，FDR=0.176）。C9的iCAF程序达到校正后显著（n=14，P=0.000366，FDR=0.00659），矿化程序为名义显著（P=0.00671，FDR=0.0671）。两群在myCAF内部的成骨程序差异均未显著。',
table(m,['cluster','program','n_patients','positive_patients','median_delta','P','FDR'],['群','程序','患者数','正向患者','中位差','P','FDR']),'',
'C7与C9在myCAF内部可直接比较的患者为8例：7/8例C7的myCAF程序更高，中位差+0.361；8/8例C9的iCAF程序更高，中位差−0.305。该部分按小样本探索呈现，完整P与FDR如下。',
table(D[D.population=='within_myCAF'],['program','n_patients','C7_higher_patients','C9_higher_patients','median_C7_minus_C9','P','FDR'],['程序','患者数','C7较高','C9较高','C7减C9中位差','P','FDR']),'',
'## 敏感性分析支持主方向并显示适用范围','',
'将每位患者群内/外至少10个细胞提高到至少20个，C7的myCAF、成骨、矿化，以及C9的iCAF、矿化仍达到各自分层家族的FDR显著。逐一去掉一个主研究后，C9的iCAF与矿化在每次剩余10例中仍为FDR显著；C7方向保持，但去掉GSE152048后只剩7例，其myCAF和矿化保留名义显著，FDR未过。',
'分布假设更弱的精确符号检验另列：C7对其余候选的myCAF为P=0.00342，但C7对C9直接myCAF比较只有8/10例同向，符号检验P=0.109，表明该直接差异的显著性依赖效应大小信息。C9的iCAF对C7为10/10例同向，两种检验均支持。完整敏感性结果保留，不以更有利检验替换主检验。','',
'## Methods','',
'保持Phase17的三层群标签与Phase18的客户程序定义。每位患者只贡献一个群内减群外的程序差值；主纳入要求每侧至少10个细胞、面板覆盖至少80%。主层为supported，严格和广义拓展层为敏感性。全52个群及5个程序均计算，未只挑选C7/C9进行多重校正。',
'主统计池使用GSE152048、GSE162454和GSE318840。原始论文分别报告上海第六人民医院11例、广西医科大学6例、河内Vinmec8例的独立招募；依据这些独立采集来源及项目内患者映射，构成本轮主统计范围。不同比较的实际可评估患者数分别计数。其余来源保留既有描述，未将尚未解决复用或别名对应的来源加入本轮合并N。',
'采用双侧精确Wilcoxon配对符号秩检验，差值保留12位小数，对并列绝对秩以整数化秩分布动态规划求精确P；零差剔除。无有效生物学重复时仅保留描述，不作推断检验。效应为患者差值的中位数、方向一致数及按研究分层5000次bootstrap置信区间。另计算精确双侧符号检验。',
'主supported层的18群×5程序构成90项BH校正家族；strict为50项、expanded为120项的敏感性家族。不可评估项以P=1占位参与家族校正，交付表显示NA。myCAF内部比较另构成各层同样大小的敏感性家族。C7与C9直接比较在每个细胞范围内按5个程序校正。描述性marker候选基因使用另一个按层划分的家族，不与程序检验混合。','',
'## Discussion','',
'这些统计说明已选发现群的特征差异获得探索性统计支持。聚类、选群、marker发现和本轮检验来自同一表达数据，P/FDR未校正聚类与候选选择过程，不能当作预先指定标签的盲验证，也不存在一个通用P值可证明“新群本身成立”。独立性依据为公开招募来源和患者映射，而非实名核验；完整全8来源去重仍未完成。',
'当前最合适的主结论为：发现队列识别出特征不同的候选状态，C7偏myCAF表达，C9偏iCAF表达，并存在矿化相关特征；在myCAF内部仍保留部分差异。当前数据不支持把C7/C9直接命名为不同的成骨或矿化亚型。外部转移/复发结果是独立扩展问题，不作为本轮发现结论的必要条件；旧外部结果原样保留。','',
'## 结果位置与来源','',
'- `C7_C9_primary_results.tsv`、`C7_vs_C9_direct_tests.tsv`：主候选结果。',
'- `discovery_program_tests.tsv`：全部260项程序检验；`within_myCAF_program_tests.tsv`：myCAF内部检验。',
'- `discovery_marker_tests.tsv`：候选基因检验，保留原始P和独立家族FDR。',
'- `primary_patient_program_differences.tsv`、`within_myCAF_patient_differences.tsv`：患者效应，可逐行复核。',
'- `within_cohort_program_tests.tsv`、`leave_one_cohort_out.tsv`、`minimum20_cell_sensitivity.tsv`：敏感性及全部阴性结果。','',
'[GSE152048原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/)、[GSE162454原始研究](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full)、[GSE318840原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13330631/)提供本轮招募来源依据；[SciPy符号秩检验说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)说明检验与零值/并列秩处理。']
(O/'delivery/发现队列主结论统计评估.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
gene=G[(G.layer=='supported')&G.cluster.isin([7,9])];summary=R.groupby(['layer','result']).size().reset_index(name='tests');summary.to_csv(O/'tables/test_result_counts.tsv',sep='\t',index=False)
qa=json.loads((O/'delivery/QA.json').read_text());qa.update(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),primary_feature_conclusion_supported=True,direct_C7_C9_myCAF_FDR=float(D[(D.population=='all_candidate_cells')&(D.program=='myCAF_customer')].FDR.iloc[0]),direct_C7_C9_iCAF_FDR=float(D[(D.population=='all_candidate_cells')&(D.program=='iCAF_customer')].FDR.iloc[0]),novel_lineage_or_metastasis_claim=False,external_validation_required=False)
(O/'delivery/QA.json').write_text(json.dumps(qa,indent=2));(O/'STATUS.json').write_text(json.dumps({'status':'DISCOVERY_FEATURE_STATISTICS_COMPLETE','QA':qa,'latest_report':'delivery/发现队列主结论统计评估.md'},ensure_ascii=False,indent=2),encoding='utf-8')
inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [B/'CAF_evidence_20260921/tables/within_donor_program_effects.tsv',B/'CAF_evidence_20260921/tables/gene_candidate_donor_effects.tsv']}
(O/'SOURCE_SHA256.json').write_text(json.dumps(inputs,indent=2))
manifest={str(p.relative_to(O)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['tables','delivery','figures','code'] for p in (O/folder).glob('*') if p.is_file() and p.name!='SHA256.json'}
(O/'delivery/SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');(O/'ANALYSIS_COMPLETE').write_text('DISCOVERY_FEATURE_STATISTICS_COMPLETE\n')
print(json.dumps(qa),flush=True);print(summary.to_string(index=False))
