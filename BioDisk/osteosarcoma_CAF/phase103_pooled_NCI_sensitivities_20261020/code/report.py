import pathlib,json,base64,html,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase103_pooled_NCI_sensitivities_20261020'
d=pd.read_csv(P/'tables/all_14gene_results.tsv',sep='\t')
s=pd.read_csv(P/'tables/scenario_summary.tsv',sep='\t').fillna('')
coverage=pd.read_csv(P/'tables/cohort_coverage.tsv',sep='\t')
genes=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text());order=[g['gene'] for g in genes]
scenarios=[
('expanded','with_NCI','pooled_Welch','Expanded + NCI'),
('expanded','without_NCI','pooled_Welch','Expanded - NCI'),
('expanded','with_NCI','source_adjusted_OLS_HC3','Expanded / source'),
('singlets','with_NCI','pooled_Welch','Singlets + NCI'),
('conservative','with_NCI','pooled_Welch','Conservative + NCI'),
('conservative','without_NCI','pooled_Welch','Conservative - NCI'),
('conservative','with_NCI','source_adjusted_OLS_HC3','Conservative / source'),
('expanded_matched_strict_donors','with_NCI','source_adjusted_OLS_HC3','Matched / source')]
fig,axes=plt.subplots(1,2,figsize=(18,9.5),constrained_layout=True)
for ax,ct,title in zip(axes,['tumor_only','normal_min20'],['Tumor Lung minus Bone','After normal-background subtraction']):
 vals=np.zeros((14,len(scenarios)));rr={}
 for j,(cell,co,mo,lab) in enumerate(scenarios):
  z=d[d.cellset.eq(cell)&d.cohort.eq(co)&d.model.eq(mo)&d.contrast.eq(ct)].set_index('gene').loc[order]
  vals[:,j]=z.effect.to_numpy()
  for i,g in enumerate(order):rr[i,j]=z.loc[g]
 im=ax.imshow(np.ma.masked_invalid(vals),cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-4,vcenter=0,vmax=4),aspect='auto')
 ax.set_xticks(range(len(scenarios)),[v[3] for v in scenarios],rotation=45,ha='right',fontsize=9)
 ax.set_yticks(range(14),order);ax.set_title(title,fontsize=13)
 for (i,j),r in rr.items():
  text='NA' if not np.isfinite(r.effect) else f'{r.effect:.2f}'+('**' if r.FDR<.05 else '*' if r.P<.05 else '')
  ax.text(j,i,text,ha='center',va='center',fontsize=8,color='white' if np.isfinite(r.effect) and abs(r.effect)>2.5 else 'black')
 ax.set_xticks(np.arange(-.5,len(scenarios),1),minor=True);ax.set_yticks(np.arange(-.5,14,1),minor=True);ax.grid(which='minor',color='white',linewidth=.5);ax.tick_params(which='minor',bottom=False,left=False)
fig.colorbar(im,ax=axes,label='Effect on log2(CP10k + 0.1) scale; color clipped at +/-4',shrink=.65)
fig.suptitle('Patient-level pooled and sensitivity analyses | * P < 0.05; ** candidate-family FDR < 0.05',fontsize=14)
fig.savefig(P/'figures/scenario_effects.png',dpi=180);plt.close(fig)
def table(df):return df.to_html(index=False,border=0,na_rep='不可评估',float_format=lambda x:format(x,'.5g'),escape=True)
main=d[d.cellset.eq('expanded')&d.cohort.eq('with_NCI')&d.model.eq('pooled_Welch')]
m=main[main.contrast.eq('tumor_only')][['gene','effect','CI_low','CI_high','P','FDR','same_direction']]
names={'expanded':'扩展候选','singlets':'排除预测双细胞','conservative':'成纤维保守子集','expanded_matched_strict_donors':'扩展候选／保守患者名单','with_NCI':'加入NCI','without_NCI':'去掉NCI','pooled_Welch':'合并Welch','source_adjusted_OLS_HC3':'来源调整HC3（有效来源df）','source_adjusted_OLS_HC3_global_df':'来源调整HC3（全局df诊断）','tumor_only':'肿瘤部位差异','normal_min20':'正常背景扣除'}
core=s[~s.model.eq('source_adjusted_OLS_HC3_global_df')].copy()
for col in ['cellset','cohort','model','contrast']:core[col]=core[col].map(names)
core=core[['cellset','cohort','model','contrast','Bone_N','Lung_N','same_direction_nominal_n','same_direction_FDR_n','same_direction_nominal_genes','opposite_nominal_genes']]
core.columns=['细胞集合','是否加入NCI','模型','对比','骨患者','肺患者','同向P<0.05数','同向FDR<0.05数','同向名义显著基因','反向名义显著基因']
parts=['''<h1>骨—肺联合分析：加入 NCI 后的结果与敏感性检查</h1><p>内部完整分析 · 2026-10-20 · phase103</p>
<div class="lead"><b>骨16例、肺23例的扩充比较中，14个候选全部保持原效应方向，9个同向且候选家族FDR&lt;0.05。排除预测双细胞后这9个不变。严格成纤维子集与来源调整联合使用时支持明显减弱，当前不能称为跨全部方案稳健。</b></div>
<h2>扩充比较支持9个同向候选，正常背景扣除后保留6个</h2>
<p>联合比较的9个基因为 NEBL、NDUFS7、PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1、AL096854.1，均为肺高。RGS2-AS1 效应为2.104，P=1.21×10⁻⁷，FDR=1.70×10⁻⁶；PNKD 为1.074，P=0.000256，FDR=0.000911；GRID2 为1.716，P=0.000260，FDR=0.000911。效应单位为患者 log2(CP10k+0.1) 表达均值差，不是未经伪计数处理的直接倍数变化。</p>
<p>保留每位至少20个正常参考细胞后，固定使用正常骨5位、正常肺3位。在12个可评估基因中，PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1 保持原方向且FDR&lt;0.05。EBF1 在扣除背景后反转为肺高（效应4.634，P=0.000331，FDR=0.000795）；它不是原骨高信号的复现。PEF1-AS1、AL096854.1 缺少正常参考特征，未进行背景扣除检验。</p>
<h2>去掉NCI后仍有9个同向名义信号，但候选集合发生变化</h2>
<p>原两队列合并的骨16例、肺6例中，RTKN、NDUFS7、PNKD、TOR1AIP1、GRID2、BEST1、RGS2-AS1、AL096854.1、EBF1 同向且P&lt;0.05，其中5个FDR&lt;0.05。加入NCI后，NEBL和TAF10进入同向名义显著集合，RTKN和EBF1退出；两种方案共有7个同向名义信号。这显示补充肺患者后结果发生实际变化，不能只凭显著基因数量判断一致性。</p>
<h2>双细胞排除稳定，严格身份筛选使结果减弱</h2>
<p>排除预测双细胞后，原主队列纳入细胞未变化；NCI由1530降为1320，患者仍为17位。扩充比较仍为骨16例、肺23例，上述9个候选均保持同向FDR显著；正常背景扣除后同向FDR显著仍为6个。</p>
<p>严格集合采用原冻结 conservative_fibro 且排除预测双细胞，NCI采用作者CAF单细胞子集。实际为骨13例、肺22例，共2471个细胞。未调整来源时，NEBL、NDUFS7、PNKD、GRID2、BEST1、RGS2-AS1、AL096854.1 共7个同向名义显著，其中除NDUFS7外6个通过FDR。去掉NCI后仅剩骨13例、肺5例，未检出同向名义显著基因。</p>
<p>严格集合再调整来源后，肿瘤部位比较仅 NDUFS7 名义显著（效应1.393，P=0.0187，FDR=0.262），没有FDR显著基因。正常背景扣除后 TOR1AIP1 与 TAF10 同向名义显著，只有 TAF10 通过FDR（效应2.333，P=0.00159，FDR=0.0191）。EBF1在该背景对比仍是反向名义显著。</p>
<h2>来源调整结果不能归功于新增NCI的部位信息</h2>
<p>扩展候选的来源调整模型中，12个基因同向且FDR&lt;0.05（NEBL、CALHM5未通过）。这一结果在去掉NCI后同样成立：NCI只有肺，其来源内部没有骨—肺变化，对来源调整后部位系数的信息权重为零。因此调整模型支持来自原研究内部比较，不能称为NCI独立验证。</p>
<p>初次尝试采用模型全局残差自由度，这会让NCI增加自由度、在部位效应和标准误不变时改变P值。完整保留该尝试作为诊断；主来源调整结果改为仅使用具有骨、肺两组的来源计算推断自由度。该约定下加入前后调整效应与P值一致。此前全局自由度下严格子集RTKN达到名义显著；有效来源自由度下P=0.0551，未达到阈值，故没有将其计入主结论。</p>
<h2>患者名单变化与细胞定义变化都影响结果</h2>
<p>将扩展表达限制在严格集合保留的患者名单，仍用骨13例、肺22例，未调整比较有8个同向名义信号。来源调整后有7个同向名义信号，但均未通过FDR；进一步换成严格细胞表达后只剩NDUFS7名义显著。患者减少并不能完全解释严格筛选后的减弱，细胞集合定义也影响效应及不确定性。</p>
<h2>实际纳入6633个细胞，不是全部9704个</h2>
<p>8174是跨多个队列的完整冻结候选池。本次按已约定的两个原主队列及非配对名单，使用原细胞5103个，加NCI作者CAF1530个，共6633个。8174中另有2446个来自其他队列、162个为这两个主队列的未知/混合部位、463个来自原方案排除的配对患者，本次均未纳入。没有把细胞数量当作患者数。</p>
<p>NCI按用户确认合并16178/16178b、16183/16183b，同患者文库先汇总原始计数再归一化。公开患者标识去除来源前缀后未发现跨队列重复；没有全体患者基因型或跨研究作者对应表，未将这一核查写成遗传身份确认。严格子集中SCPCP000023仅剩2位骨患者、共7个细胞，完整保留其覆盖度和小样本解释。</p>
<h2>方法</h2>
<p>目标为预先指定14个候选基因；原两主队列SCPCP000017、SCPCP000023的患者与独立设计名单沿用phase81，配对患者保持排除。每位患者、每个部位汇总原始计数，以完整可用基因UMI总量归一化后取log2(CP10k+0.1)。原扩展患者表达308个值与冻结版本逐项一致。NCI使用phase102作者CAF全集或其singlet子集。无目标基因缺失；原逐文库特征覆盖掩码用于原队列有效性检查。</p>
<p>四种细胞集合分别为扩展候选、排除预测双细胞、成纤维保守集合、按保守集合患者名单匹配的扩展表达。每种集合分别加入或去掉NCI，计算非配对Welch比较及来源固定效应OLS（截距、肺指示变量、来源哑变量）目标肺系数，HC3异方差稳健协方差。主要来源调整推断使用具有部位组内变化的来源的残差自由度；全局残差自由度结果另存诊断。</p>
<p>正常背景对比为肿瘤肺减骨对比减去正常肺减骨均值差。固定正常骨5位、肺3位，每位至少20个参考细胞；参考方差/n加至肿瘤对比方差，并用Welch–Satterthwaite自由度计算双侧P及95%CI。20细胞阈值不施加于肿瘤患者。P&lt;0.05为名义显著，FDR&lt;0.05为候选家族FDR显著。BH分别在细胞集合×队列组合×模型×背景对比的同类14候选内计算；正常背景12个有效检验。原全转录组FDR保持原版，不与本次候选家族FDR混用。</p>
<h2>讨论</h2>
<p>新增肺患者使非配对联合分析可行，并为若干候选提供扩充样本证据。NCI为单核测序且采用作者CAF标签，原扩展集合含身份未决细胞；研究来源、测序方式和细胞注释差异需要与部位效应一并解释。严格子集及患者匹配对照显示身份筛选与患者覆盖均影响结果，当前数据不支持把14基因整体描述为跨全部方案稳健的CAF部位标志。</p>
<p>本次为既有候选的探索性再分析，不是新的全转录组发现或独立验证。正常背景扣除服务于组织背景解释，不作为剔除骨—肺区分候选的唯一门槛。模型、子集和背景分析分别报告，未通过挑选某一最小P值合成“稳健”结论。原来源内精确置换结果与本次合并Welch/HC3检验不同，不能把检验方法改变造成的差异全部归因于新增样本。</p>
<h2>方案效应图</h2><p>红色为肺高，蓝色为骨高；格内为效应。*：原始P&lt;0.05；**：本方案候选家族FDR&lt;0.05。颜色限制于±4，数值未截断；NA表示正常参考不可评估。“Matched”为扩展表达使用保守患者名单。</p>''']
img=base64.b64encode((P/'figures/scenario_effects.png').read_bytes()).decode()
parts+=['<img style="width:100%" alt="各方案效应与显著性" src="data:image/png;base64,'+img+'">','<h2>扩充主比较14基因完整统计</h2>'+table(m),
'<h2>32项主要方案汇总</h2><div class="scroll">'+table(core)+'</div>',
'<h2>患者与细胞覆盖</h2>'+table(coverage),
'<details><summary>展开全部672行统计（含全局自由度诊断）</summary><div class="scroll">'+table(d)+'</div></details>',
'<details><summary>核查记录与方案变更</summary><pre>'+html.escape((P/'resources/STATISTICAL_QA.json').read_text()+'\n'+(P/'resources/PROTOCOL_SOURCE_DF_AMENDMENT.json').read_text())+'</pre></details>']
doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CAF扩充联合分析</title><style>body{font:16px/1.8 "Microsoft YaHei",sans-serif;color:#20303c;max-width:1240px;margin:36px auto;padding:0 25px}h1{font-size:29px}h2{font-size:22px;margin-top:32px}.lead{background:#e8f3f3;border-left:5px solid #237c83;padding:18px}table{border-collapse:collapse;font-size:12px;width:100%;margin:15px 0}th,td{border:1px solid #d6e0e4;padding:7px;text-align:left}th{background:#edf3f6}.scroll,details{overflow-x:auto}summary{cursor:pointer;font-weight:bold}pre{white-space:pre-wrap;font-size:12px}</style><body>'+''.join(parts)+'</body></html>'
f=P/'report/CAF_pooled_NCI_sensitivities_20261020.html';f.write_text(doc,encoding='utf-8')
print('REPORT',len(doc),hashlib.sha256(f.read_bytes()).hexdigest())
