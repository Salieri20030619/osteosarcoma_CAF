import pathlib,json,html,hashlib
import numpy as np,pandas as pd
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase102_NCI_author_CAF_20261020'
x=pd.read_csv(P/'tables/NCI_author_CAF_external_contrasts.tsv',sep='\t')
a=pd.read_csv(P/'tables/author_CAF_patient_expression.tsv',sep='\t')
old=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/original_tumor_patient_expression.tsv',sep='\t')
nb=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/bone_reference_with_coverage_sensitivity.tsv',sep='\t');nb=nb[nb.specification.eq('min20')]
nl=pd.read_csv(R/'phase98_reference_thresholds_20261020/tables/lung_reference_pseudobulk.tsv',sep='\t');nl=nl[nl.specification.eq('min20')]
errors=[]
for _,r in x[x.status.eq('evaluated')].iterrows():
 arrays=[a[a.stratum.eq(r.stratum)&a.gene.eq(r.gene)].expression.to_numpy(),old[old.source.eq(r.bone_source)&old.site.eq('Bone')&old.gene.eq(r.gene)].expression.to_numpy()]
 c=np.array([1.,-1.])
 if r.contrast=='normal_min20':arrays.extend([nl[nl.gene.eq(r.gene)].expression.to_numpy(),nb[nb.gene.eq(r.gene)].expression.to_numpy()]);c=np.array([1.,-1.,-1.,1.])
 y=np.concatenate(arrays);X=np.zeros((len(y),len(arrays)));s=0
 for j,z in enumerate(arrays):X[s:s+len(z),j]=1;s+=len(z)
 inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;h=np.einsum('ij,jk,ik->i',X,inv,X)
 v=inv@(X.T@(((y-X@beta)**2/(1-h))[:,None]*X))@inv
 effect=float(c@beta);se=float(np.sqrt(c@v@c))
 errors.extend([abs(effect-r.effect),abs(se-r.SE)])
assert max(errors)<1e-10
q=json.loads((P/'resources/STATISTICAL_QA.json').read_text());q.update(independent_OLS_HC2_checks=len(errors)//2,OLS_HC2_max_error=max(errors))
(P/'resources/STATISTICAL_QA.json').write_text(json.dumps(q,indent=2),encoding='utf-8')
def table(d):return d.to_html(index=False,escape=True,border=0,na_rep='不可评估',float_format=lambda v:format(v,'.5g'))
main=x[x.stratum.eq('author_all')].copy()
rows=[]
for gene,z in main.groupby('gene',sort=False):
 d={'基因':gene}
 for source,label in [('SCPCP000017','骨队列17'),('SCPCP000023','骨队列23')]:
  for contrast,lab in [('tumor_only','肿瘤差异'),('normal_min20','背景扣除')]:
   r=z[z.bone_source.eq(source)&z.contrast.eq(contrast)].iloc[0]
   d[label+' '+lab+' 效应']=r.effect;d[label+' '+lab+' P']=r.P;d[label+' '+lab+' FDR']=r.FDR
 rows.append(d)
coverage=a[a.gene.eq('EBF1')].pivot(index='patient_id',columns='stratum',values='cells').reset_index()
summ=json.loads((P/'resources/RESULT_SUMMARY.json').read_text())
intro='''<h1>NCI 作者 CAF 分析与 NCC 补充材料核查</h1><p>内部分析版本 · 2026-10-20 · phase102</p>
<div class="lead"><b>作者标注 CAF 为 1,530 个，覆盖 17 位患者。外部比较中，7 个基因在两套骨参考队列下均与原方向一致且名义显著；正常背景扣除及单细胞敏感性支持其中一部分信号。</b></div>
<p>作者 CAF 全集按 Ann_Level0=Host 且 Ann_Level2=CAFs 提取，不使用此前的 CopyKAT、参考细胞验证或成纤维评分门槛。此前 478 个是作者标签与扩展筛选的交集，本次使用全部 1,530 个。19 个文库按用户确认合并为 17 位患者：16178/16178b、16183/16183b 分别合并；计数先相加再归一化。</p>
<h2>外部比较保留 7 个同向名义显著基因</h2>
<p>NCI 肺 CAF（17 位患者）分别与原 SCPCP000017 骨样本（12 位）和 SCPCP000023 骨样本（4 位）比较。NEBL、PNKD、TOR1AIP1、GRID2、BEST1、RGS2-AS1、AL096854.1 在两项比较中均同向且 P&lt;0.05。其中 TOR1AIP1、GRID2、BEST1、RGS2-AS1、AL096854.1 在各自 14 基因检验家族内均 FDR&lt;0.05。两个比较共享同一批 NCI 患者，不能视为两次独立外部验证。</p>
<p>RGS2-AS1 的肺减骨效应分别为 1.864（P=3.11×10⁻⁵，FDR=0.000436）和 2.474（P=0.00171，FDR=0.00478）。效应单位为 log2(CP10k+0.1) 的组均值差，并非未经伪计数修正的倍数变化。</p>
<h2>正常背景扣除后，5 个基因同向且通过候选家族 FDR</h2>
<p>固定正常参考每位至少 20 个细胞：骨 5 位、肺 3 位。12 个可评估基因中，TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1 在两套骨参考比较下均同向且 FDR&lt;0.05。PEF1-AS1 与 AL096854.1 因正常参考缺失，不能进行背景扣除。</p>
<p>排除预测双细胞后保留 1,320 个作者 CAF，仍为 17 位患者。肿瘤比较的 7 个同向名义显著基因保持不变；背景扣除后，TOR1AIP1、GRID2、TAF10、RGS2-AS1 仍在两套比较中同向且 FDR&lt;0.05，BEST1 不再同时满足。因此 TOR1AIP1、GRID2、RGS2-AS1 同时获得未扣除背景及背景扣除结果的支持，且通过单细胞子集检查；TAF10 的双队列支持主要来自背景扣除分析。</p>
<h2>新增数据不能支持全部 14 个基因一并复现</h2>
<p>RTKN、ZNF526、CALHM5、PEF1-AS1、EBF1 在全部作者 CAF 的两套比较中未同时保持原方向；NDUFS7 与 TAF10 虽在肿瘤比较中保持同向，但未同时达到名义显著。EBF1 在背景扣除后两套比较均变为肺高，与原骨高方向相反。这些结果保留在完整统计表中。</p>
<h2>来源调整后的原骨—肺效应没有获得新增组内信息</h2>
<p>控制研究来源的共同部位效应模型中，NCI 仅有肺样本，其来源内部不存在骨—肺对照，故新增 NCI 对部位系数的信息权重为零。数值核查确认加入前后系数相同。本报告没有把这种“系数不变”作为新增验证成功。原来源内部精确置换检验的 14 个候选结果保持原版本，本次另列共同效应模型与外部比较。</p>
<h2>方法</h2>
<p>每位患者汇总作者 CAF 的原始计数，以全部可用基因计数为分母，计算 log2(CP10k+0.1)。主要细胞集合为作者 CAF 全集；作者 CAF 中 scDblFinder=singlet 为敏感性集合。NCI 患者不另设 20 细胞门槛；20 细胞仅用于既定正常参考筛选。以患者为独立单位，NCI 肺减原骨用 Welch 两组检验；背景扣除用四组独立均值线性对比及 Welch–Satterthwaite 自由度。BH 按骨队列×细胞集合×对比类型分别计算；肿瘤比较 14 个有效检验，正常背景对比 12 个。提供原始 P、FDR、效应、95%CI、SE 和自由度。</p>
<p>来源调整共同效应以每个来源的 n骨×n肺/(n骨+n肺) 为权重，计算来源内部肺减骨对比，方差允许组间异方差；NCI 权重为零。普通最小二乘独立核查部位系数，组别饱和模型 HC2 独立核查两组与四组对比的效应和标准误。</p>
<h2>讨论</h2>
<p>NCI 为单核测序，且本次使用作者 CAF 标签；原骨比较保留原扩展候选定义。外部比较同时包含部位、来源、测序方式及注释口径差异，属于探索性支持，不能独立归因于骨—肺部位。正常参考来自其他队列，背景扣除不等同于充分消除来源效应。SCPCP000023 骨组仅 4 位，相关比较按小样本探索性解释。候选家族 FDR 与原全转录组 FDR 不同，不可直接互换。作者合并对象包含 14 个目标特征，但原始逐文库测量掩码不可得。</p>
<h2>NCC 补充材料与下载</h2>
<p>上次为服务器直连且总时限 180 秒。本次通过已配置 sgs10-codex 转发使用本机 127.0.0.1:7897 代理，允许 1,200 秒并断点续传。补充文件 1 已完整下载：21,165,445 字节，MD5=734d6b81da4646c2ecd11ced4071f554，与出版社一致。完整渲染核查显示其为单页细胞注释、比例及表达图，包含 BS 样本标识，未提供逐样本实际取材部位映射。</p>'''
fullcols=['gene','stratum','bone_source','contrast','NCI_N','Bone_N','effect','SE','df','CI_low','CI_high','P','FDR','family_valid_tests','same_direction','status']
body=intro+'<h2>主分析完整效应与 P/FDR</h2><div class="wide">'+table(pd.DataFrame(rows))+'</div><h2>每位患者作者 CAF 数</h2>'+table(coverage)+'<details><summary>展开全部主分析及单细胞敏感性统计</summary>'+table(x[fullcols])+'</details><details><summary>展开来源调整对比</summary>'+table(pd.read_csv(P/'tables/source_adjusted_site_contrasts.tsv',sep='\t'))+'</details><p>统计核查：'+html.escape(json.dumps(q,ensure_ascii=False))+'</p>'
doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>NCI 作者 CAF 分析</title><style>body{font:16px/1.75 "Microsoft YaHei",sans-serif;color:#20303c;max-width:1180px;margin:36px auto;padding:0 26px}h1{font-size:28px}h2{font-size:21px;margin-top:32px}.lead{padding:18px;background:#e8f3f3;border-left:5px solid #237c83}table{border-collapse:collapse;font-size:13px;width:100%;margin:14px 0}th,td{border:1px solid #d6e0e4;padding:7px;text-align:left}th{background:#edf3f6}.wide,details{overflow-x:auto}summary{cursor:pointer;font-weight:bold}</style><body>'+body+'</body></html>'
(P/'report/NCI_author_CAF_analysis_20261020.html').write_text(doc,encoding='utf-8')
print('REPORT',len(doc),q)
