from pathlib import Path
import json,pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase60_F3_review_20261005'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def tab(d):return d.to_csv(sep='\t',index=False,float_format='%.6g',na_rep='NA')
assert (O/'QA_COMPLETE').exists();c=read(O/'tables/F3_cells_review_QC.tsv.gz');s=read(O/'tables/program_statistics.tsv');e=read(O/'tables/patient_program_effects.tsv');sens=read(O/'tables/F3_threshold_sensitivity_summary.tsv');bg=read(O/'tables/background_QC.tsv');roster=read(O/'tables/patient_site_roster.tsv');qa=json.loads((O/'delivery/QA_ACCEPTANCE.json').read_text());main='SCPCP000017';m=c[c.dataset==main];donors=['SCPCP000017|17-439_00067','SCPCP000017|CCG1106_010','SCPCP000017|CCG1106_077'];short={d:d.split('|')[-1] for d in c.donor_record.unique()};short['GSE318840|GSE318840__Osteosarcoma_Metastasis_1']='GSE318840_Met1';pal=['#E07A5F','#3D85A7','#63A375','#B595CE','#D6AC4D'];colors=dict(zip(sorted(c.donor_record.unique()),pal));plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,aa=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
for ax,col in zip(aa[:2],['donor_record','lineage_stratum']):
 for j,(label,g) in enumerate(c.groupby(col)):
  ax.scatter(g.UMAP1,g.UMAP2,s=20,alpha=.8,label=short.get(label,label),color=colors[label] if col=='donor_record' else pal[j])
 ax.set_title('Frozen F3 UMAP by '+('patient' if col=='donor_record' else 'expression lineage'),weight='bold');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(fontsize=8,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2)
ct=pd.crosstab(c.round1_state,c.donor_record);bottom=np.zeros(len(ct))
for d in ct.columns:aa[2].bar(range(len(ct)),ct[d],bottom=bottom,label=short[d],color=colors[d]);bottom+=ct[d]
aa[2].set_xticks(range(3),['C0','C1','C2']);aa[2].set_ylabel('Cells');aa[2].set_title('Internal states are patient concentrated',weight='bold');fig.savefig(O/'figures/01_F3_composition.png',dpi=180);plt.close(fig)
fig,aa=plt.subplots(1,2,figsize=(12,4.8),layout='constrained');progs=['MYC_TARGETS_V1','OXIDATIVE_PHOSPHORYLATION','TRANSLATION','MHCII','cycling_markers']
for ax,ctx in zip(aa,['F3_vs_rest','F3_Fibro_vs_rest_Fibro']):
 d=e[(e.dataset==main)&(e.matrix=='corrected')&(e.context==ctx)&e.program.isin(progs)].pivot(index='donor_record',columns='program',values='effect').reindex(index=donors,columns=progs);im=ax.imshow(d,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto');ax.set_yticks(range(3),[short[d] for d in donors]);ax.set_xticks(range(len(progs)),['MYC V1','OXPHOS','Translation','MHC-II','Cycle markers'],rotation=30,ha='right');ax.set_title(ctx.replace('_',' '),weight='bold')
 for i in range(3):
  for j in range(len(progs)):ax.text(j,i,f'{d.iloc[i,j]:+.2f}',ha='center',va='center',color='white' if abs(d.iloc[i,j])>1.8 else 'black')
fig.colorbar(im,ax=aa,label='Within-patient mean gene log2(CP10k+0.1) effect',shrink=.8);fig.supxlabel('Corrected counts. Three patients; effects are exploratory, not independent cell-level evidence.',fontsize=9);fig.savefig(O/'figures/02_patient_program_effects.png',dpi=180);plt.close(fig)
fig,aa=plt.subplots(1,3,figsize=(12,4.5),layout='constrained');panels=[('MYC_TARGETS_V1','MYC V1'),('OXIDATIVE_PHOSPHORYLATION','OXPHOS'),('MHCII','MHC-II')];contexts=['F3_Fibro_vs_rest_Fibro','F3_Fibro_vs_BG_Fibro','F3_Fibro_vs_BG_Fibro_UMI_matched'];labels=['Other F Fibro','Broad Fibro','UMI-matched Fibro']
for ax,(prog,title) in zip(aa,panels):
 for d in donors:
  vals=[e[(e.dataset==main)&(e.matrix=='raw')&(e.context==ctx)&(e.program==prog)&(e.donor_record==d)].effect.iloc[0] for ctx in contexts];ax.plot(range(3),vals,'o-',label=short[d],color=colors[d],lw=1)
 ax.axhline(0,color='#999',lw=.8);ax.set_xticks(range(3),labels,rotation=30,ha='right');ax.set_title(title,weight='bold');ax.set_ylabel('Raw-count within-patient effect')
aa[-1].legend(frameon=False,fontsize=8);fig.supxlabel('Reference sensitivity: all comparisons match source/patient/site. UMI matching is a diagnostic, not a causal adjustment.',fontsize=9);fig.savefig(O/'figures/03_reference_sensitivity.png',dpi=180);plt.close(fig)
fig,aa=plt.subplots(1,3,figsize=(12,4.5),layout='constrained');pm=m.groupby('donor_record').agg(UMI=('UMI','median'),features=('detected_features','median'),MT=('mitochondrial_fraction','median'));b=bg[(bg.dataset==main)&(bg.background=='Fibro_UMI_matched')].set_index('donor_record')
for ax,metric,bc in zip(aa,['UMI','features','MT'],['median_UMI','median_features','median_MT']):
 for i,d in enumerate(donors):ax.plot([0,1],[pm.loc[d,metric],b.loc[d,bc]],'o-',label=short[d],color=colors[d])
 ax.set_xticks([0,1],['F3 all','Matched BG-Fibro']);ax.set_title('Median '+metric,weight='bold')
 if metric=='MT':ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
aa[-1].legend(frameon=False,fontsize=8);fig.supxlabel('Matching uses the 10th-90th raw-library percentiles of same-sample F3-Fibro. Values describe cells; no cell-level P values.',fontsize=9);fig.savefig(O/'figures/04_QC_context.png',dpi=180);plt.close(fig)
fig,aa=plt.subplots(1,2,figsize=(11,4.5),layout='constrained');rules=[('raw','at_least_1','Raw >=1'),('corrected','nonzero','Corrected >0'),('corrected','at_least_0.5','Corrected >=0.5'),('corrected','at_least_1','Corrected >=1')]
for ax,key,title in [(aa[0],'MHCII_fibro_co','MHC-II + fibro co-detection'),(aa[1],'cycle_ge2','At least two cycling markers')]:
 for i,d in enumerate(donors):
  vals=[sens[(sens.donor_record==d)&(sens.matrix==mode)&(sens.threshold_rule==rule)][key].sum() for mode,rule,l in rules];ax.plot(range(4),vals,'o-',label=short[d],color=colors[d])
 ax.set_xticks(range(4),[x[2] for x in rules],rotation=25,ha='right');ax.set_title(title,weight='bold');ax.set_ylabel('Number of cells');ax.legend(frameon=False,fontsize=8)
fig.supxlabel('Fractional corrected counts: a raw count of 1 can become 0.98 without meaningful signal removal. All thresholds retained.',fontsize=9);fig.savefig(O/'figures/05_detection_thresholds.png',dpi=180);plt.close(fig)
summary=json.loads((O/'delivery/SUMMARY.json').read_text());summary['QC']=bg.to_dict('records');summary['threshold_sensitivity']=sens.to_dict('records');(O/'delivery/SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
qcf=m.groupby(['donor_record','site']).agg(cells=('cell_id','size'),median_UMI=('UMI','median'),median_features=('detected_features','median'),median_MT=('mitochondrial_fraction','median'),median_removed_fraction=('fraction_removed_from_library','median')).reset_index();pview=s[(s.dataset==main)&(s.matrix=='corrected')&(s.context=='F3_vs_rest')];threshold=sens.groupby(['matrix','threshold_rule'])[['cells','MHCII_ge2','MHCII_fibro_co','cycle_ge2','PTPRC_detected']].sum().reset_index()
text='''CAF项目 F4表达分层与F3深入核查
2026-10-05  phase60

F4可以按两种主要表达谱系分层推进。F3目前应保留为患者集中、身份混合的状态标签，当前数据不支持将其认定为稳定的独立CAF功能亚型，也不足以把整群命名为apCAF。

本轮F3的关键发现不仅是样本量少，还在于原始F3对其余四群的参照细胞极少：两位贡献主要F3细胞的患者，其余四群在匹配部位仅有2个和3个细胞。扩大到同患者Fibro背景后，原先的MYC/氧化磷酸化优势明显减弱；按相近RNA量选取Fibro背景后，三个患者均不再表现该优势。因此F3的高程序排序不能直接解释为独立细胞类型。

F4分层的实施结果

已导出F4_Fibro_like 909个、F4_Osteo_like 3300个、F4_other_unresolved 1095个，三者合计5304个原始F4细胞。Osteo_like包括原Osteo 3297个和Chondro 3个；残余层保留47个Mural/SMC表达细胞及1048个其他/未定细胞。原始F4标签和细胞索引不变。

这两种主表达层的差异已有phase59支持，可以分别做后续比较；它们并非本轮重新聚类得到的两个稳定亚群。“Fibro_like”不等于已经确证CAF，“Osteo_like”不等于已经病理确诊恶性细胞。未将1095个残余细胞强行并入二类。

方法

沿用原始8174对象中的251个F3、冻结的C0/C1/C2状态和UMAP，叠加既有表达复核、RNA-CNV、双细胞及质量记录。使用原始和DecontX计数、完整library分母及既有canonical基因映射。分析32个预设程序/诊断面板，并对MYC、氧化磷酸化、翻译集合增加去除MT-/RPL/RPS基因的敏感性。具体集合及版本保存在resources/GENE_SETS.json。

患者级表达效应为同来源同患者同部位pseudobulk的log2(CP10k+0.1)差值，在集合内平均、部位内配对后对患者等权汇总。每项至少覆盖80%定义基因，且至少min(3,集合大小)个被测基因。缺测保留NA，不补生物学零。比较包括F3对其余四群、Fibro/Osteo限定比较、冻结状态对F3其余状态，以及同患者原始8174之外的Fibro/Osteo/Myeloid背景。另用同样本F3-Fibro library的10%–90%范围选择相近RNA量的Fibro背景。

使用患者级双侧精确符号翻转检验，N<2或全零效应不输出P。按来源、原始/去污染版本、比较情境及完整程序/诊断面板/基因剔除敏感性分别定义FDR家族；标志基因另设家族。原始P<0.05为名义候选标准。细胞共检出只作描述，基线要求每套面板至少2个基因且各计数≥1；同时报告校正计数>0及≥0.5的敏感性，以处理分数计数的阈值效应。所有原始细胞及低覆盖记录保留。

F3主要由两位患者构成，三个内部状态缺少跨患者支撑

全体F3来自3个来源、5个来源内患者记录。SCPCP000017贡献245/251个；其余两个来源仅1个和5个，不能作为可靠的外部复制。
主来源的患者分别贡献124、118和3个F3细胞。前两位合计242/251（96.4%）。C0为119个细胞，其中115个来自CCG1106_077右肺样本；C1为42个，其中40个来自17-439_00067椎旁样本；C2为90个，其中84个来自同一椎旁样本。C0与C1/C2主要反映患者/取材差异。C1与C2可以描述为同一患者内的表达状态差别，不能据此建立跨患者稳定亚型。

主来源逐患者组成与质量
'''+tab(qcf)+'''
F3对其余四群的对照规模不足以支撑强分型结论

主来源实际配对分别为124对2、118对3、3对162个细胞；双方都要求至少10细胞时，没有可评估患者对。限定Fibro后，两位主要患者的参照更分别只有1个细胞。这里不是把细胞数当作统计N，而是说明每位患者表达估计的稳定性受其实际细胞覆盖影响。

实际配对记录
'''+tab(roster[(roster.dataset==main)&(roster.matrix=='raw')&roster.context.isin(['F3_vs_rest','F3_Fibro_vs_rest_Fibro'])])+'''
MYC 翻译和能量信号存在方向，但优势依赖参照

原始F群比较中，DecontX后MYC V1、氧化磷酸化、翻译平均效应分别为+1.017、+0.847、+1.072，均为2/3患者同向；两个主要患者正向，只有3个F3细胞的患者反向。Fibro限定后仍有相同方向结构。去掉MT-/RPL/RPS基因后，MYC为+0.960、翻译为+0.743，说明相对旧参照的信号并非完全由这些基因组成造成。冻结氧化磷酸化集合在本版本与该剔除集合相同，应视作无变化敏感性，不能当作独立验证。

与同患者广义Fibro背景比较，F3-Fibro的MYC/OXPHOS平均效应变为−0.180/−0.300；与相近RNA量Fibro背景比较分别为−0.419/−0.466，三位患者均为负向。由于N=3，这些方向也不作确定的降低结论，但足以说明原F群内部的排序优势不具有参照独立性。

MHC-II是稀疏局部信号，不能把整群叫做apCAF

F3内共有166个Fibro表达细胞、69个Osteo/Chondro表达细胞、16个其他/未定细胞。MHC-II面板在原F群比较和Fibro限定比较中为3/3正向，但相近RNA量Fibro背景下三位患者均未显示升高。

原始计数中，仅12/251（4.8%）细胞同时检出至少2个MHC-II和2个Fibro面板基因；其中Fibro表达层7个、Osteo层4个、未定层1个。去污染后按>0或≥0.5计数仍为12个，按≥1仅为1个。这主要反映单UMI变为略小于1的分数计数，不能把12降到1直接解释成环境RNA已被去除。

原始HLA-DRA只有13个细胞检出、HLA-DRB1为4个，CD74为55个；多数不是成套HLA-II共表达。PTPRC仅1个细胞检出，当前没有证据将整群解释为明显的白细胞混入；但少量MHC-II/Fibro共表达也不足以证明抗原呈递功能。该部分应保留为局部探索线索，不做apCAF整群命名。

局部增殖信号保留，同时存在明显质量差异

原始计数57/251细胞检出至少2个增殖面板基因；去污染后>0及≥0.5口径仍为57个，≥1口径为17个。57个中，Fibro层25、Osteo层30、未定层2，主要集中在两位大贡献患者。可以记录局部增殖表达状态，不能把F3整群等同于增殖CAF。

两位主要患者F3的线粒体比例中位数为17.6%和16.2%；相近RNA量Fibro背景约为7.4%和3.1%，显示技术质量/细胞状态差异应纳入解释。其DecontX library去除比例中位数仅约0.8%和2.3%，核心数据并没有在去污染后整体消失。GSE318840的5个外围细胞则去除中位数约96.6%，不宜承担外部验证结论。另一来源只有1个细胞，并有既有双细胞阳性记录。

完整主来源程序与诊断面板结果（去污染后）
'''+tab(pview[['program','N','mean_delta','positive_patients','P','FDR']])+'''
阈值敏感性汇总
'''+tab(threshold)+'''
讨论与研究定位

本轮定位为：F3是一个由少数患者主导、含Fibro与Osteo表达成分、具有局部增殖和稀疏MHC-II线索的历史状态群。MYC/能量相关方向不能直接升级为稳定功能亚型，其对照稀疏、参照依赖及线粒体比例差异均会影响解释。相近RNA量匹配改变参照组成，属于敏感性分析，不能证明全部差异都是技术伪影。

本轮患者N最大为3，双侧精确检验最小P=0.25，所以没有原始P<0.05的患者级候选是数学约束，不等于没有表达方向。本轮将方向与效应保留，并完整输出P/FDR；没有用细胞数扩充N，也没有把既有GSEA模型P替代患者级证据。多个来源患者不能假定全局独立后合并。

RNA谱系复核和功能程序来自相关表达信息，不能算作独立实验验证。既有CNV筛选条件下未检出异常不能证明CAF身份；MHC-II表达不等于真实抗原呈递，增殖基因不等于增殖速率实测。不能仅凭本轮结果把全部F3删除为低质量细胞，也不能将其整体重命名为myCAF、mCAF、apCAF或成骨CAF。

后续数据使用口径

F4采用两个主表达层加残余层，保持索引可追溯。F3暂保留原标签；探索时限定Fibro表达层166个，报告每患者构成，使用同患者覆盖充分的广义Fibro参照。局部增殖/MHC-II作为状态维度记录，现阶段不追加新的聚类标签。空间壳层、功能实验和治疗应答不在本轮实施范围。

验收与交付

独立R复核634项有效程序P、1014项标志基因P，最大误差0；FDR最大误差1.11e-16。与phase58的22项Hallmark原始/去污染结果逐项一致。所有F3/F4原标签及数量均核对通过。

tables/F4_two_main_strata_with_residual.tsv.gz：5304个F4的两个主表达层及残余索引。
tables/F3_cells_review_QC.tsv.gz：251个F3的原标签、状态、UMAP、谱系、QC与污染信息。
tables/program_statistics.tsv、marker_statistics.tsv：完整新检验，含不可评估记录。
tables/patient_program_effects.tsv、patient_site_roster.tsv：逐患者效应与配对细胞数。
tables/F3_threshold_sensitivity_summary.tsv、F3_patient_state_marker_expression.tsv：共检出敏感性和逐状态基因表达。
tables/background_QC.tsv、UMI_matching_ledger.tsv：最终修正的背景QC及匹配名单口径。
tables/reused_F3_full_transcriptome_statistics.tsv.gz及reused_*GSEA.tsv.gz：原全转录组结果与模型结果，保留原检验家族。
resources、code、objects、QA：完整集合、实施脚本、计数与背景缓存、独立验收材料。

图件01显示患者/谱系构成及冻结状态；02显示逐患者程序效应；03显示参照依赖；04显示质量背景；05显示分数计数阈值敏感性。复现顺序为analyze.py、summarize.py、sensitivity.py、report.py。sensitivity.py修正背景MT基因映射并保留共检出阈值敏感性，最终background_QC.tsv为交付口径。
'''
(O/'delivery/CAF_F3深入核查与F4表达分层报告_20261005.txt').write_text(text,encoding='utf-8-sig');print('REPORT_AND_FIGURES_COMPLETE',flush=True)
