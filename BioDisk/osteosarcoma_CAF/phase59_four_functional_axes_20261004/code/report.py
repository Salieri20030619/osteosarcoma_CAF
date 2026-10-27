from pathlib import Path
import json,pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase59_four_functional_axes_20261004')
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
s=read(O/'tables/axis_statistics.tsv');ob=read(O/'tables/frozen_cells_with_identity.tsv.gz');sets=json.loads((O/'resources/GENE_SETS.json').read_text());ds='SCPCP000017';four=['inflammation','contraction','matrix','osteogenic']
def tab(d,cols):return d[cols].to_csv(sep='\t',index=False,float_format='%.6g',na_rep='NA')
def pick(context='main',mode='corrected',programs=None,contrasts=None):
 d=s[(s.dataset==ds)&(s.context==context)&(s.matrix==mode)]
 if programs is not None:d=d[d.program.isin(programs)]
 if contrasts is not None:d=d[d.contrast.isin(contrasts)]
 return d
cols=['contrast','program','N','mean_delta','positive_patients','P','FDR']
# Additional figure directly exposes dependence on program definition.
fig,aa=plt.subplots(1,2,figsize=(11,5.5),layout='constrained')
for ax,cl,progs in zip(aa,['F2','F4'],[['contraction','contractile_markers','matrix','matrix_markers'],['osteogenic','osteogenic_effectors','Liu_C2_OMD_OGN','mineralization']]):
 for j,(ctx,color,shift) in enumerate([('main','#276F9B',-.12),('Fibro_only','#C65156',.12)]):
  d=pick(ctx,programs=progs,contrasts=[cl+'_vs_rest']).set_index('program').loc[progs]
  ax.scatter(d.mean_delta,np.arange(len(progs))+shift,c=color,label=ctx.replace('_',' '),s=45)
  for k,r in enumerate(d.itertuples()):
   marker='**' if r.FDR<.05 else '*' if r.P<.05 else '';ax.text(r.mean_delta+.07,k+shift,f'N={r.N}{marker}',fontsize=9,va='center')
 ax.axvline(0,color='#999999',lw=.8);ax.set_yticks(range(len(progs)),[p.replace('_',' ') for p in progs]);ax.invert_yaxis();ax.set_xlabel('Mean patient effect, corrected counts');ax.set_title(cl+' full sets and reference panels',weight='bold');ax.margins(x=.35);ax.legend(frameon=False,loc='best')
fig.supxlabel('Marker panels are prespecified secondary evidence. * P < 0.05; ** FDR < 0.05 within their separate hypothesis families.',fontsize=9);fig.savefig(O/'figures/07_program_definition_comparison.png',dpi=180);plt.close(fig)
v=read(O/'tables/patient_four_axis_overlap.tsv');co=[]
for (source,mode,contrast),g in v.groupby(['dataset','count_matrix','contrast']):
 ev=g[g.all_four_evaluable];co.append({'dataset':source,'matrix':mode,'contrast':contrast,'N_complete':len(ev),'contraction_and_matrix_positive':int((ev.contraction.gt(0)&ev['matrix'].gt(0)).sum()),'matrix_and_osteogenic_positive':int((ev['matrix'].gt(0)&ev.osteogenic.gt(0)).sum()),'two_or_more_positive_axes':int(ev.n_positive_axes.ge(2).sum())})
co=pd.DataFrame(co);co.to_csv(O/'tables/patient_axis_coactivity_summary.tsv',sep='\t',index=False)
# Technical detectability description; patients, not cells, remain inferential units.
cd=read(O/'tables/cell_axis_codetection.tsv.gz');cd['coverage_eligible']=cd.assayed_panel_genes.ge(np.ceil(.8*cd.definition_genes));cv=cd[cd.matrix.eq('corrected')].pivot(index='cell_id',columns='axis',values='coverage_eligible');cp=cd[cd.matrix.eq('corrected')].pivot(index='cell_id',columns='axis',values='detected_panel_genes').rename(columns={'matrix':'matrix_detected'});cp=cp.join(ob.set_index('cell_id')[['dataset','class_id','donor_record','lineage_stratum']]);cp['contractile_ECM_evaluable']=cv.contraction&cv['matrix'];cp['ECM_osteogenic_evaluable']=cv['matrix']&cv.osteogenic;cp['contractile_ECM_codetected']=cp.contraction.ge(2)&cp.matrix_detected.ge(2)&cp.contractile_ECM_evaluable;cp['ECM_osteogenic_codetected']=cp.matrix_detected.ge(2)&cp.osteogenic.ge(2)&cp.ECM_osteogenic_evaluable;cp.groupby(['dataset','class_id','donor_record','lineage_stratum']).agg(cells=('class_id','size'),contractile_ECM_evaluable=('contractile_ECM_evaluable','sum'),ECM_osteogenic_evaluable=('ECM_osteogenic_evaluable','sum'),contractile_ECM_codetected=('contractile_ECM_codetected','sum'),ECM_osteogenic_codetected=('ECM_osteogenic_codetected','sum')).reset_index().to_csv(O/'tables/patient_lineage_codetection_summary.tsv',sep='\t',index=False)
name='CAF四功能轴与F2F4信号来源分析报告_20261004.txt'
text='''骨肉瘤 CAF 四功能轴与 F2 F4 信号来源分析报告
2026-10-04  phase59

本轮保留原始F1–F5，完成炎症、收缩、基质、成骨相关四个可重叠功能轴的患者级检验，并对F2/F4进行表达谱系分解。结果显示：F2经典收缩标志升高，但完整基质程序未整体升高；F4具有成骨/矿化相关转录信号，主要成骨转录来自Osteo表达层，Fibro层仅保留部分OMD/OGN线索。四轴结果不支持按客户四分类直接重命名F1–F5。

方法

研究对象沿用phase58冻结的8174个细胞：F1 695、F2 225、F3 251、F4 5304、F5 1699，来自9个来源。复用原始全基因矩阵、已有DecontX校正矩阵、phase55逐细胞复核注释以及107样本体系中的既有RNA-CNV证据。没有新聚类、重跑CNV或改变原细胞标签。本轮同患者谱系背景对照实际纳入原始8174以外的204985个细胞。

主轴定义来自phase58冻结的完整canonical基因集，未经按本轮结果挑选：炎症为HALLMARK_INFLAMMATORY_RESPONSE（200基因）；收缩为REACTOME_SMOOTH_MUSCLE_CONTRACTION（43）；基质为NABA_CORE_MATRISOME（275）；成骨相关为GOBP_OSTEOBLAST_DIFFERENTIATION（267）。八项次级完整通路、十一项文献/谱系参考面板，以及四个去除主轴间共享基因后的敏感性集合，共27个程序、1300个联合基因。成骨轴衡量相关基因集的相对表达，不把集合中每个调节因子的升高都解释为成骨功能增强。

按来源、独立donor_record和部位聚合，使用相应矩阵的完整library分母计算CP10k。效应为每个基因log2(CP10k+0.1)的目标组减参照组差值，再在预定集合内平均；同患者多部位等权平均后，以患者等权计算效应与检验。主比较为每个F群对其余四群，身份限制比较在两侧同时限于Fibro。外部背景为原始8174之外同患者同部位的Fibro/Mural/Osteo表达层，仅使用原始计数。

基因需在组成对比的全部样本中被测量，缺测不补零；程序至少覆盖定义的80%，且至少min(3,集合大小)个基因才评估。正常测量但未表达保留为0。主分析要求对比两侧各至少1细胞，另登记每侧至少10细胞的敏感性分析。

P为患者效应双侧符号翻转检验：N≤18枚举全部符号，N>18使用99999次随机符号并加入观测项；N<2或全零效应不检验。每来源、计数版本和比较情境内，四主轴、次级完整通路、次级面板、非重叠敏感性分别定义FDR家族；五个F群主比较共同进入相应家族。基因检验另设家族，所有预设不可评估行保留，BH分母按预设行数处理。以目标效应原始P<0.05筛选名义候选，并同时报告FDR。DecontX校正与FDR校正分别记录。

RNA来源以同患者同部位F群内每个基因的转录本贡献计算，先按部位平均、再按患者等权平均；不存在的谱系贡献为0，整群该基因全零时记不可评估。该分母与F群对其余四群的配对分母不同。图中N随基因可检测表达变化。另保留汇总RNA总量比例。细胞共检出仅作描述：每个参考面板至少两个基因、各计数≥1，不做细胞级推断。

结果

F2 保留经典收缩表型，完整基质程序未整体升高

在主要可配对来源SCPCP000017中，F2与其余四群比较N=12。DecontX后，经典收缩面板效应+2.330，11/12患者同向，P=0.002441、FDR=0.014920；文献基质面板效应+0.766，10/12同向，P=0.005859、FDR=0.021484。这确认特定收缩装置及基质标志的升高。

完整集合给出更窄的解释边界：43基因收缩通路效应+0.270，8/12同向，P=0.275879、FDR=0.459798；275基因核心基质集效应−0.048，4/12同向，P=0.719727、FDR=0.833333。当前数据不支持这两个完整集合在F2整体上升。此结果与phase58少量客户标志构成的myCAF/基质面板不同，必须同时保留，不能用短面板阳性代替完整程序结论。

限制到Fibro表达层后N=10，收缩标志仍呈+1.902，9/10同向，P=0.009766、FDR=0.107422，属于名义显著；基质标志+0.649、P=0.138672、FDR=0.346680。说明F2的收缩标志不能全部归因于Mural混入，但Fibro层的基质差异未获名义显著支持。更严格的Fibro筛选在该主来源未进一步改变此处结果，不能将二者算为两次独立复制。

来源分解也支持Fibro参与：DecontX后，F2内ACTA2和MYL9的患者等权RNA贡献分别有67.1%和64.9%来自Fibro；COL1A1和POSTN分别为82.3%和94.7%。对应可评估N分别10、10、13、9。Mural亦有收缩转录贡献，不能由这些比例确认CAF身份。

同患者外部背景对照中，F2-Fibro完整收缩轴低于BG-Mural（N=12，效应−0.858，P=0.000488、FDR=0.004688）；完整核心基质集也没有显示超越广义Fibro背景的优势。TGFβ Hallmark在F2为下降（效应−1.002，P=0.000977、FDR=0.003005）；LOX家族与LRRC15也没有获得名义上升支持。F2适合保留为“收缩标志较高的候选表达群”，现阶段不定义为全局最强基质生产型或LRRC15/LOX依赖亚型。

F4 的主要成骨转录贡献来自 Osteo 表达层

全体F4在主来源的完整成骨轴呈名义上升：N=15，效应+0.322，10/15同向，P=0.038025、FDR=0.095062。成骨效应标志面板+1.302，11/15同向，P=0.004089、FDR=0.017301；OMD/OGN面板+0.949，12/15同向，P=0.003174、FDR=0.015869。

次级完整GOBP_BONE_MINERALIZATION集合（127基因）在F4同样升高：效应+0.436，12/15同向，P=0.003967、FDR=0.009918。这一转录信号通过其预设FDR家族，必须保留；它与phase58较短矿化面板的结果口径不同，不能相互替代，也不能直接称为矿化沉积增加。限制于Fibro后，该完整矿化集合为+0.266（N=11，P=0.248047、FDR=0.522204），没有名义显著支持。

按既有表达复核，F4全来源共有Osteo层3300细胞（含3个Chondro）、Fibro层909、Mural层47，其余1048。只按原始Osteo标签计为3297。DecontX后主来源F4内RUNX2、ALPL、SP7、IBSP的患者等权RNA贡献分别有71.7%、77.9%、76.0%、75.4%来自Osteo层；对应可评估N为15、15、13、14。原始与校正后图形方向相近。

在F4内部同患者配对的Fibro对Osteo比较中，完整成骨轴在Fibro更低（N=12，效应−0.310，P=0.000488、FDR=0.001465）；Fibro的文献基质面板更高（+0.860，P=0.000488、FDR=0.006714）。这支持F4内部的基质样与成骨样表达分量可以区分。

将F4对其余四群的比较双方均限制于Fibro后，完整成骨轴为+0.178（N=11，7/11同向，P=0.472656、FDR=0.833333），当前数据不支持该层具有完整成骨程序优势；OMD/OGN仍为+1.045（7/11同向，P=0.023438、FDR=0.117188），保留为名义显著候选线索。OMD的RNA贡献也以Osteo为主，OGN则在Fibro/Osteo间更分散，所以这两个基因不能用作排他的CAF身份标签。

F4-Fibro的完整成骨轴也低于同患者BG-Osteo（N=12，效应−0.547，P=0.000488、FDR=0.001172）。本轮确定的是表达来源构成，不把Osteo层直接等同于病理确诊恶性细胞，也不将F4整群命名为osteogenic CAF。

F1 F3 F5 的四轴结果保持各自定位

F1的炎症轴未获得名义上升支持；收缩与成骨完整集较低。F3主来源只有3位配对患者，四轴均可报告效应方向，但本检验最小双侧P=0.25，不能建立稳定分型。F5四主轴均较低，当前数据不支持活跃基质生产或成骨型解释。

以下为主来源DecontX后完整四轴结果。mean_delta为患者平均效应；positive_patients为效应大于0的人数。
'''
text+=tab(pick(programs=four),cols)
text+='''
来源差异与敏感性分析改变证据强度

九个来源单独保留。F2收缩完整集在部分来源正向、部分来源负向，多数外部来源有效N为1–4；F4完整成骨集也出现跨来源方向差异。因此本轮没有把主来源结论写成全部队列一致复现。

每侧至少10细胞的敏感性中，F2由12位缩至3位、F4由15位缩至7位。F2收缩/基质正向但无法取得P<0.05；F4成骨方向保留。F4矿化集在该子集出现名义上升（N=7，效应+0.385，7/7同向，P=0.015625、FDR=0.208333），作为探索性信号完整保留，不能升级为物理矿化证据。

四个去重叠基因集、逐患者效应及留一患者后的效应范围均已输出。共同基因被剔除后的结果用于判断结论对集合重叠的依赖，不作为重新筛选较好P值的机会。

完整基因集的跨来源结果（DecontX后）
'''
text+=tab(s[(s.matrix=='corrected')&(s.context=='main')&(((s.contrast=='F2_vs_rest')&s.program.isin(['contraction','matrix']))|((s.contrast=='F4_vs_rest')&(s.program=='osteogenic')))],['dataset']+cols)
text+='\n主来源非重叠集合敏感性\n'+tab(pick(programs=[k+'_nonoverlap' for k in four],contrasts=['F2_vs_rest','F4_vs_rest']),cols)
text+='\n主来源其他完整程序\n'+tab(pick(programs=['collagen_formation','ECM_degradation','ECM_regulators','mineralization','ossification','TGF_beta'],contrasts=['F2_vs_rest','F4_vs_rest']),cols)
text+='\n主来源逐患者四轴共现描述\n'+tab(co[(co.dataset==ds)&(co.matrix=='corrected')],list(co.columns))
text+='''
讨论

本轮结果将“某些标志升高”与“完整功能程序一致增强”区分开。F2的收缩装置表型得到支持，ECM特征依赖基因集与身份范围；F4的成骨方向主要由Osteo表达分量贡献，Fibro层OMD/OGN线索值得保留，但不能据此认定成骨CAF。

原始8174对象已经过历史CAF/CNV规则筛选，其中的RNA-CNV“未检出异常”分布具有选择条件，不能再次充当独立非恶性证明。Fibro/Osteo/Mural是先前RNA复核标签，与当前表达不是独立实验验证；低RNA细胞、双细胞未全面评估及去污染模型均影响解释。保守Fibro筛选保留未评估双细胞，未将其当作阴性。外部背景使用原始计数，不能与候选细胞DecontX结果混为同一矩阵的直接对比。

同一患者多个样本或部位不增加N；来源之间的全局患者去重关系未确立，因此不合并为更大独立队列。少于10位患者的结果按探索定位。有限细胞数造成的稀疏效应与阈值敏感性已完整保留，未以细胞量扩充统计样本数。

本轮不涉及真实收缩力、ECM蛋白沉积、矿化测量、空间壳层几何或化疗反应，转录轴不能替代这些终点。所列背景细胞并非正常组织金标准。后续若推进细胞来源/蛋壳问题，应使用本轮分层索引保留F2-Fibro与Mural的对照，以及F4-Fibro、F4-Osteo和未定层的对照。

验收与文件索引

独立R从逐患者效应重新计算全部有效P：5837项程序检验、11350项标志基因检验，最大P误差0；504个程序FDR家族和126个基因FDR家族最大FDR误差1.11e-16。主来源四轴患者效应另有独立代数重建，RNA贡献与细胞比例均通过守恒检查。

tables/axis_statistics.tsv：全部18468个预设程序结果，含不可评估记录与原始P/FDR。
tables/marker_statistics.tsv.gz：全部38304个预设标志基因结果。
tables/patient_axis_effects.tsv.gz：逐患者效应；paired_patient_site_roster.tsv.gz：实际配对名单。
tables/frozen_cells_with_identity.tsv.gz：8174细胞原标签及复核层。
tables/F2F4_lineage_RNA_contribution.tsv.gz：逐患者部位的RNA归属；lineage_RNA_patient_equal_summary.tsv和lineage_RNA_pooled_summary.tsv：两种汇总口径。
tables/patient_four_axis_vectors.tsv、patient_axis_coactivity_summary.tsv：连续四轴及共现描述。
tables/patient_lineage_codetection_summary.tsv：细胞面板共检出描述，无细胞级P。
resources/GENE_SETS.json、gene_set_membership.tsv、primary_axis_overlap.tsv：完整基因集、来源及重叠；PROTOCOL.json：预设口径；INPUT_MANIFEST.json：输入指纹。
objects/selected_gene_counts.npz：本轮1300基因计数、完整分母和测量掩码；background_pseudobulks.npz：外部背景聚合对象（内部可信numpy对象格式）。
QA/：独立R、患者效应重建及验证日志。

图件
01_four_axes.png：完整四轴原始/去污染对照。
02_F2_F4_sensitivity.png：细胞数与身份分层敏感性。
03_F2_RNA_sources.png：F2逐基因RNA来源。
04_F4_RNA_sources.png：F4逐基因RNA来源。
05_source_consistency.png：完整集的来源差异与N。
06_external_competitors.png：同患者外部谱系竞争。
07_program_definition_comparison.png：完整集与文献面板的区别。

研究决策

保留F1–F5。F2按收缩标志较高的候选表达群推进，并保留ECM广谱结论不成立及Mural竞争；F4分开呈现Osteo主导信号与Fibro层OMD/OGN候选，不使用整体成骨CAF命名。本轮授权的四轴检验与F2/F4来源解析完成。
'''
# The narrative constants are checked against the saved statistical table below.
check=pick(programs=['osteogenic_effectors'],contrasts=['F4_vs_rest']).iloc[0]
assert abs(check.P-0.00408935546875)<1e-12
(O/'delivery'/name).write_text(text,encoding='utf-8-sig');print('REPORT_COMPLETE',name,flush=True)
