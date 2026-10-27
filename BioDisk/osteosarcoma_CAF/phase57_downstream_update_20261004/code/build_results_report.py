from pathlib import Path
import json,hashlib,datetime
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase55_annotation_freeze_20261004';B=P/'phase56_patient_programs_20261004';O=P/'phase57_downstream_update_20261004'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
q=O/'independent_QA';assert (q/'independent_R_QA.exit_code').read_text().strip()=='0' and (q/'INDEPENDENT_R_QA_COMPLETE').is_file()
stats=read(B/'tables/all_patient_program_statistics.tsv.gz');source=read(A/'tables/source_candidate_distribution.tsv');dataset=source.drop_duplicates('source_id').set_index('source_id').dataset.to_dict() if 'source_id' in source else read(B/'tables/pseudobulk_sample_metadata.tsv').drop_duplicates('source_id').set_index('source_id').dataset.to_dict()
raw=stats[stats.scope.eq('agreed11')&stats.reference.eq('fixed_background')&stats.normalization.eq('raw')];corr=stats[stats.scope.eq('agreed11')&stats.reference.eq('fixed_background')&stats.normalization.eq('corrected')];reg=raw.merge(corr,on=['source_id','program'],suffixes=('_raw','_corrected'),validate='one_to_one');reg['raw_nominal_up']=reg.P_raw.lt(.05)&reg.mean_effect_raw.gt(0);reg['raw_nominal_down']=reg.P_raw.lt(.05)&reg.mean_effect_raw.lt(0);reg['nominal_same_direction_after_correction']=reg.P_raw.lt(.05)&reg.P_corrected.lt(.05)&(np.sign(reg.mean_effect_raw)==np.sign(reg.mean_effect_corrected));reg['direction_preserved_after_correction']=np.sign(reg.mean_effect_raw)==np.sign(reg.mean_effect_corrected);save(reg,'program_source_sensitivity_registry.tsv')
cs=reg.groupby('program',as_index=False).agg(sources_evaluable=('source_id','nunique'),raw_nominal_up_sources=('raw_nominal_up','sum'),raw_nominal_down_sources=('raw_nominal_down','sum'),nominal_same_direction_after_correction_sources=('nominal_same_direction_after_correction','sum'),direction_preserved_sources=('direction_preserved_after_correction','sum'));save(cs,'program_cross_source_descriptive_summary.tsv')
main=stats[stats.reference.eq('fixed_background')&stats.program.isin(['NABA_CORE_MATRISOME','NABA_COLLAGENS','F3_SELECTED_REMODELING','CUSTOM_mCAF_matrix','CUSTOM_iCAF_customer','CUSTOM_myCAF_customer','HALLMARK_TGF_BETA_SIGNALING'])].copy();save(main,'main_program_results_all_scopes.tsv')
gene_summary=read(O/'tables/gene_statistics_summary.tsv');gene_summary['dataset']=gene_summary.dataset.fillna(gene_summary.source_id.map(dataset));save(gene_summary,'gene_statistics_summary.tsv')
# Explicitly document HGNC-preferred mappings and preserved unresolved source symbols.
hg=read(P/'phase41_full_transcriptome_setup_20261001/resources/hgnc_complete_set.txt');approved=set(hg.loc[hg.status.eq('Approved'),'symbol']);mapping=[]
for f in (B/'resources').glob('D*_gene_feature_mapping.tsv'):
 m=read(f);mapping.append({'source_id':f.name.split('_')[0],'fixed_gene_units':len(m),'HGNC_approved_symbols':int(m.gene.isin(approved).sum()),'retained_other_source_symbols':int((~m.gene.isin(approved)).sum()),'policy':'Original source symbols retained when unresolved; no phenotype/P-based representative selection.'})
save(pd.DataFrame(mapping),'gene_mapping_summary.tsv')
# Export patient-level source plots, using separate sources rather than pooling their donors.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
F=O/'figures';F.mkdir(exist_ok=True);order=sorted(raw.source_id.unique());colors={'raw':'#276b8b','corrected':'#d28b37'}
fig,axes=plt.subplots(1,2,figsize=(13,5.5),sharey=True)
for ax,program,title in zip(axes,['NABA_CORE_MATRISOME','F3_SELECTED_REMODELING'],['Core matrisome','COL1A1 / LOX / MMP14 program']):
 z=stats[stats.scope.eq('agreed11')&stats.reference.eq('fixed_background')&stats.program.eq(program)]
 for i,src in enumerate(order):
  zz=z[z.source_id.eq(src)].set_index('normalization');ax.plot([zz.loc['raw','mean_effect'],zz.loc['corrected','mean_effect']],[i-.1,i+.1],color='#999999',lw=1)
  for norm,offset in [('raw',-.1),('corrected',.1)]:
   rr=zz.loc[norm];ax.scatter(rr.mean_effect,i+offset,s=45,color=colors[norm],marker='o' if rr.P<.05 else 'x',label=norm if i==0 else None)
 ax.axvline(0,color='#777777',lw=.8);ax.set_title(title);ax.set_xlabel('Mean patient log2 pseudobulk effect');ax.grid(axis='x',alpha=.15)
labels=[]
for src in order:
 n=int(raw[raw.source_id.eq(src)].N.iloc[0]);labels.append(f'{dataset[src]} (N={n})')
axes[0].set_yticks(range(len(order)),labels);axes[0].invert_yaxis();axes[1].legend(title='Count version',loc='lower right');fig.suptitle('Agreed 11-cluster candidates vs fixed background',y=.98);fig.text(.22,.025,'Filled circles: raw P < 0.05. Crosses: P >= 0.05. Complete P and FDR in statistical tables.',fontsize=9);fig.tight_layout(rect=[0,.06,1,.95]);fig.savefig(F/'patient_program_raw_corrected.png',dpi=180);fig.savefig(F/'patient_program_raw_corrected.pdf');plt.close(fig)
scopes=['agreed11','agreed11_no_flag','all15','all15_no_flag'];allsource=sorted(dataset)
x=stats[stats.reference.eq('fixed_background')&stats.normalization.eq('raw')&stats.program.eq('NABA_CORE_MATRISOME')].pivot(index='source_id',columns='scope',values='mean_effect').reindex(index=allsource,columns=scopes)
fig,ax=plt.subplots(figsize=(8,6));cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#dddddd');im=ax.imshow(x.to_numpy(),aspect='auto',cmap=cmap,vmin=0,vmax=1)
for i in range(len(x)):
 for j in range(len(scopes)):
  v=x.iloc[i,j];ax.text(j,i,'NA' if pd.isna(v) else f'{v:.3f}',ha='center',va='center',color='white' if np.isfinite(v) and v>.65 else 'black')
ax.set_yticks(range(len(x)),[dataset[s] for s in allsource]);ax.set_xticks(range(4),['Agreed 11','Agreed 11\nwithout CNV flag','All 15','All 15\nwithout CNV flag']);ax.set_title('Core matrisome: candidate-scope sensitivity');fig.colorbar(im,ax=ax,label='Mean patient log2 effect');fig.tight_layout();fig.savefig(F/'scope_sensitivity.png',dpi=180);fig.savefig(F/'scope_sensitivity.pdf');plt.close(fig)
shares=read(B/'tables/source_target_gene_RNA_shares.tsv');s=shares[shares.normalization.eq('raw')&shares.scope.isin(['agreed11','all15'])&shares.gene.isin(['COL1A1','LOX','MMP14'])].copy();s['column']=s.scope+'|'+s.gene;columns=[sc+'|'+g for sc in ['agreed11','all15'] for g in ['COL1A1','LOX','MMP14']];x=s.pivot(index='source_id',columns='column',values='median_RNA_share').reindex(index=allsource,columns=columns)
fig,ax=plt.subplots(figsize=(10,6));im=ax.imshow(x.to_numpy()*100,cmap='Blues',vmin=0,vmax=60,aspect='auto')
for i in range(len(x)):
 for j in range(len(columns)):
  v=x.iloc[i,j];ax.text(j,i,'NA' if pd.isna(v) else f'{v*100:.1f}%',ha='center',va='center',color='white' if v>.36 else 'black')
ax.set_yticks(range(len(x)),[dataset[s] for s in allsource]);ax.set_xticks(range(6),[z.replace('|','\n') for z in columns]);ax.set_title('Candidate contribution to total measured RNA\nMedian of source-scoped patient shares');fig.colorbar(im,ax=ax,label='RNA share (%)');fig.tight_layout();fig.savefig(F/'target_gene_RNA_contribution.png',dpi=180);fig.savefig(F/'target_gene_RNA_contribution.pdf');plt.close(fig)
# Report numerical results, including material negative/sensitivity findings.
def pp(v):return 'NA' if pd.isna(v) else f'{v:.6g}'
core=raw[raw.program.eq('NABA_CORE_MATRISOME')].sort_values('source_id');sh=shares[shares.normalization.eq('raw')&shares.scope.eq('all15')&shares.gene.isin(['COL1A1','LOX','MMP14'])]
lines=['CAF复核后分层注释、患者程序与下游更新报告','版本：phase55–57 / 2026-10-04','',
'一、分层注释已冻结，全部细胞保留，四种候选范围完整验收',
'已生成all_cell_annotation_v2.tsv.gz，覆盖9来源、107份样本、307簇、876,648细胞。保留原标签、盲审表达谱系、42簇深度复核结论、混合标记及CNV证据；未生成任何已确认CAF标签。',
'原11个成纤维样簇agreed11含45,721细胞；暂不纳入3个已有明确CNV风险标记簇的agreed11_no_flag含41,193；加入4个候选的all15含91,012；其no_flag版本含86,484。3个标记簇为D01_K017、D01_K036、D09_K008。no_flag仅是风险敏感性版本，不等于非恶性已获证明。患者特异CNV线索另列，不新设数值门槛。',
'13个盲审混合程序簇未进入上述成纤维样候选范围。全部细胞仍保留在注释主表；表达比较的固定背景为all15之外且不属于13个混合程序簇的细胞。四种范围使用同一背景，不将暂不纳入的风险簇移进对照。骨源样Osteo群作为独立命名的第二种同患者参照。',
'样本、患者×部位、患者和来源分布均已导出。不同来源的患者未合并为全局独立样本量。GSE299023没有这15个候选簇；SCPCP000023仅在扩展范围有候选；GSE162454的原11簇候选全部属于上述CNV标记簇，no_flag版本在该来源不可评估，不能解释为该来源生物学上不存在基质程序。','',
'二、核心基质程序的平均方向跨来源保留，校正后统计强度存在来源差异',
'84个既有基因程序在四种候选范围、两个参照和原始/DecontX两种表达版本中分析。144个预设情境中112个可评估，形成9,408条患者级程序结果。每名患者先在相同部位内汇总伪bulk，再对该患者的部位效应等权平均；统计单位是来源内患者。',
'以下为agreed11相对固定背景的NABA_CORE_MATRISOME原始表达主结果。效应是成员基因患者log2表达差的平均，不是细胞比例或单一基因倍数。']
for r in core.itertuples(index=False):
 level='FDR显著' if r.FDR<.05 else ('名义显著' if r.P<.05 else '未达名义显著')
 if r.N<10:level+='；小样本探索性'
 lines.append(f'{r.dataset}：N={r.N}，平均效应={r.mean_effect:.3f}，正向患者{r.positive_patients}/{r.N}，P={pp(r.P)}，FDR={pp(r.FDR)}；{level}。')
lines += ['七个有原11簇候选的来源平均方向均为正，各来源逐一留一患者后平均方向也保持为正。SCPCP000018与GSE152048的主结果在校正和CNV标记簇敏感性中继续保留。GSE318840为N=8的探索性结果，虽通过该家族FDR，仍不写成大样本确证。',
'SCPCP000017是重要敏感性结果：N=20，核心基质平均效应由原始0.285降至校正后0.199，P由0.0016升至0.0549，FDR由0.00217升至0.06232。方向保留，但校正后不满足P<0.05，不能列入校正后名义显著信号。相对Osteo参照时，该程序原始/校正效应0.433/0.380，P=0.00015/0.0009，表明参照选择也影响证据强度。',
'GSE250015的N=2和GSE270231的N=5给出方向及效应；全正向时二侧精确置换最小P也分别只能达到0.5和0.0625。其未显著不用于抹去观察方向，也不据此扩充有效样本量。','',
'三、扩展候选改变效应幅度，新增来源并未自动产生显著基质信号',
'SCPCP000018扩展至all15后核心基质效应0.833→0.809，N=10、方向10/10，P=0.001953，FDR=0.007813；主体方向保留。',
'SCPCP000017扩展后核心基质效应0.285→0.336，正向患者17/20→19/20，原始P=0.00005、FDR约0.000168；同一来源的三基因程序效应却从2.555降至1.137。这是候选组成改变的敏感性，不是新增独立重复。',
'SCPCP000023扩展候选可比较N=25，但核心基质效应仅0.073，17/25正向，P=0.3061、FDR=0.37264；三基因程序效应0.327，P=0.1645。当前数据不支持该来源这两个预设程序达到名义显著。低计数及CNV覆盖不足的边界继续保留。',
'GSE152048暂不纳入D09_K008后核心基质效应0.597→0.559，N=11，仍11/11正向，P=0.0009766、FDR约0.003729。其他来源没有对应的3个CNV标记簇时，no_flag与包含版本相同；不能把这些重复结果当作独立验证。','',
'四、三基因重塑程序和TGFβ程序均具有来源差异，不能代表统一CAF机制',
'COL1A1/LOX/MMP14程序在SCPCP000018、SCPCP000017和GSE152048的原始与校正后均保持正向且P<0.05，并通过相应FDR；GSE162454为N=6的名义显著探索性信号，FDR约0.0597。',
'GSE318840三基因平均效应-0.087，P=0.7656；GSE270231为-0.108，P=0.875。校正后方向也未转为普遍正向。不能把全体成纤维样候选统一称为这一重塑程序的高表达群。',
'TGFβ程序在SCPCP000017的agreed11主比较中为负向（效应-0.581，P=0.00005、FDR约0.0001），其他来源方向不同。完整84程序的正、负及不显著结果均保留。',
'在原始与校正后均保留正向P<0.05的来源数中，既有iCAF_customer面板、mCAF_matrix面板和EMT Hallmark各覆盖5个来源。这些是预定义表达程序的重叠富集，不据此重新命名稳定iCAF/myCAF亚型或证明发生EMT；完整逐来源P、FDR和效应见program_source_sensitivity_registry.tsv。','',
'五、新候选范围的差异与通路已更新，旧F3/F4结果只保留原操作标签含义',
'已对四种候选范围相对固定背景、原始表达开展来源内全共同测量基因/feature的患者检验，共28个可评估情境；另8个预设情境无候选或无匹配患者，状态表完整保留。基因映射在查看P值前确定，优先HGNC对应feature，无法解析的源符号保留并单列覆盖。',
'原11簇主比较的完整基因结果概览如下，重复候选范围不作为新增发现数累加。']
for r in gene_summary[gene_summary.scope.eq('agreed11')].sort_values('source_id').itertuples(index=False):
 if r.N==0:lines.append(f'{r.dataset}：无匹配候选患者，未作推断检验。')
 else:lines.append(f'{r.dataset}：N={r.N}，有效检验{int(r.valid_tests)}，原始P<0.05且平均效应正向{int(r.nominal_up)}、负向{int(r.nominal_down)}，FDR<0.05共{int(r.FDR_genes)}。')
lines += ['所有原始P<0.05基因均进入候选列表，不要求FDR显著。按患者平均效应方向分别进行Hallmark、NABA Matrisome和Reactome过度代表分析；三种库、方向及来源/候选范围分别形成FDR家族。背景为预先固定映射后有效检验的基因单位，不在代表feature不可评估时改选一个更显著feature。全结果包含无名义基因时的不可评估通路行。',
'患者级程序检验与通路ORA分开。ORA是基因集合富集模型统计，不是第二批独立患者的生物效应验证。旧F3/F4的差异、TF、LR、临床和外部映射结果保留为历史操作标签结果，不直接转移到新候选集合。','',
'六、空间数值保留，来源解释更新为基质/RNA关联，不能归为已确认CAF功能',
'固定空间数据、预设基因集和已验收模型均未变化，因此phase46/50的现有相关数值保留，记录原文件SHA256。三基因程序与CD8 RNA的调整后关联N=19、平均效应-0.02970、P=0.003723、FDR=0.009308；邻域CD8为-0.03289、P=0.0005646、FDR=0.001694。',
'核心基质程序的CD8调整后关联却为-0.007635，P=0.3711、FDR=0.3976；邻域关联P=0.5949。当前数据不支持把三基因的负相关推广为整个核心基质程序的普遍CD8抑制。',
'新候选集合的单细胞RNA贡献进一步限制谱系归因。以下为all15在各来源患者等权后的中位RNA份额，分母是该患者/部位全部输入细胞中对应基因的测得RNA：']
for src in sorted(sh.source_id.unique()):
 g=sh[sh.source_id.eq(src)].set_index('gene');lines.append(f'{dataset[src]}：COL1A1 {g.loc["COL1A1","median_RNA_share"]:.1%}，LOX {g.loc["LOX","median_RNA_share"]:.1%}，MMP14 {g.loc["MMP14","median_RNA_share"]:.1%}。')
lines += ['例如SCPCP000017加入扩展候选后，三基因中位RNA份额从10.7%/3.8%/5.7%升至40.3%/30.8%/30.0%，说明来源解释对候选范围敏感；其余RNA仍来自候选之外。零份额代表该来源在当前候选定义下没有相应候选贡献，不代表该基因未表达。','',
'七、Methods',
'数据：phase53原始匿名H5AD及对应键、phase54复核附加表、冻结盲审及揭盲对照；全部107份样本和876,648细胞进入版本化主注释。复用已成功完成且细胞/基因顺序核对一致的DecontX计数。原始及校正文库逐样本核对，四个范围的总数与注释验收一致。',
'基因空间：每来源使用该来源所有样本共同测量的feature；不存在全来源交集限制。每符号优先HGNC批准Ensembl对应feature，其后按feature ID确定；保留未解析源符号并记录。完整文库分母包含该细胞全部测得feature，即使某feature不在共同分析集合中。',
'程序：冻结的50 Hallmark、22 NABA、4个既有功能补充程序及8个客户面板，共84。每个来源的成员覆盖与缺失名单完整记录，无新增覆盖比例门槛。',
'伪bulk：同来源同患者同部位内，分别合并目标及固定参照的计数与全文库，计算log2(CP10k+0.1)差；患者多部位对效应等权平均。程序分数为该来源实际测量成员的患者基因效应均值。组间不共享目标细胞，四个候选范围的背景保持固定。',
'检验：二侧患者平均效应符号置换；N<=15枚举全部2^N，N>15用固定种子20261004的19,999次随机置换并加入观测项，P=(极端次数+1)/20,000。不足2名患者、全零效应或不可测结果记NA。原始P<0.05为名义筛选硬标准；全部P、FDR、患者数、方向和效应保留。',
'多重检验：程序在来源×候选范围×参照×原始/校正版本内的84个预定义程序家族进行BH；基因每个对应原始表达对比为独立家族；ORA按来源×范围×库×方向分别BH。不可评估项保留NA，校正时占预定义家族位置。未将不同类型分析混入同一校正池。',
'敏感性：原始/校正版本、四种候选范围、固定背景/Osteo两种参照、留一患者；候选范围的效应变化另在双方共同患者上比较，避免把患者组成变化当作单纯细胞范围效应。',
'验收：独立R实现复算全部9,408程序P值及504个预定/可复现抽查基因P值，最大差0；全部308个程序、基因、ORA家族共1,137,474行FDR独立核对，最大浮点差2.22e-16。包括已知精确P值、N=1、全零效应、RNA贡献范围与注释细胞集合检查。','',
'八、Discussion与使用边界',
'本轮建立了可追溯候选注释及患者级表达程序结果。表达定义的候选群与同一数据中的基因/程序统计存在选择依赖，因此P与FDR用于候选发现和范围内比较，不是独立CAF身份验证。N<10按探索性呈现；来源间患者对应关系未完全解决，未合并N或宣称多次独立外部验证。',
'原始计数、核/细胞技术、部位和细胞捕获组成影响RNA贡献及相对表达。DecontX是无空液滴背景的模型敏感性分析；未纳入3个风险簇也不建立其余群的非恶性身份。单细胞RNA份额不能直接等同于空间spot的谱系构成。',
'固定基因集的空间相关支持的是局部RNA关联，不能写成CAF造成免疫排斥、真实保护/抑制或临床可用性。程序资源名称中的iCAF/myCAF、EMT、其他癌种及转移用语是原基因集名称，不据此确立对应亚型、转分化、来源癌种或转移功能。',
'本次未重新计算整批CNV，也未以新集合自动替换历史临床、TF、通讯、外部映射或分选结论；这些模块在结果索引中明确标为历史操作标签结果。当前新注释、患者程序、差异及通路结果以本次版本为准。','',
'九、交付与当前使用版本',
'phase55/tables/all_cell_annotation_v2.tsv.gz：全细胞复核后注释；cluster_annotation_v2.tsv：307簇并行标签；sample/patient_site/patient/source_candidate_distribution：分布。',
'phase56/tables/all_patient_program_statistics.tsv.gz：全部程序效应/P/FDR；patient_program_effects：逐患者值；scope_sensitivity_common_patients、correction_sensitivity：范围和校正敏感性；patient/source_target_gene_RNA_shares：来源贡献。',
'phase57/tables/genes与pathways：新候选差异及通路完整结果；main_program_results_all_scopes和program_source_sensitivity_registry：主结果及敏感性索引；DOWNSTREAM_RESULT_INDEX：当前与历史模块的使用说明。',
'三个阶段的冻结资源、代码、退出码、独立R验收和图件均纳入交付。旧正式报告与统计包按原版本保留，本报告为复核后新的分析补充，不作已确认CAF或临床验证结论。']
report=O/'delivery/CAF复核后分层注释与患者分析报告_20261004.txt';report.write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
(O/'REPORT_READY').write_text('SUCCESS\n');print(json.dumps({'report':str(report),'figures':3,'program_registry':len(reg),'report_ready':True},ensure_ascii=True))