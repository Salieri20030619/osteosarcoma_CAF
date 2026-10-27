from pathlib import Path
import json,datetime,shutil
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase8_original_plan_20260919';D=P/'phase7_review_states_20260917'
assert (O/'BULK_COMPLETE').exists() and (O/'VALIDATION_COMPLETE').exists()
inv=pd.read_csv(O/'tables/analysis_sample_inventory.tsv',sep='\t')
program=pd.read_csv(O/'tables/strict_programs_raw_corrected.tsv',sep='\t')
r=program[(program.group=='all_strict')&(program.version=='raw')].copy()
comparisons=[]
for (dataset,name),g in r.groupby(['dataset','program']):
 lung=g[g.site.str.lower()=='lung']
 bone=g[g.site.str.lower().isin(['bone','femur','tibia','ulna','humerus'])]
 if len(lung) and len(bone):
  comparisons.append({'dataset':dataset,'program':name,'lung_samples':len(lung),'bone_samples':len(bone),'lung_mean':lung.mean_log1p_CP10K.mean(),'bone_donor_median':bone.mean_log1p_CP10K.median(),'lung_minus_bone_median':lung.mean_log1p_CP10K.mean()-bone.mean_log1p_CP10K.median(),'inference':'descriptive_only; one lung CAF donor per dataset; treatment not separately identified'})
pd.DataFrame(comparisons).to_csv(O/'tables/bone_lung_descriptive_programs.tsv',sep='\t',index=False)
# Exact genotype mapping and source labels remain separate by cohort.
old=pd.read_csv(D/'review/active_patient_tissue_sample_library_master.tsv',sep='\t')
new=pd.read_csv(P/'geo_preprocess_20260918/manifest.tsv',sep='\t')
pd.concat([old,new],ignore_index=True,sort=False).to_csv(O/'tables/source_patient_tissue_library_master.tsv',sep='\t',index=False)
my=pd.read_csv(D/'analysis/mycaf_subcluster_patient_counts.tsv',sep='\t').set_index('patient_id')
myreview=pd.DataFrame({'subcluster':my.columns,'cells':my.sum(0).to_numpy(),'patients_with_ge20_cells':(my>=20).sum(0).to_numpy(),'patients_with_ge10_cells':(my>=10).sum(0).to_numpy()})
myreview.to_csv(O/'tables/myCAF_subcluster_acceptance.tsv',sep='\t',index=False)
ab=[
 {'target':'FAP','clone':'427819','evidence':'manufacturer validates human flow cytometry; assess fresh OS tissue separately','url':'https://www.rndsystems.com/products/human-fibroblast-activation-protein-alpha-fap-apc-conjugated-antibody-427819_fab3715a'},
 {'target':'THY1/CD90','clone':'5E10','evidence':'human cell-surface flow quality tested; not fibroblast-specific','url':'https://www.biolegend.com/en-gb/products/brilliant-violet-711-anti-human-cd90-thy1-antibody-15886'},
 {'target':'PDGFRA/CD140a','clone':'alphaR1','evidence':'used in primary study of human bone marrow stromal cell identification; tissue-specific titration required','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8252541/'},
 {'target':'CD45/CD31/EpCAM','clone':'HI30 / WM59 / EBA-1','evidence':'official human antibody panel lists clones; no assumption that AbSeq conjugates directly substitute a flow reagent','url':'https://www.bdbiosciences.com/en-de/products/reagents/single-cell-multiomics/bd-omics-one-protein-panels/bd-omics-one-tumor-protein-panel'}]
pd.DataFrame(ab).to_csv(O/'tables/antibody_evidence.tsv',sep='\t',index=False)
bulk=pd.read_csv(O/'tables/GSE21257_primary_associations.tsv',sep='\t')
cv=pd.read_csv(O/'tables/transfer_donor_CV.tsv',sep='\t')
effects=pd.read_csv(O/'tables/full_customer_gene_effects.tsv',sep='\t')
rank=pd.read_csv(O/'tables/surface_combination_ranking.tsv',sep='\t')
lung=pd.read_csv(O/'tables/lung_reference_phenotype.tsv',sep='\t')
fig,ax=plt.subplots(figsize=(10,5))
for i,row in bulk.iterrows():
 ax.errorbar(row.OR_per_SD,i,xerr=[[row.OR_per_SD-row.CI_low],[row.CI_high-row.OR_per_SD]],fmt='o',color='#15786b' if row.FDR_BH<.05 else '#566579',capsize=3)
ax.set_yticks(range(len(bulk)),bulk.program);ax.invert_yaxis();ax.axvline(1,color='gray',ls='--');ax.set_xscale('log');ax.set_xlabel('Odds ratio per 1 SD (95% CI)');ax.set_title('GSE21257 | 53 patients | original metastasis endpoint\nGreen: BH FDR < 0.05 across 8 frozen programs')
fig.tight_layout();fig.savefig(O/'figures/bulk_association.png',dpi=160,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4))
ax.barh(cv.heldout_donor.str.split('__').str[-1],cv.balanced_accuracy,color='#566579')
ax.axvline(.6,color='#a5483f',ls='--',label='Predefined acceptance 0.60');ax.set_xlim(0,1);ax.set_xlabel('Balanced accuracy');ax.set_title('Discovery donor-held-out transfer check | mean 0.561');ax.legend(loc='lower right');fig.tight_layout();fig.savefig(O/'figures/transfer_validation.png',dpi=160);plt.close(fig)
mat=effects.pivot(index='state',columns='gene',values='median_log2CPM_difference')
fig,ax=plt.subplots(figsize=(15,3.2));v=max(abs(mat.to_numpy()).max(),1)
im=ax.imshow(mat,cmap='RdBu_r',vmin=-v,vmax=v,aspect='auto');ax.set_xticks(range(len(mat.columns)),mat.columns,rotation=90);ax.set_yticks([0,1],['H0 (4 paired donors)','H2 (3 paired donors)']);ax.set_title('Full customer gene panels | within-donor state minus other CAF');fig.colorbar(im,ax=ax,label='Median log2(CPM+1) difference');fig.tight_layout();fig.savefig(O/'figures/customer_gene_effects.png',dpi=160);plt.close(fig)
top=rank[(rank.target=='CAF_H0')&rank.positive_genes.str.contains('+',regex=False)].head(6).copy()
fig,ax=plt.subplots(figsize=(9,4));x=np.arange(len(top))
ax.bar(x-.18,top.median_RNA_target_coverage,width=.36,label='RNA H0 coverage',color='#258573')
ax.bar(x+.18,top.median_RNA_known_CAF_fraction_gate,width=.36,label='Known CAF fraction among RNA gate',color='#536f92')
ax.set_xticks(x,top.positive_genes,rotation=20,ha='right');ax.set_ylim(0,1);ax.set_ylabel('Median donor fraction');ax.set_title('Candidate surface pairs | RNA screen, not measured FACS purity');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(O/'figures/marker_pairs.png',dpi=160);plt.close(fig)
def tab(df):
 z=df.copy()
 for c in z:
  if pd.api.types.is_float_dtype(z[c]):z[c]=z[c].map(lambda v:'' if pd.isna(v) else f'{v:.4g}')
 return '| '+' | '.join(map(str,z.columns))+' |\n| '+' | '.join(['---']*len(z.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in z.itertuples(index=False,name=None))
sampleview=inv[inv.strict_CAF>0][['sample','dataset','site','strict_CAF']].copy();sampleview['sample']=sampleview['sample'].str.split('__').str[-1]
key=effects[effects.gene.isin(['RUNX2','SP7','ALPL','IBSP','SPP1','MGP','ENPP1','STC1'])][['state','gene','patients','median_log2CPM_difference','positive_patients','negative_patients']]
bv=bulk[['program','OR_per_SD','CI_low','CI_high','P_value','FDR_BH']]
lv=lung[lung.site.str.lower()=='lung'][['sample','program','cells','two_gene_positive_cells','two_gene_positive_pct']]
checklist=[
 ('1-2','公开数据和首版队列','完成','33个组织矩阵证据；GSE21257完整矩阵及48701行平台注释验收'),
 ('1-2','患者和建库映射','完成主表并保留未知','源患者ID、部位、治疗及建库缺项逐列保留；不合并跨队列N'),
 ('1-2','正常参考来源','完成签名替代路径','已发表肺泡/外膜标志及Gong人肺共表达面板；未加载完整HLCA分类器'),
 ('3-4','QC 双细胞 谱系 CNV','完成','沿用已验收阶段；保留失败/分歧/未评估；新增GEO身份报告另附'),
 ('3-4','环境RNA和参考敏感性','完成','发现阶段17份DecontX和6份参考敏感性；新增GEO12份DecontX'),
 ('5-6','CAF与myCAF细分','完成并采用计划内替代','H0/H2为探索性状态；myCAF细分缺跨患者支持，交付连续程序'),
 ('5-6','完整成骨矿化钙面板','完成','补齐16个成骨和15个钙调节基因；62个状态-基因效应'),
 ('7-8','独立单细胞验证','完成且保留验证失败','供者留出均值0.561未通过；外部状态标签不予接受，表达程序证据完整保留'),
 ('7-8','bulk转移关联','完成','53例、8个冻结程序、BH及年龄性别/基线转移/逐例剔除敏感性'),
 ('7-8','表面marker与抗体','完成候选设计','膜定位、RNA覆盖和误入、厂家/原始研究抗体证据'),
 ('9-10','稳健性和实验设计','完成','历史敏感性纳入；本轮供者留出、阈值及bulk敏感性；提供候选门控和验证'),
 ('增强项','A/B1与同步骨肺配对','资料依赖未执行','缺转移时间线、配对时间/治疗；不虚构分组'),
 ('增强项','后到ScPCA扩展','未纳入首版','按原计划第6天后大队列归增强版；不阻塞公开数据主线'),
 ('配置项','仪器通道和消化方案','待实验条件','先交付不绑定荧光通道的门控和抗体候选，不虚构配置')]
pd.DataFrame(checklist,columns=['阶段','工作项','状态','验收说明']).to_csv(O/'tables/original_plan_acceptance.tsv',sep='\t',index=False)
report=f"""# 骨肉瘤肺转移相关CAF原计划首版分析报告

本报告完成原定十天计划的公开数据主线、计划内替代分析和实验设计交付。研究目标是识别肺转移相关CAF候选状态，解析myCAF与成骨、矿化和钙调节程序，并提出可检验的表面marker组合。

## 核心结论

现有数据支持收缩相关和炎症相关的探索性CAF表达状态，但不支持将其命名为已稳定复现的肺转移特异新亚群。myCAF内7个细分簇均仅有1位患者达到每簇20细胞，独立供者留出验证平均平衡准确率0.561，低于预设0.60；本次交付连续程序与逐患者基因证据。

GSE21257的53例组织表达验证显示，抗原呈递程序与较低转移发生关联：每升高1 SD的OR为0.192，95% CI 0.074–0.493，P=0.000617，BH FDR=0.00494。该结果属于组织层面关联。H0/H2、成骨及钙调节程序未通过主分析FDR校正。

FAP与THY1可作为探索性富集组合的实验起点，但RNA筛查显示明显非CAF混入。本版不将任何组合称为已经验证的纯CAF或新亚群分选方案。

## Methods

沿用既有逐样本QC、双细胞、谱系、CNV与环境RNA审计。汇总33个原患者组织分析单元，保留原始样本、患者、取材、治疗和捕获映射字段。GSE152048的2047个严格候选来自6位患者；新增GEO严格集合为2065个细胞，分队列使用，154个阈值敏感候选单列。原始计数、身份阈值和既有分群均不修改。

发现阶段以已复核Harmony H0/H2为探索性候选状态。补齐客户16个成骨和15个钙调节基因。每位患者状态内与其余CAF各至少20细胞时，分别从原始counts构建pseudobulk，计算log2(CPM+1)差值；跨患者报告中位效应和方向。一项比较至少3位患者时，计算精确符号翻转P值，按每状态的31个客户基因家族进行BH校正。S100A8/S100A9作为髓系相关控制单列。

冻结发现队列H0/H2各30个候选基因及预定义程序后开展验证。状态映射使用跨矩阵可测且不属于后续膜标志/客户验证基因的49个特征，采用标准化和L2多分类logistic模型，按供者平衡样本权重并平衡类别。以逐供者留出平均平衡准确率≥0.60作为接收门槛，外部细胞还须预测置信度≥0.70。本轮同时保留原始及DecontX计数的逐样本程序值、覆盖和多基因检出比例。

肺参考采用原计划允许的公开标志组合路径。肺泡与外膜参考分别使用5和3个已发表标志；Gong参考使用论文明确列出的人肺PTGS2及IL6、CCL2、CXCL1、PTX3、TNFAIP6共表达基因。报告至少2个参考基因检出的细胞比例及连续表达。各参考可重叠。[肺参考来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC7704697/)；[Gong等原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC9830653/)。

GSE21257使用官方53列表达矩阵及GPL10295的48701行探针注释。沿用作者lumi方差稳定变换和robust spline normalization值，不再次取log或进行RNA计数标准化；同基因多探针取中位数，得到24998个基因。程序先逐基因z标准化后求均值，再标准化为1 SD。主终点沿用原记录34例转移、19例未转移；8个冻结程序分别做logistic回归，报告OR、95% CI、P及同族BH FDR。敏感性分析包含年龄性别调整、剔除诊断时已有转移者、随访记录核查和逐患者剔除。[GEO数据来源](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21257)。

对冻结的16个膜候选及其两两组合进行RNA共检出筛查，使用PTPRC、PECAM1、EPCAM作为排除信号，分别考察原始计数≥1及≥2的口径。记录目标覆盖、已知CAF占比、CNV高间质细胞和未决细胞混入。按供者中位数排序，结合UniProt蛋白定位及厂家/原始研究的抗体证据设计实验原型。

## Results

### 身份甄别形成发现与独立表达验证集合

33个组织的细胞顺序和原始计数核对通过。严格CAF主要分析对象如下，零候选组织完整保留在样本主表。GSE250015作为患者重叠未解决的补充表达证据单列。

{tab(sampleview)}

### myCAF细分主要受患者来源影响

{tab(myreview)}

所有细分簇仅1位患者达到20细胞门槛。按照原计划的后备规则，当前数据不支持命名稳定的新myCAF亚群，保留较高层级状态与连续程序。发现阶段H0仍显示RUNX2和ALPL偏高的探索性方向，但完整成骨/钙调节程序并非所有基因一致变化。

{tab(key)}

完整62项效应、P和FDR在full_customer_gene_effects.tsv；逐患者结果在full_customer_patient_effects.tsv。小样本的方向和效应不能因离散检验FDR阴性而被抹去，也不能提升为已验证机制。

![完整客户基因效应](../figures/customer_gene_effects.png)

### 独立验证未通过状态迁移门槛

六次供者留出的平均平衡准确率为0.561，单个供者为0.408–0.637。未达到预设门槛后，7份新增严格CAF样本的主状态映射全部保持未分配；诊断性预测单独保存。本次没有以验证队列重新挑基因、调门槛或强制分群。独立队列的逐样本连续程序与marker表达仍在完整统计表中，用于评估表达支持和技术差异。

![供者留出检验](../figures/transfer_validation.png)

### 组织抗原呈递程序与转移终点存在校正后关联

{tab(bv)}

抗原呈递程序在年龄性别调整后OR=0.194；排除14例诊断时已有转移者后，39例中的OR=0.234，95% CI 0.090–0.605。逐例剔除53次均保留负向关联。炎症程序仅名义显著（P=0.0210，FDR=0.0839）；其余程序未达到主分析校正后显著性。

![组织层面关联](../figures/bulk_association.png)

### 肺参考表型可量化但不构成来源或迁移证明

{tab(lv)}

这里的比例表示参考标志共检出，不是经正常肺参考分类器确定的驻留细胞比例。GSE152048和GSE318840各仅1份肺组织保留足量严格候选；骨肺比较按队列提供描述性差值，不建立A/B1或同步配对标签。

### 表面组合可进入预实验但尚不能可靠排除非CAF

{tab(top[['positive_genes','donors','median_RNA_target_coverage','median_RNA_known_CAF_fraction_gate','worst_RNA_CNV_high_fraction_gate']])}

FAP+THY1在发现队列5位具备至少20个H0细胞的供者中，RNA目标覆盖中位数为39.1%，RNA门内已知严格CAF占比中位数为19.1%，个别样本CNV高间质细胞混入可达16.0%。这些数字是RNA筛查指标，不是蛋白分选纯度。未决细胞完整保留，不能当作已证实CAF。

![候选组合筛查](../figures/marker_pairs.png)

## 流式分选和分选后验证设计

主预实验按单细胞、活细胞、CD45阴性/CD31阴性/EpCAM阴性、FAP与CD90二维分区的顺序保留各象限。FAP+CD90+作为候选富集池，与单阳和双阴对照共同分选，不直接命名为myCAF或肺转移特异CAF。另设CD140a比较支路，检验是否改善富集；不依赖CD45/CD31排除门来排除恶性骨肉瘤间质细胞。

FAP可评估427819克隆，厂家有人细胞流式验证；CD90可评估5E10克隆，厂家提供人细胞表面流式质控。CD140a的alphaR1用于已发表的人骨髓基质细胞研究。具体货号的偶联形式、保存液及批次需匹配活细胞用途和实际仪器，不将AbSeq或胞内试剂直接替代流式试剂。[FAP证据](https://www.rndsystems.com/products/human-fibroblast-activation-protein-alpha-fap-apc-conjugated-antibody-427819_fab3715a)；[CD90证据](https://www.biolegend.com/en-gb/products/brilliant-violet-711-anti-human-cd90-thy1-antibody-15886)；[CD140a原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC8252541/)。

预实验记录消化前后回收、活率和表位保留，并设置未染、单染、补偿/解混及关键marker的FMO对照。保留未分选组织、排除门阳性细胞和各阳性象限作身份对照。分选后用独立DNA/CNV或与原肿瘤异常谱对照的检测评估恶性混入，结合DCN/LUM/PDGFRA等成纤维表达，以及ACTA2/TAGLN和完整成骨/钙调节程序复核。涉及培养时记录培养前后状态漂移。

ACTA2、TAGLN及S100A4作为胞内/表达表型验证对象，不作为本版活细胞表面分选门。SR1004的问询另有独立编号核查，不将其与S100A4混用。

## Discussion

本研究完成现有公开数据能支持的身份、状态、独立验证、组织结局和实验设计主线。核心阴性结果是离散myCAF亚群缺乏跨患者支持、状态映射未通过供者检验，以及候选RNA分选组合的非CAF混入。它们决定本版应交付可验证的程序和富集原型，而不是已确认新亚群或可直接使用的纯化方案。

现有单细胞有效独立样本有限；取材、治疗、平台和临床分组不能完全分离。正常肺结果属于公开标志参考，未加载完整HLCA对象，也不能证明骨CAF向肺迁移。跨队列患者去重、A/B1时间线和同步配对条件尚未建立，不把细胞数或GSM数作为额外独立患者N。

GSE21257反映组织表达，抗原呈递信号可来自免疫等多类细胞，不能归因为apCAF，也不能等同CAF比例或因果保护作用。其转移终点未细化为肺转移；主分析包含基线已有转移者，已按既定敏感性路径剔除核查。

DecontX校正和计数阈值会影响低表达标志，原始、校正和不同检出规则均保留。RNA丰度不等于膜蛋白可及性；厂家流式证据不能替代本项目新鲜骨肉瘤组织上的消化、抗体、门控和分选后验证。

## 原计划验收与资料依赖

{tab(pd.DataFrame(checklist,columns=['阶段','工作项','状态','说明']))}

原计划的可执行主线按实际结果完成，资料依赖和增强项保持开放，不将未完成项目改写成已完成。细胞通讯、转录因子与拟时序是原计划按生物问题择优的可选项，本轮未在缺少稳定独立亚群的情况下额外扩展。

## 交付内容

样本与来源主表、逐细胞身份索引、完整客户基因效应、独立验证及失败记录、bulk临床映射和八程序关联、敏感性结果、肺参考表型、膜定位与抗体证据、候选RNA门控组合、方法代码和图件均随交付包提供。原始矩阵与较大中间对象保留在服务器版本目录。

本版分析不采用后续CD36/MYH11/C7/FNDC1参考论文分群替换原计划。完成本版后，另行评估其方法可用性。
"""
(O/'delivery/original_plan_report.md').write_text(report,encoding='utf-8')
# Supplementary linked paths are made self-contained in the package.
for src in [D/'analysis/state_top40_effect_markers.tsv',D/'analysis/mycaf_top30_effect_markers.tsv',D/'analysis/mycaf_subcluster_support_programs.tsv',D/'integration_review/validation/conclusions.json',D/'remediation_20260917/STATUS.json']:
 shutil.copy2(src,O/'delivery'/('prior_'+src.parent.name+'_'+src.name))
(O/'CORE_REPORT_READY').write_text('CORE_ANALYSES_AND_REPORT_READY_VISUAL_QA_PENDING\n')
print('CORE_REPORT_READY',len(report),flush=True)
