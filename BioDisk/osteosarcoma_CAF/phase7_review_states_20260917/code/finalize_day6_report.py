from pathlib import Path
import json,hashlib,datetime
import pandas as pd,numpy as np
from PIL import Image
p=Path.cwd();a=p/'analysis';d=p/'delivery'
def read(n):return pd.read_csv(a/n,sep='\t')
def table(df):
 def v(x):
  if pd.isna(x):return 'NA'
  if isinstance(x,(float,np.floating)):return f'{x:.3f}'
  return str(x).replace('|','／').replace('\n',' ')
 return '\n'.join(['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(v(x) for x in row)+' |' for row in df.itertuples(index=False,name=None)])
s=read('state_support_and_programs.tsv');my=read('mycaf_subcluster_support_programs.tsv');summary=json.loads((a/'analysis_summary.json').read_text())
assert (p/'logs/day5_6.exit_code').read_text().strip()=='0'
assert s.cells.sum()==2047 and my.cells.sum()==1287
stats=[]
for fn,key,n in [('state_gene_effects_all.tsv.gz','state','eligible_paired_patients'),('mycaf_gene_effects_all.tsv.gz','mycaf_state','paired_patients')]:
 x=read(fn);assert x[['P_value','FDR_BH']].isna().all().all();assert x[n].max()<3
 for state,z in x.groupby(key):stats.append({'analysis':fn,'state':state,'eligible_paired_patients':int(z[n].iloc[0]),'gene_rows':len(z),'P_FDR_status':'not tested; insufficient independent pairs'})
pd.DataFrame(stats).to_csv(d/'差异分析有效患者与统计口径.tsv',sep='\t',index=False)
labels=['收缩与基质程序，BC10主导','基质及外膜参考标志，BC10主导','收缩伴细胞周期程序，BC2主导','收缩与基质程序，BC2主导','基质重塑程序，BC21主导','基质重塑及外膜参考标志，BC22主导','炎症与应激程序，BC21单患者','炎症与抗原呈递候选程序，BC5单患者']
ann=s[['state','cells','patients_at_least_10_cells','max_patient_fraction','dominant_patient','cross_patient_support']].copy();ann['interpretation']=labels;ann['accepted_general_CAF_subtype']=False
ann.to_csv(d/'CAF计算簇生物学复核.tsv',sep='\t',index=False)
my['accepted_general_myCAF_subtype']=False;my.to_csv(d/'myCAF细分复核.tsv',sep='\t',index=False)
figs=[]
for f in sorted((p/'figures').glob('*.png')):
 with Image.open(f) as im:sz=im.size;im.verify()
 figs.append({'file':f.name,'width':sz[0],'height':sz[1],'file_integrity':'PASS','visual_layout_review':'PENDING'})
assert len(figs)==21
pd.DataFrame(figs).to_csv(d/'图件完整性与视觉验收状态.tsv',sep='\t',index=False)
meta=read('patient_metadata.tsv');genes=read('CAF_vs_malignant_gene_specificity.tsv');sens=read('extended_set_sensitivity.tsv');qc=read('quality_sensitivity_counts.tsv')
sel=genes[genes.gene.isin(['DCN','LUM','COL1A1','MGP','ENPP1','ALPL','SP7','IBSP','BGLAP','STC1'])][['gene','paired_patients','median_log2CPM_difference','patients_CAF_higher']]
med=read('GSE152048_strict_sample_program_medians.tsv');lung=med[['patient_id','program_Lung_alveolar_reference','program_Lung_adventitial_reference']]
report='''# 骨肉瘤 CAF 第1至6天中期报告

日期：2026年9月17日。范围：已实际接管队列的质量控制、CAF身份复核、状态与程序分析，以及原计划第1至6天遗漏审计。第7至8天独立验证和临床结局分析尚未开展。

## 严格CAF集合已经冻结，表达程序比新增离散亚型更适合作为本阶段结论

完成17个新增组织矩阵的下载验收、逐组织QC、谱系与CNV分析。125,212个保留细胞中，2,047个通过严格CAF标准，来自GSE152048的6位患者；扩展集合共3,167个细胞，其中GSE152048为3,055个、GSE299023为112个。既有GSE270231的188个严格候选单独用于探索性投影。上述集合和统计单位不合并计算样本量。

本阶段显示CAF具有收缩、基质重塑及炎症等表达程序。8个计算簇和7个myCAF细分簇均不能确认为普遍的新亚型；分簇主要受患者组成影响。3例患者内CAF与恶性细胞比较中，MGP均在CAF更高，而ALPL、SP7、IBSP均在恶性细胞更高。当前数据不支持以这些结果确认稳定的“成骨型CAF”亚型。

## Methods

### 数据与临床

新增输入包含GSE299023的6个非PDX患者组织和GSE152048的11个组织；原始160,541个barcodes，经QC保留131,609个，排除6,397个双细胞后得到125,212个细胞。另使用已接管GSE270231严格CAF集合进行独立存储的探索性投影。原始counts、组织来源、治疗记录、病理类型、肿瘤尺寸、坏死率及Ki67按来源保存。统计重复单位为患者。所有捕获建库与患者映射不确定项均在主表显式标记。

### 身份与参考复核

CAF纳入采用多基因谱系、两个CNV方法一致低CNV、至少500个检测基因、至少2个成纤维标志及竞争谱系筛查。正常参考按预先冻结的免疫与内皮训练/留出角色评估，每个原始分层至少10个可评估细胞且正常比例至少90%。inferCNV为第一方法，CopyKAT或通过参考门槛的SCEVAN为第二方法。保留原始分歧、失败和待定分类。

### 状态与程序

对严格CAF的原始counts进行每细胞10,000总量标准化及log1p转换，采用患者感知HVG、PCA、近邻图及多分辨率Leiden。聚类HVG排除线粒体、核糖体、性染色体及预设应激/周期基因。预设收缩、炎症、抗原呈递、基质重塑、成骨、矿化/钙调节等程序，使用表达匹配对照评分，并保存基因覆盖、共检测和平均表达。跨患者支持门槛为至少3位患者各贡献10个细胞；同时检查最大患者占比。myCAF候选满足至少3个收缩基因检测及正评分，在至少150个细胞且3位患者各20个细胞时执行条件性细分。

### 患者层面统计与敏感性

保存患者原始counts pseudobulk和CAF内部状态比例。状态与其余CAF比较要求同一患者两侧各至少20个细胞；至少3个合格患者时预设精确配对符号翻转检验，按每个状态的全基因比较家族单独计算BH FDR，否则仅报告适用的表达效应。无配对患者的汇总表达差异不解释为患者内效应。CAF与恶性细胞比较使用同患者每类至少20个细胞的描述性log2(CPM+1)差异。另进行分辨率、提高基因检测门槛、去除高复杂度标记及严格/扩展集合敏感性核查。

正常肺比对采用文献来源的肺泡/外膜成纤维细胞候选标志组合。外部投影使用冻结的发现集特征及近邻距离筛查，不重新定义发现集状态。

## Results

### 补充CNV复核完成，五个重点样本仍未通过严格纳入门槛

BC6及另外两个追加SCEVAN任务均已正常结束。参考审计通过的组织仍为12/17，严格/扩展CAF数量未改变。OS.C.21存在第二方法失败或留出参考不足；OS.C.32的inferCNV基因数门槛失败；OS.C.33和BC3参考分层细胞不足；BC6两个第二方法均与冻结正常参考明显矛盾。因此未将这些样本强行纳入严格CAF。完整数值见“五个重点样本复核结论.tsv”和“工具参考检查完整结果.tsv”。

### 六位患者贡献严格CAF，肺病灶来自一位患者

'''+table(meta[['patient_id','site','CAF_cells']])+'''

共有5位骨病灶患者、1位肺病灶患者。其余组织严格CAF为零，不能据此推断组织生物学上没有CAF。状态比例的分母为该患者通过标准的严格CAF；该比例不是全肿瘤CAF浸润率。

### 计算簇主要由患者驱动，myCAF细分未形成跨患者亚型

'''+table(ann[['state','cells','patients_at_least_10_cells','max_patient_fraction','interpretation']])+'''

只有C1达到至少3位患者各10个细胞的最低计数门槛，但91.4%的细胞来自BC10，仍不能作为通用亚型。1,287个myCAF候选满足执行细分的条件，得到7个细分簇，全部未达到跨患者支持门槛。C2的细胞周期程序较高；C6应激程序较高。PC1与检测基因数及线粒体比例的描述性Spearman rho分别为0.563和0.649，提示技术质量与患者差异必须共同考虑。结果以程序及患者内描述为主，不强制划分iCAF、apCAF或新的myCAF亚型。

### 成骨标志与矿化调节基因呈现不同的细胞类型分布

下表为3位合格配对患者的描述性比较，正值表示CAF高于恶性细胞；效应为患者内log2(CPM+1)差异的中位数。

'''+table(sel)+'''

MGP在3/3位患者的CAF中更高，ALPL、SP7、IBSP、BGLAP在3/3位患者的CAF中更低。ENPP1方向为2/3位患者CAF更高。COL1A1和LUM在两类细胞中均广泛检测，表达富集不等于特异标志。上述结果将矿化调节候选线索与完整成骨程序区分开来。

### 肺外膜参考标志并非肺病灶专属

'''+table(lung)+'''

外膜参考组合在部分骨病灶和肺病灶CAF均有较高评分，肺泡参考组合未形成一致的高评分；当前数据不支持将这一比对解释为肺来源或骨肺迁移的证明。来源为公开的正常肺细胞标志研究：[Travaglini等](https://pmc.ncbi.nlm.nih.gov/articles/PMC7704697/)及[Tsukui等](https://pmc.ncbi.nlm.nih.gov/articles/PMC7174390/)。

### 临床组间差异仅能作为单病例与骨病例分布的描述

BC10相对5例骨病灶患者的中位程序评分差：收缩+0.304、基质重塑−0.933、成骨+0.010、矿化/钙调节+0.034。未进行肺与骨组间推断检验。虽然治疗和病理协变量设计矩阵满秩，单一肺患者仍不构成可估计的生物学组内重复。

### 完整差异表已输出，但有效配对患者不足以进行预设检验

状态差异表262,912行，myCAF差异表230,048行；各比较有效配对患者分别为0–2例和0–1例，均未达到预设推断门槛，所有P值/FDR保持NA。这代表未作推断检验，不是“校正后不显著”。有合格配对时报告方向和效应；没有合格配对时只保留描述性表达记录及明确标记。

### 集合敏感性显示基质程序较稳定，矿化及成骨排序较敏感

'''+table(sens)+'''

提高基因检测门槛并排除高复杂度标记后保留'''+str(int(qc['sum'].sum()))+'''个严格CAF。Leiden分辨率0.3/0.5及0.5/0.8的ARI分别为0.778、0.662。既有GSE270231的188个细胞完成探索性投影，特征覆盖99.9%，未触发距离拒绝；已确认患者只有1位，另80个细胞患者映射未确定，因此不构成独立生物学重复的验证成功。

## 第1至6天工作覆盖与补漏

已补录公开临床表遗漏字段，形成24个实际分析单元、89个GEO记录的四层主表框架；已输出17张最终谱系/身份图，完成周细胞与血管平滑肌表达证据分型、CAF与其他谱系比较、肺参考标志分析、患者比例、counts pseudobulk、程序与敏感性分析。完整25项审计见“原计划第1至6天工作覆盖审计.md”。

仍未闭环：原定GSE318840、GSE162454、GSE21257的完整输入可用性验收；部分捕获建库/患者映射；A/B1所需肺转移史和时间；部分ScPCA扩展信息；健康正常骨参考。GSE169396来自骨关节炎伴骨量减少患者，GSE233676来自骨折不愈合，均未当作健康正常骨。[骨来源论文](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2023.1107511/full)、[GSE233676记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233676)。

## Discussion

本阶段确立了经过双方法和谱系复核的CAF分析集合，并显示收缩/基质及矿化调节相关表达。患者驱动、单一肺患者、CNV参考适用性和集合敏感性限制了通用亚型与临床组间结论。缺少空液滴配套输入，未完成环境RNA污染率估计或校正；跨谱系筛查不能替代该步骤。apCAF和混合表达保持候选定位。原始参考或CNV失败不代表肿瘤组织不存在CAF；严格选择本身也会影响可观察程序。表达相关性、候选特异性和外部近邻投影均不证明因果、起源、蛋白定位或临床可用性。

后续应优先完成预定独立队列的实际输入验收与患者去重，核查CAF信息量，再验证患者层面程序。第7至8天分析尚未运行。本次按用户要求在第5至6天计算及数值复核后暂停扩展分析。

## 交付与验收状态

完整计算产物位于同级analysis目录；delivery包含身份复核、状态解释、统计口径和图件审计。21张PNG通过文件完整性检查，但尚未完成逐图视觉排版验收，本报告为可复核的文字与数值中期版本。既有自动审批未允许远端图件下载到本地，也未允许本机ScPCA映射表上传，未绕过这些限制。图件不得标记为已完成视觉QA。
'''
(d/'骨肉瘤CAF_第1至6天中期报告.md').write_text(report,encoding='utf-8')
review={'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'compute_exit_code':0,'strict_CAF':2047,'myCAF_candidates':1287,'general_new_CAF_subtypes_confirmed':0,'general_myCAF_subtypes_confirmed':0,'inferential_state_tests':0,'inferential_site_tests':0,'figure_integrity_pass':21,'visual_QA':'pending','status':'NUMERIC_REVIEW_COMPLETE_INTERIM_TEXT_REPORT_READY','remaining_early_plan_gaps':'see original_day1_6_checklist.tsv'}
(p/'review/day6_scientific_review.json').write_text(json.dumps(review,indent=2))
print(json.dumps(review))
