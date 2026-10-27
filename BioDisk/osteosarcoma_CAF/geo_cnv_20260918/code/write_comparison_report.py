from pathlib import Path
import pandas as pd,json,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918')
s=pd.read_csv(P/'review/input_summary.tsv',sep='\t');t=pd.read_csv(P/'review/initial_vs_reviewed_labels.tsv',sep='\t');f=pd.read_csv(P/'review/FLEX_author_code_marker_support.tsv',sep='\t')
ready=s[s.cnv_eligible];hold=s[~s.cnv_eligible]
def table(df):
 return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(str(x) for x in row)+' |' for row in df.itertuples(index=False,name=None)])
z=f[(f.current_label=='FibLike')&f['sample'].isin(['GSE318840__GSM9504303','GSE318840__GSM9504304'])&f.author_code_panel.isin(['OS','Fibroblast','Pericyte'])].copy();z.pct_two_detected=z.pct_two_detected.round(2)
text=f'''# GEO 注释对比与后台 CNV 启动记录

日期：2026-09-18。

## 15份样本具备CNV输入条件，1份因正常参考不足暂缓

完成了本项目初轮与复核注释的逐细胞对比，以及原作者公开论文、标志集合和GSE318840注释代码的对照。15份样本通过counts、基因坐标及正常参考检查，安排inferCNV与CopyKAT共30个后台任务。覆盖 {int(ready.run_cells.sum()):,} 个细胞；所有保留细胞都进入待评估集合，仅将筛选出的正常细胞分配为训练参考或留出校准细胞，不只分析FibLike群。

GSE318840__GSM9504299暂缓：严格筛选后只能分配18个训练正常细胞和10个留出细胞，未达到30/20的门槛。没有放宽参考阈值，也没有将其标为无CAF或CNV阴性。

## 原作者对照支持保留OB/Fib歧义，并用CNV继续判定

GSE318840原研究及注释代码把骨样/成纤维样表达重叠群作为OB/Fib处理，随后使用CNV判断恶性与低CNV信号。作者公开代码中的骨样、成纤维和周细胞标志集合在本项目肺样本候选中均有重叠支持。下表为cluster内至少检出该作者panel两个基因的细胞比例，仅作描述性核查：

{table(z[['sample','cluster','cells','author_code_panel','pct_two_detected']])}

因此，1,025和562个FibLike仍是待CNV的表达候选，不是作者原标签的逐细胞复现，也不是已确认CAF。上述重叠群全部作为CNV目标，未混入正常参考。

GSE162454按原文及勘误后的标志口径对照。原文CAF标志与骨样表达存在重叠；本项目保留骨样、成纤维样和未定群。GSE250015原文将MYF5/CHRNA1群称为osteocyte，本项目相关肌源性/歧义群仍保留为目标，未按名称将其视为正常骨细胞参考。病例转移状态与实际取材部位分开保留。

## 对比层级和数据范围

- 本项目逐细胞对比覆盖 {int(t.cells.sum()):,} 个最终保留细胞，完整转换表为 `initial_vs_reviewed_labels.tsv`。这是初轮自动标签对复核标签的对比，不是作者标签一致率。
- 在已下载GEO文件和已核查作者仓库中未获得完整的作者barcode→cell type表或注释RDS/H5AD。GSE318840仓库提供注释规则与分析代码，代码引用的RDS位于作者本地路径。作者整合后的cluster编号不能映射成本项目逐样本cluster编号。
- 原作者对照使用论文标志、作者代码panel及细胞群定义，记录完整panel覆盖和表达支持；未虚构逐细胞混淆矩阵或一致率。
- 原文处理后的细胞数与本项目最终保留量不同，不能当成注释错误率：GSE318840原文40,401，本项目46,935；GSE162454原文29,278，本项目38,424；GSE250015公开输入18,179，本项目再次QC后16,128。这里记录处理口径差异，不声称复现作者全部过滤流程。

## CNV输入及运行设置

{table(s[['sample','run_cells','ref_immune','ref_endo','holdout','cnv_eligible']])}

正常参考仅来自同一样本的标志支持免疫/内皮群，沿用项目的基因数、UMI、竞争谱系表达及双细胞评分筛选；FibLike、骨样、周细胞、破骨、肌源及未定群均不能作参考。每份至少30个训练正常细胞、20个独立留出校准细胞；禁止不足时退回训练细胞校准。基因映射保留全counts，不能映射的symbol单独标识，由工具按坐标筛选。15份导出矩阵与accepted对象的counts逐元素一致，细胞顺序、基因顺序和每细胞总counts均核对通过。

使用项目既有inferCNV/CopyKAT方法参数，inferCNV cutoff=0.1、101基因平滑窗、denoise开启、HMM关闭；CopyKAT窗口25、KS.cut=0.1、同样本显式正常参考。方法参数不等于对作者所有参数的原样复现。环境预检确认infercnv 1.26.0、copykat 1.2.5，15份样本的CopyKAT基因映射满足输入要求。

后台同时运行两条方法队列，每个任务4线程、64 GiB虚拟内存上限，优先肺样本。队列原生等待子进程结束，并检查退出码、COMPLETE标记、参数和逐细胞结果完整性；失败任务单独记录。完成后仅产生CNV计算结果与完整性摘要，不自动确认CAF或增加独立患者N。FLEX适用性、留出正常信号、跨方法一致性及环境RNA边界需随结果复核。

## 文件与恢复入口

- `STATUS.json`：当前后台状态；`logs/cohort.pid`：队列PID。
- `logs/cohort.log`、`logs/infercnv_jobs.json`、`logs/copykat_jobs.json`：方法队列与任务记录。
- `review/annotation_gate.json`：本轮启动依据及暂缓样本。
- `review/author_marker_support.tsv`、`review/FLEX_author_code_marker_support.tsv`：作者标志对照。
- `review/export_integrity.json`：新CNV输入的counts一致性记录。
- `logs/cohort.exit_code`及`JOBS_FINISHED`：最终退出与完成记录；进程消失不代表成功。

## 来源

[GSE318840原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13330631/)；[作者注释代码](https://github.com/khanhnguyendok/osteosarcoma/tree/4305dae30e226ea46bbeccf05b5e7ebdeef91076/06_major_celltype_annotation)；[GSE162454原研究](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full)；[GSE162454勘误](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.1077067/full)；[GSE250015原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11467127/)。作者代码固定于提交4305dae30e226ea46bbeccf05b5e7ebdeef91076，下载文件SHA256已保存。
'''
(P/'delivery/GEO_annotation_comparison_CNV_launch_20260918.md').write_text(text)
import shutil
for name in ['input_summary.tsv','initial_vs_reviewed_labels.tsv','author_marker_support.tsv','FLEX_author_code_marker_support.tsv','author_marker_coverage.tsv','annotation_gate.json','export_integrity.json','CNV_environment.json','copykat_mapping_preflight.tsv']:
 shutil.copy2(P/'review'/name,P/'delivery'/name)
print('ANNOTATION_COMPARISON_REPORT_COMPLETE')
