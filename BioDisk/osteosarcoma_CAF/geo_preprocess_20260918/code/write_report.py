from pathlib import Path
import json,datetime,html,hashlib,zipfile
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918')
s=pd.read_csv(P/'delivery/GEO_sample_acceptance_and_QC.tsv',sep='\t')
p=pd.read_csv(P/'delivery/GEO_cluster_lineage_review.tsv',sep='\t').fillna('')
assert len(s)==16
q=json.loads((P/'review/visual_QA.json').read_text());assert q['status']=='PASS'
for n in ['pipeline','review','doublet','accepted']:
 assert (P/'logs'/f'{n}.exit_code').read_text().strip()=='0'
agg=s.groupby('dataset',sort=False).agg(samples=('sample','size'),input_cells=('input_cells','sum'),initial_retained=('retained','sum'),analysis_ready_cells=('analysis_ready_cells','sum'),initial_doublets=('doublets','sum'),extra_doublets=('extra_doublets_excluded','sum'),fibroblast_like=('fibroblast_like_after_doublet_review','sum'))
def table(df):
 cols=list(df.columns);out=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for row in df.itertuples(index=False,name=None):out.append('| '+' | '.join(str(v).replace('|','/') for v in row)+' |')
 return '\n'.join(out)
lung=s[s.site=='lung'];markers=int(s.proposed_fibroblast_cells.sum())
text=f'''# 骨肉瘤 CAF：新增 GEO 队列验收与预处理

处理日期：2026-09-18。目录：`{P}`。

## 16份矩阵完成验收与逐样本预处理

三个 GEO 队列共 44 个下载文件全部通过 SHA256 和文件大小核对，矩阵维度、非负整数 counts、唯一基因标识和细胞条形码检查通过。共读入 {int(s.input_cells.sum()):,} 个细胞，初轮逐样本 QC 与双细胞评估后保留 {int(s.retained.sum()):,} 个细胞；8份Scrublet阈值复核样本再经scDblFinder评估排除 {int(s.extra_doublets_excluded.sum()):,} 个预测双细胞，最终形成 {int(s.analysis_ready_cells.sum()):,} 个细胞的分析输入。所有保留细胞的 counts 已与输入逐项比较一致；16份输出均通过细胞数量、元数据与有限值检查，计算、初轮自动复核、补充双细胞评估和最终对象组装退出码均为0；最终子集counts与初轮对象一致。

{table(agg.reset_index())}

表中 fibroblast_like 为标志基因复核提出的成纤维样细胞数，不能视为已确认 CAF 数。自动谱系提出的完整证据保留在 `GEO_cluster_lineage_review.tsv`。本轮没有执行新队列的 CNV、环境 RNA 校正或跨队列状态整合。

## 两份肺取材样本提供新的候选输入

{table(lung[['sample','input_cells','analysis_ready_cells','fibroblast_like_after_doublet_review','immune_reference_after_doublet_review','endothelial_reference_after_doublet_review']])}

上述参考细胞数是当前标志与QC规则下的可选参考数量，尚未经过CNV结果和留出参考校准。GSE318840 的两份肺组织对应原研究报告的两位患者；有效 CAF 患者数仍须根据身份验收确定。GSE250015 的 Metastasis 是患者转移状态，两份材料均来自原发灶，未标作肺组织。

## Methods

数据来自已下载的 GSE318840（8份10x FLEX）、GSE162454（6份10x 3′ v3）与 GSE250015（2份Singleron GEXSCOPE）矩阵。按样本建立包含原始counts的AnnData对象，Ensembl标识与基因symbol分列保存；Singleron原始symbol标识保留。标准化为每细胞10,000 counts后进行log1p转换，原始counts保存为独立layer。

沿用既有项目逐样本阈值：log10(genes+1)与log10(counts+1)低于中位数减3个未缩放MAD，分别以200 genes和300 counts为下限；线粒体比例上限为中位数加3个MAD，限制于5%–20%。高复杂度只标记，不直接删除；同时记录不使用线粒体过滤时的细胞保留情况。FLEX的实测线粒体基因覆盖和表达单独记录，采用逐样本阈值。

Scrublet使用0.06的既有预期双细胞率先验及固定随机种子270231。自动阈值有限、QC后至少100个细胞且预测比例在0.1%–25%时排除预测双细胞，其余情况保留并标记复核。高变基因上限3,000，嵌入中排除线粒体、核糖体和指定性染色体技术驱动基因；PCA、邻居图、Leiden 0.4/0.8/1.2与UMAP逐样本计算，主要分辨率0.8。

谱系标志使用已记录的广泛细胞类型panel。cluster的前50个描述性标志按正向log2倍数变化乘检出率差排序，要求cluster内检出率至少10%；未计算细胞级P值或将细胞数当作生物学重复。谱系复核保存panel覆盖、检出比例和表达证据。原始表达输出与谱系复核表分开保存。软件版本、随机种子和脚本溯源见results与review目录。

8份Scrublet触发低检出率复核的样本采用scDblFinder 1.24.10补充评估：各样本独立运行，随机人工双细胞、dbr=0.06、dbr.sd=1、固定种子270231及串行BiocParallel。宽dbr.sd设置使阈值由分类误差决定，不强制删除6%的细胞。补充模型输出及原始阈值均保留，最终对象排除其预测双细胞。其余8份保留通过既定阈值检查的Scrublet处理结果。完整方法分配及前后数量见逐样本表。

标志图的点面积表示该cluster中检出基因的细胞比例，颜色表示该基因的cluster平均log1p表达在cluster间标准化后的z值。最终对象保存复核后的候选谱系标签。展示使用初轮UMAP坐标上的最终保留细胞子集；未在双细胞排除后重算邻居图或重新定义cluster，最终对象移除旧邻居图并明确标记初轮坐标。标志图在最终保留细胞上重新汇总。

## 本轮处理边界与下一步

1. 新队列的成纤维样、间充质歧义和骨样细胞保留为待鉴定对象；CAF与恶性身份仍需双方法CNV、参考/留出校准及必要的环境RNA敏感性验收。
2. GSE318840为探针型FLEX，其基因覆盖记录不等同于CNV方法已通过验证；不同技术的线粒体比例不作为生物学效应比较。
3. GSE162454原研究明确采样前未接受化疗或放疗；转移史仍不按推测分组。GSE250015与既有上海队列的患者交叉去重尚未解决，不计入合并后新增独立患者N。各队列内部样本ID和跨队列身份状态均保留。
4. 双细胞率0.06是分析先验，不是实测建库上样率；缺少独立capture映射时保留相应字段。未提供空液滴矩阵，未宣称完成环境RNA校正。
5. ScPCA由客户下载，本轮未接入。原有第1至6天中期报告与既有CAF身份集合保持原版本。

## 交付文件

- `GEO_sample_acceptance_and_QC.tsv`：逐样本输入、QC、双细胞、保留量及谱系候选数量。
- `GEO_cluster_lineage_review.tsv`：每个cluster的候选标签、panel证据与前列标志。
- `../accepted/<sample>/analysis_ready.h5ad`：最终保留细胞的标准化表达、counts layer、双细胞处理记录、复核候选谱系及初轮坐标。初轮完整对象另存于results目录。
- `../results/<sample>/all_cells_qc.tsv.gz`：所有输入细胞的保留/排除原因。
- `../review/`：文件完整性、counts一致性、panel覆盖、参考可用性及图形QA记录。

## 来源

[GSE318840原研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13330631/)；[GSE162454原研究](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full)；[GSE162454图形标签勘误](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.1077067/full)；[GSE250015官方元数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250015)；[scDblFinder官方方法文档](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html)。本轮采用counts重新注释，未沿用原文勘误前的CAF/骨样细胞图形标签。
'''
(P/'delivery/GEO_preprocessing_report_20260918.md').write_text(text)
htmltext='<!doctype html><meta charset="utf-8"><title>GEO QC review</title><style>body{font:16px Arial;margin:30px}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:6px}img{max-width:100%}</style><h1>GEO sample QC and provisional lineage review</h1>'+'<p>QC panels show the initial filtering. UMAP uses initial coordinates with the final retained subset and marker-reviewed candidate labels. CAF identity remains pending CNV.</p>'+s[['sample','input_cells','analysis_ready_cells','fibroblast_like_after_doublet_review']].to_html(index=False)
for sm in s['sample']:
 htmltext+='<h2>'+html.escape(sm)+'</h2>'
 for fn in ['QC_initial.png','UMAP.png','Markers.png']+(['Doublet_review.png'] if (P/'accepted'/sm/'Doublet_review.png').exists() else []):htmltext+=f'<img src="figures/{sm}/{fn}"><br>'
(P/'delivery/GEO_QC_gallery.html').write_text(htmltext)
import shutil
for sm in s['sample']:
 dst=P/'delivery/figures'/sm;dst.mkdir(parents=True,exist_ok=True)
 for fn in ['QC_initial.png','UMAP.png','Markers.png']+(['Doublet_review.png'] if (P/'accepted'/sm/'Doublet_review.png').exists() else []):shutil.copy2(P/'accepted'/sm/fn,dst/fn)
with zipfile.ZipFile(P/'GEO_preprocessing_delivery_20260918.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in (P/'delivery').rglob('*'):
  if f.is_file():z.write(f,f.relative_to(P/'delivery'))
print('GEO_REPORT_READY',flush=True)
