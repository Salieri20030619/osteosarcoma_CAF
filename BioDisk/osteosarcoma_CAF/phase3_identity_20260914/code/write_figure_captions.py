"""SYNTAX_CHECK: python -m py_compile write_figure_captions.py
RUN_COMMAND: python write_figure_captions.py <phase3>
EXPECTED_SUCCESS_MARKER: FIGURE_CAPTIONS_COMPLETE
"""
import sys
from pathlib import Path
import pandas as pd
p=Path(sys.argv[1]);d=p/'identity'
s=pd.read_csv(d/'sample_identity_summary.tsv',sep='\t')
sel=pd.read_csv(d/'method_selection.tsv',sep='\t').set_index('short')
qc=pd.read_csv(d/'tool_reference_QC.tsv',sep='\t')
plots=pd.read_csv(d/'CNV_plot_manifest.tsv',sep='\t')
genes=pd.read_csv(d/'identity_marker_expression.tsv',sep='\t');rows=[]
for r in s.itertuples():
    short=r.short
    rows.append({'file':f'{short}_identity_UMAP.png','caption':f'{short} 的{r.retained}个细胞/核。横纵轴为复用的UMAP1/UMAP2，左图按整合身份着色，中图为inferCNV，右图为{sel.loc[short,"selected_source"]}。本样本严格CAF候选{r.caf_strict}个，扩展候选{r.caf_extended}个，身份待定{int(r.Identity_pending)}个。每点为一个细胞/核，不代表独立患者。FibCand为身份候选类别，严格集合另见CAF_sets图。'})
    rows.append({'file':f'{short}_CAF_sets.png','caption':f'{short} 的CAF集合位置。横纵轴为UMAP1/UMAP2，灰色为全部{r.retained}个保留细胞/核，绿色分别标记严格候选{r.caf_strict}个和扩展候选{r.caf_extended}个；扩展集合包含严格集合。统计为本样本细胞数量。'})
    fib=genes[(genes.short==short)&(genes.identity=='FibCand')].set_index('gene')
    if len(fib):
        result=f'FibCand组有{int(fib.cells.iloc[0])}个细胞/核；'+', '.join(f'{g}检测比例{fib.loc[g,"pct_expressed"]:.1f}%' for g in ['DCN','LUM','COL1A1','ACTA2','ALPL','CD36'])+'。'
    else:result='本样本没有整合为FibCand的细胞，图中保留其他实际存在的身份组。'
    rows.append({'file':f'{short}_identity_markers.png','caption':f'{short} 的谱系及功能marker证据。横轴为基因，纵轴为实际存在的整合身份组，括号是该组细胞/核数；点面积表示检测比例，颜色表示平均log1p标准化表达。'+result+'这些数值是细胞级描述，不是患者组间检验。'})
for r in plots.itertuples():
    x=pd.read_csv(p/r.method/r.short/'cell_results.tsv',sep='\t')
    counts='；'.join(f'{k}={v}' for k,v in x.call.value_counts().items())
    q=qc[(qc.short==r.short)&(qc.method==r.method)].iloc[0]
    rows.append({'file':f'{r.short}_{r.method}_CNV.png','caption':f'{r.short} 的{r.method} CNV表达偏离图。横轴为1–22号常染色体，纵轴为{r.displayed_cells}个展示细胞/核；每个工具实际角色与原簇最多抽取25个可评估细胞，完整计算输入为{len(x)}个。Ref/Hold/Check/Target为基线参考、留出参考、未用于本版本基线的检查细胞及目标；Type为谱系。全体输入的工具结果为{counts}。参考质量检查'+('通过。' if q.reference_check_pass else '未通过，该版本不作为严格集合的独立支持证据。')})
pd.DataFrame(rows).to_csv(d/'figure_captions.tsv',sep='\t',index=False)
print('FIGURE_CAPTIONS_COMPLETE',len(rows))
