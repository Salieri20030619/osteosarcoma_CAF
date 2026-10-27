from pathlib import Path
import json,datetime,hashlib,zipfile
import pandas as pd,numpy as np,anndata as ad
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase67_DSC1_all_cells_20261005'
assert (O/'ANALYSIS_COMPLETE').exists() and (O/'logs/analysis.exit_code').read_text().strip()=='0'
r=json.loads((O/'resources/SUMMARY.json').read_text())
d=pd.read_csv(O/'tables/DSC1_all_cells.tsv.gz',sep='\t',low_memory=False)
pos=pd.read_csv(O/'tables/DSC1_positive_cells_with_native_QC.tsv.gz',sep='\t',low_memory=False)
# Direct original-H5AD row/column checks, using a read path independent of the streaming extraction.
checks=[]
for ds,g in d.groupby('dataset'):
 t=pd.concat([g[g.DSC1_assayed].nlargest(2,'DSC1_raw'),g[g.DSC1_raw.eq(0)].head(1)]).drop_duplicates('original_cell_id')
 for row in t.itertuples(index=False):
  f=B/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad'
  a=ad.read_h5ad(f,backed='r');ix=a.obs_names.get_loc(row.original_cell_id);j=a.var_names.get_loc('ENSG00000134765')
  v=float(a.X[ix:ix+1,j:j+1].toarray()[0,0]);assert v==row.DSC1_raw
  a.file.close();checks.append({'cell_id':row.original_cell_id,'source':ds,'raw_count':v,'matched':True})
pd.DataFrame(checks).to_csv(O/'QA/direct_count_checks.tsv',sep='\t',index=False)
native=[]
for c in pos:
 if c.startswith('native_'):
  for v,n in pos[c].fillna('not_available').astype(str).value_counts().items():native.append({'field':c,'value':v,'positive_cells':int(n)})
pd.DataFrame(native).to_csv(O/'tables/DSC1_positive_native_QC_summary.tsv',sep='\t',index=False)
sources=pd.read_csv(O/'tables/DSC1_by_source.tsv',sep='\t')
lineages=pd.read_csv(O/'tables/DSC1_by_lineage.tsv',sep='\t')
pool=pd.read_csv(O/'tables/DSC1_by_candidate_pool.tsv',sep='\t')
samples=pd.read_csv(O/'tables/DSC1_by_sample.tsv',sep='\t')
dist=pd.read_csv(O/'tables/DSC1_positive_count_distribution.tsv',sep='\t')
tot=r['overall'];outside=r['outside_candidate8174']
lines=['DSC1（SR1004抗体对应靶点）全细胞表达核查','2026-10-05 | Phase 67','',
'结论',
f"全细胞输入共{tot['cells']:,}个、107份样本、9个来源；DSC1有覆盖{tot['assayed_cells']:,}个、缺测{tot['missing_cells']:,}个。原始计数阳性{tot['positive_cells']:,}个，检出率{tot['positive_pct']:.6f}%；总原始计数{tot['raw_total']:.0f}，单细胞最大{tot['max_raw']:.0f}，阳性细胞中位数{tot['positive_median_raw']:.2f}。",
f"原8174候选池仍为3个阳性；池外{outside['positive_cells']}个阳性、检出率{outside['positive_pct']:.6f}%。候选池外阳性占全体阳性{100*outside['positive_cells']/tot['positive_cells'] if tot['positive_cells'] else 0:.2f}%。",
f"原始计数>=2的细胞{tot['ge2_cells']}个，>=5的细胞{tot['ge5_cells']}个；这些阈值仅描述计数分布，不定义生物学高表达。",
f"覆盖细胞平均原始计数{tot['mean_raw_all_assayed']:.8f}；全库CP10k逐细胞均值{tot['mean_CP10k_all_assayed']:.8f}；汇总计数/汇总全库计数的CP10k为{tot['pseudobulk_CP10k']:.8f}。跨来源汇总仅描述本数据集合，不作患者级推断。",'',
'数据与方法',
'使用107个样本的full_called_counts.h5ad原始X，按唯一DSC1/ENSG00000134765映射逐块读取，完整重算每个细胞的全库计数与检出feature数。缺少该feature标记缺测，未转成零。未使用HVG、整合矩阵或仅marker列表。',
'细胞类型继承phase55全细胞审阅标签与置信度；原8174候选池及双联体风险由phase66当前清单识别。全部source-called输入均保留，其中包含提供方标为Remove/质控不通过的细胞；不能把876648称为全部合格singlet。',
'仅报告描述性表达与检出指标，不把细胞当独立生物学重复，不进行差异检验。样本及来源内患者/部位汇总另表保存，跨来源患者重叠未解决，不合并作为独立样本量。','',
'按数据来源',
sources[['dataset','cells','assayed_cells','positive_cells','positive_pct','raw_total','max_raw','positive_samples','positive_donor_records']].to_string(index=False),'',
'按审阅谱系',
lineages.sort_values('positive_cells',ascending=False)[['review_lineage_code','review_expression_lineage','assayed_cells','positive_cells','positive_pct','raw_total','max_raw']].to_string(index=False),'',
'原候选池内外',
pool[['candidate_pool','assayed_cells','positive_cells','positive_pct','raw_total','max_raw']].to_string(index=False),'',
'原始计数分布',dist.to_string(index=False),'',
'阳性细胞最多的样本（完整表另附）',
samples.sort_values('positive_cells',ascending=False).head(15)[['dataset','sample','site','technology','assayed_cells','positive_cells','positive_pct','raw_total','max_raw']].to_string(index=False),'',
'原始计数最高的细胞（完整阳性表另附）',
pos.sort_values('DSC1_raw',ascending=False).head(15)[['original_cell_id','dataset','sample','review_expression_lineage','lineage_confidence','candidate_pool','DSC1_raw','DSC1_CP10k','full_library','detected_features']].to_string(index=False),'',
'解释与湿实验对照',
'本次回答全细胞范围内的DSC1 RNA检出及来源分布。客户报告的SR1004湿实验高信号作为待对齐证据保留；厂家货号、方法、组织/培养/分选细胞范围尚未取得，不能据本结果判定蛋白不存在或客户实验错误。',
'SR1004是抗体克隆号；已核对产品NBP3-22640靶向Desmocollin-1/DSC1。原始RNA计数不等同蛋白含量，谱系注释不等同蛋白定位或病理确诊。所有阳性细胞的原生双联体/提供方QC字段在附表保留，缺测未当阴性。','',
'验收',
f"876648个唯一细胞与全细胞注释逐ID一致；107样本/9来源完整；原8174细胞的DSC1原始值逐细胞重现Phase66；分组总数/阳性数/计数总量一致；独立按原H5AD直接读取{len(checks)}个细胞与提取结果一致。",
'分析退出码0，ANALYSIS_COMPLETE与QA/ACCEPTANCE.json通过。','',
'参考身份映射：https://www.rndsystems.com/jp/products/desmocollin-1-antibody-sr1004_nbp3-22640']
report=O/'delivery/DSC1全细胞表达核查_20261005.txt';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
qa=json.loads((O/'QA/ACCEPTANCE.json').read_text());qa['independent_direct_H5AD_checks']=len(checks);qa['independent_direct_H5AD_checks_passed']=True
(O/'QA/ACCEPTANCE.json').write_text(json.dumps(qa,indent=2)+'\n')
manifest=[]
for sub in ['code','tables','resources','QA','logs']:
 for p in sorted((O/sub).rglob('*')):
  if p.is_file():manifest.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for p in [report,O/'ANALYSIS_COMPLETE',O/'analysis.log']:
 manifest.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'delivery/FILE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
z=O/'delivery/DSC1_all_cells_20261005.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for v in manifest:f.write(O/v['path'],v['path'])
 f.write(O/'delivery/FILE_MANIFEST.json','delivery/FILE_MANIFEST.json')
with zipfile.ZipFile(z) as f:assert f.testzip() is None
receipt={'report':str(report),'zip':str(z),'zip_sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'bytes':z.stat().st_size,'files':len(manifest)+1,'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'delivery/DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
print(json.dumps(receipt,ensure_ascii=False))
print(sources[['dataset','assayed_cells','positive_cells','positive_pct','raw_total','max_raw']].to_string(index=False))
print(lineages[['review_lineage_code','assayed_cells','positive_cells','positive_pct','raw_total','max_raw']].to_string(index=False))
print(dist.to_string(index=False))
