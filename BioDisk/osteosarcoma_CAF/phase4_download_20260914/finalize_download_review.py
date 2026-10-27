"""SYNTAX_CHECK: python -m py_compile finalize_download_review.py
RUN_COMMAND: python finalize_download_review.py --root /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914
EXPECTED_SUCCESS_MARKER: NON_PDX_REVIEW_READY
"""
from pathlib import Path
import argparse,csv,json,hashlib,zipfile
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
assert (root/'COMPLETE').read_text().strip()=='NON_PDX_DOWNLOAD_VALIDATION_COMPLETE'
rows=json.loads((root/'matrix_validation.json').read_text())
clinical=list(csv.DictReader((root/'clinical_source_audit.tsv').open(encoding='utf-8-sig'),delimiter='\t'))
lookup={(r['dataset'],r['sample']):r for r in clinical}
for r in rows:
 assert r['validation']=='PASS' and r['direct_patient_verified'] and r['is_xenograft_verified'] is False
 r['clinical_group']=lookup[(r['dataset'],r['sample'])]['clinical_group_verified']
 r['patient_history_note']=lookup[(r['dataset'],r['sample'])]['unresolved']
with (root/'analysis_entry_manifest.tsv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
assert len(rows)==17
table='\n'.join('| '+ ' | '.join(str(r[k]) for k in ['dataset','sample','genes','input_barcodes','clinical_group'])+' |' for r in rows)
downloads=json.loads((root/'download_checks.json').read_text())
text=f'''# 非 PDX 队列下载与输入核查完成

GSE299023的6个样本及GSE152048的11个组织矩阵已实际下载并读取，共17个组织、{sum(r['input_barcodes'] for r in rows):,}个输入barcode。这里的细胞数为下载矩阵列数，不是本地重新QC后的细胞数或CAF数量。

下载文件共{len(downloads)}个，完整文件共{sum(r['bytes'] for r in downloads):,}字节。每个文件均记录SHA256，压缩数据通过读取校验。矩阵与基因/barcode维度一致，counts非负有限整数，barcode在样本内唯一；基因重复单独记录，不静默丢弃。基因ID未检测到ENSMUSG前缀仅为辅助核查，非PDX身份依据取材方法而非人基因注释。

| Dataset | Sample | Genes | Input barcodes | Clinical group |
| --- | --- | ---: | ---: | --- |
{table}

## 来源判定

两队列原始论文明确为患者直接取材组织并直接解离测序，无PDX传代环节。GSE299023样本中心和编号与论文补充表一致；GSE152048的76个GEO记录对应11个组织矩阵。官方缺少独立is_xenograft字段时保留field_absent，文献核实结果单独保存为false。PDX、PDX衍生培养、来源未确认样本不在本轮分析入口中。

参考来源：[GSE299023论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12248827/)，[GSE152048论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/)。论文XML、补充临床表、样本清单和对应证据已保存。来源审核与临床分组的原始值、核实值分开记录于clinical_source_audit.tsv。

## 临床分组与下一步

BC10、BC17为明确肺转移灶，进入B2。其余样本保留Unassigned，不能凭原发灶或治疗前标签归入A。BC11、BC20是复发组织，不按首次骨原发处理。

GSE152048的11例均在补充表记录术前化疗，BC17另有抗PD-L1暴露；GSE299023为治疗前骨原发样本。因此两个队列之间的治疗与技术差异需保留，不能直接解释为转移效应。

GSE299023作者的QC后数量中，OS.C.21.BX.OCT和OS.C.32.BX分别只有126和144个细胞。这是作者结果，不能替代本地QC；后续先评估CAF数量和信息量。不能预设17个组织都能提供可用CAF。

analysis_entry_manifest.tsv仅表示通过非PDX来源和表达输入检查、可进入质控/身份分析。CAF注释、双细胞评估、两种CNV及严格/扩展CAF集合尚未在这17个样本上运行。

## 文件与运行

原始下载位于downloads/，解包表达文件位于matrices/，大型矩阵独立保存。download_checks.json提供文件校验值；matrix_validation.json和analysis_entry_manifest.tsv提供逐矩阵检查结果。运行及传输命令见运行与监控.md。

早期FTP直连速度过低，保留原日志和部分下载，最终采用GEO官方备用下载入口。最终状态以COMPLETE、download_geo.log末尾成功标记以及逐样本PASS共同确认。
'''
(root/'确认报告.md').write_text(text,encoding='utf-8')
include=['确认报告.md','运行与监控.md','analysis_entry_manifest.tsv','clinical_source_audit.tsv','download_manifest.json','download_checks.json','matrix_validation.json','source_and_matrix_validation.tsv','environment.json','resource_usage_http1.txt','resource_usage_geo.txt','resource_usage.txt','download_http1.log','download_geo.log','download_validate.log','COMPLETE','started_at.txt','geo_started_at.txt','finished_at.txt','PMC7730477.xml','PMC12248827.xml','clinical_page_25.png','Supplementary Table S1 - Sample characteristics.xlsx']
include += [p.name for p in root.glob('*.py')]+[p.name for p in root.glob('*.sh')]
with zipfile.ZipFile(root/'non_PDX_download_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for name in dict.fromkeys(include):
  if (root/name).is_file():z.write(root/name,name)
(root/'PAUSED_AFTER_DOWNLOAD_VALIDATION').write_text('User requested pause after download and validation. No new cohort QC, clustering or CNV started.\n')
print('NON_PDX_REVIEW_READY',len(rows),sum(r['input_barcodes'] for r in rows))
