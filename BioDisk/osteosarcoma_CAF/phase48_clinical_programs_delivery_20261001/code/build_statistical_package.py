from pathlib import Path
import json,hashlib,zipfile,datetime,os
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');Q=P/'phase48_clinical_programs_delivery_20261001';D=Q/'delivery/client';D.mkdir(exist_ok=True)
phases=['phase41_full_transcriptome_setup_20261001','phase42_F3F4_identity_states_20261001','phase43_full_transcriptome_discovery_20261001','phase44_patient_mechanism_20261001','phase45_LR_expression_20261001','phase46_spatial_functional_20261001','phase47_external_candidate_integration_20261001','phase48_clinical_programs_delivery_20261001']
files=[]
for phase in phases:
 root=P/phase
 for f in root.rglob('*'):
  if not f.is_file() or '__pycache__' in f.parts or f.suffix=='.part':continue
  rel=f.relative_to(root);parts=rel.parts
  if any(x in parts for x in ['downloads','objects','client','pdf_render']):continue
  if f.suffix=='.zip' or f.name.startswith('report_input') or f.name=='report_update.zip':continue
  if parts[0] in ['tables','resources','code','figures','logs','delivery'] or len(parts)==1:files.append(f)
# Refer to full raw/native matrices without duplicating large matrices in the archive.
import pandas as pd
inputs=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/source_sample_gene_coverage.tsv',sep='\t');inputs[['sample','dataset','input_path','input_bytes','source_sha256_recorded','counts_sha256_recorded','gene_universe_sha256']].to_csv(D/'RAW_INPUT_REFERENCES.tsv',sep='\t',index=False)
objects=[]
for phase in phases:
 for f in (P/phase/'objects').glob('*'):
  if f.is_file():objects.append({'phase':phase,'path':str(f),'bytes':f.stat().st_size,'role':'raw/aggregated analysis object kept on configured server; path reconstructible from scripts'})
(D/'SERVER_OBJECT_REFERENCES.json').write_text(json.dumps(objects,indent=2))
readme='''CAF全转录组与空间基质研究统计包 2026-10-02

阅读顺序
1. 报告PDF与FINAL_ACCEPTANCE.json说明最终结论及模块验收。
2. phase47的prioritized_candidates.tsv、main_functional_program_patient_statistics.tsv与source_projection_statistics.tsv提供候选和来源内效应。
3. 空间主要统计为phase46/tables/spatial_patient_statistics.tsv。当前组织图及坐标验收为phase48/tables/histology_image_coordinate_audit.tsv，不能将旧的AnnData内嵌图像检查误读为没有图像。
4. 当前比例调整模型使用phase48/tables/D_fraction_HC3_current_statistics.tsv。P_t/FDR_t为t参考，P/FDR为数学可评估模型的原normal参考，original_normal_P及original_normal_FDR只用于历史审计。
5. 当前连续程序统计为phase48/tables/clinical_program_statistics.tsv。初版误将FIELD_ABSENT作为记录协变量及单位杠杆模型已修正；initial文件与失败日志保留，不能作为最终推断结果。
6. 表面新候选为phase48/tables/new_surface_discovery_statistics.tsv、surface_competing_lineage_statistics.tsv及RNA共检出表。

统计单位及结果类型
患者是生物学N；细胞和spot是测量单位。有效目标原P<0.05用于名义候选筛选，完整FDR与检验家族始终保留。N<10按探索性证据呈现，原P<0.05但FDR未通过不称校正后显著。缺重复、真实全零、未测成员与不可评估模型保持NA。
GSEA的P为有限基因集置换模型P，软件原生P=0不表示数学精确概率为0。TF为有符号靶程序响应，配体受体为RNA兼容性；它们不是直接蛋白活性、实际通讯或生物学处理效果。
冻结支持层1224与扩展层8174保留，同编号F3不是同一个细胞集合。F4身份冲突未被删除或转换成CAF确认。跨来源未知患者重叠未合并成独立N。

几何与实验边界
组织图像可用并已核查RNA坐标。没有独立可验证的肿瘤边界、癌/基质区域或CD8单细胞掩膜，环绕肿瘤的连续覆盖与内外定位终点按计划不可评估。局部RNA空间分析已经完成。蛋白分选门与真实保护作用需后续实验，本轮不声称完成实验。
OsteoCAR/Budhathoki来源获取与跨来源全局患者对应未闭合；已有来源审计保留具体状态，不把重复来源当作独立复现。

复核与输入
包内scripts、固定资源及全量统计可追溯。RAW_INPUT_REFERENCES.tsv和SERVER_OBJECT_REFERENCES.json引用sgs10项目中的大体积原始计数及分析对象。FILE_MANIFEST.tsv列出包内每个文件的字节数与SHA256。已通过的独立R核查、各当前退出码及完成标记在FINAL_ACCEPTANCE.json注册；旧控制器失败状态以该注册为准。
所有新增下载使用后台直连，禁用环境代理。历史输入、旧报告及失败日志未删除。
'''
(D/'README_阅读说明.txt').write_text(readme,encoding='utf-8');files.extend([D/'README_阅读说明.txt',D/'RAW_INPUT_REFERENCES.tsv',D/'SERVER_OBJECT_REFERENCES.json',D/'CAF客户结果报告.pdf',D/'FINAL_ACCEPTANCE.json',D/'CONDITIONAL_ENDPOINTS.tsv']);files=sorted(set(files));records=[]
def digest(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
# Keep client provenance at archive root; otherwise retain phase-relative directories.
entries={str(f.relative_to(D)) if D in f.parents else str(f.relative_to(P)):f for f in files};assert len(entries)==len(files)
archive=D/'CAF全转录组与空间基质研究_完整统计包_20261002.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
 for name,f in entries.items():
  data=f.read_bytes();z.writestr(name,data);records.append({'path':name,'bytes':len(data),'SHA256':hashlib.sha256(data).hexdigest()})
 manifest=pd.DataFrame(records);payload=manifest.to_csv(sep='\t',index=False);z.writestr('FILE_MANIFEST.tsv',payload);(D/'FILE_MANIFEST.tsv').write_text(payload)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None and len(z.namelist())==len(records)+1
status={'at':datetime.datetime.now().astimezone().isoformat(),'state':'package_integrity_accepted','archive':str(archive),'archive_bytes':archive.stat().st_size,'archive_SHA256':digest(archive),'files':len(records),'CRC_passed':True,'manifest_all_files':True};(D/'PACKAGE_INTEGRITY.json').write_text(json.dumps(status,indent=2,ensure_ascii=False));(Q/'STATISTICAL_PACKAGE_COMPLETE').write_text('SUCCESS');print(json.dumps(status),flush=True)
