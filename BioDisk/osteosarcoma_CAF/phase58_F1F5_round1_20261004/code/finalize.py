from common import *
import zipfile,shutil,importlib.metadata,platform
for stage,mark in [('freeze','FREEZE_COMPLETE'),('statistics','STATISTICS_COMPLETE'),('gsea_pipeline','GSEA_COMPLETE'),('correction','CORRECTION_COMPLETE'),('characterize','CHARACTERIZATION_COMPLETE'),('QA','INDEPENDENT_QA_COMPLETE')]:
 assert (O/(stage+'.exit_code')).read_text().strip()=='0',stage
 assert (O/mark).is_file(),mark
assert (O/'SYNTHESIS_COMPLETE').is_file()
q=json.loads((O/'delivery/QA_summary.json').read_text());assert q['P_values_checked']==25972 and q['GSEA_files']==88 and q['qa_exit_code']==0
frozen=read(O/'tables/frozen_expanded_cells.tsv.gz');assert len(frozen)==8174 and frozen.cell_id.is_unique
assert frozen.class_id.value_counts().to_dict()=={'F4':5304,'F5':1699,'F1':695,'F3':251,'F2':225}
g=read(O/'tables/all_GSEA_results.tsv.gz');assert len(g)==206889 and not g.duplicated(['dataset','contrast','library','Term']).any()
layout=json.loads((O/'delivery/PDF_layout_QA.json').read_text());assert layout['pages']==17 and not layout['missing_glyphs']
qa={'final_report_pages':17,'all_pages_visually_reviewed':True,'final_changed_pages_rechecked':[1,3,4,5,6,7,16],'Poppler_rendered_all_pages':True,'fonts_and_text_bounds_passed':True,'statistical_QA':q,'frozen_cells':8174,'full_union_features':63398,'GSEA_term_rows':len(g),'GSEA_files':88,'GSEA_new_files':44,'GSEA_verified_reused_files':44,'round1_complete':True,'round2_started':False,'CNV_rerun':False,'at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
dump(O/'delivery/FINAL_ACCEPTANCE.json',qa)
dump(O/'delivery/WORKFLOW_REPAIRS.json',{'QA_controller_first_attempt':'pidfd_open missing in batch_bench Python; no scientific calculations failed. Controller wait switched to system Python; completed dependencies retained; retry exit 0. Original error log retained.','PDF_QA':'Replaced unavailable superscript/subscript glyphs with ASCII. PDF verifier excluded project statistics.py from Python import lookup to avoid standard-library name collision. Final font coverage and bounds passed.'})
versions={'python':platform.python_version()}
for name in ['numpy','pandas','scipy','anndata','statsmodels','gseapy','matplotlib']:
 try:versions[name]=importlib.metadata.version(name)
 except importlib.metadata.PackageNotFoundError:versions[name]='not installed'
dump(O/'resources/analysis_software_versions.json',versions)
obj=[]
for f in sorted((O/'objects').rglob('*')):
 if f.is_file():obj.append({'server_path':str(f),'bytes':f.stat().st_size,'SHA256':sha(f),'included_in_zip':False})
pd.DataFrame(obj).to_csv(O/'delivery/SERVER_OBJECTS_SHA256.tsv',sep='\t',index=False)
readme='''CAF F1–F5 第一轮交付（2026-10-04）

对象：phase30 原始拓展层 8,174 个细胞，F1=695、F2=225、F3=251、F4=5304、F5=1699。
本轮：五类全转录组特征、33 个既有/原始内部状态、校正前后敏感性、表面候选及 8 类竞争、既有临床关联整理、四库 GSEA 与 ORA。
没有重新运行 CNV；第二轮蛋壳状 CAF 寻找尚未启动。

先看 delivery/CAF_F1F5_round1_report.pdf（17 页）或同目录纯文本。
完整基因统计：tables/genes/。原始和校正文件、逐项对比均保留。
程序统计：tables/all_program_statistics.tsv.gz 与 corrected_program_statistics.tsv.gz。
GSEA：tables/all_GSEA_results.tsv.gz（206,889 行），gsea/ 含 88 份独立结果；resources/*__rank.tsv 为完整排序。
ORA：tables/ora/；覆盖：tables/pathway_coverage.tsv.gz。
内部状态：tables/internal_state_composition.tsv、internal_customer_programs.tsv.gz；全量基因/程序以 kind=internal 标记。
表面候选及竞争：tables/surface_*。
临床历史与严格修订：tables/reused/；P_t/FDR_t 是 HC3 调整模型的修订统计列。
独立 QA：independent_QA/；验收：delivery/FINAL_ACCEPTANCE.json。
脚本和固定定义：code/、resources/；生成报告的脚本为 Windows 本机版本，分析脚本为 sgs10 版本。

统计解释：患者为独立样本量；N<2/全零不作检验；N<10 属探索性。以原始 P<0.05 提名，完整同类家族 FDR 同时保留。GSEA 模型 P 与患者级 P 分开。有限置换产生的 0 保留且在报告中用 0* 提醒，不能解释为真实概率为零。效应为相对表达，不代表绝对分泌量。

复用：相同旧比较核实效应和 feature 顺序后复用原始 P；本轮 FDR 按统一五类家族重算。44 份 GSEA 新计算、44 份经完整排序一致验证后复用。
可追溯性：DELIVERY_FILES_SHA256.tsv 列出包内主要文件；SERVER_OBJECTS_SHA256.tsv 给出未装入 ZIP 的大型 H5AD/NPZ 的服务器路径和校验值。原始输入路径与 SHA256 见 resources/INPUT_MANIFEST* 及代码。
服务器目录：/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase58_F1F5_round1_20261004
'''
(O/'README_交付说明.txt').write_text(readme,encoding='utf-8')
(O/'ROUND1_COMPLETE').write_text('SUCCESS\n')
files=[]
for folder in ['code','resources','tables','figures','independent_QA','logs','delivery']:
 for f in (O/folder).rglob('*'):
  if not f.is_file() or '__pycache__' in f.parts or f.suffix=='.zip' or any(s in f.parts for s in ['pdf_render','pdf_render_final']):continue
  if f.name=='DELIVERY_FILES_SHA256.tsv':continue
  files.append(f)
files.extend(f for f in O.iterdir() if f.is_file() and (f.name.endswith('_COMPLETE') or f.suffix=='.exit_code' or f.name.startswith('README_')))
files=sorted(set(files));manifest=[{'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'SHA256':sha(f)} for f in files]
pd.DataFrame(manifest).to_csv(O/'delivery/DELIVERY_FILES_SHA256.tsv',sep='\t',index=False);files.append(O/'delivery/DELIVERY_FILES_SHA256.tsv')
out=O/'delivery/CAF_F1F5_round1_complete_results_20261004.zip'
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 for f in files:z.write(f,arcname=str(f.relative_to(O)))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
res={'file':str(out),'bytes':out.stat().st_size,'SHA256':sha(out),'files':len(files),'integrity':'ZIP_CRC_PASS'}
dump(O/'delivery/ZIP_ACCEPTANCE.json',res)
# Preserve the former current pointer, then record the completed original-object first round.
current=P/'CAF_NEXT_STAGE_CURRENT.json';stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
if current.exists():shutil.copy2(current,O/'delivery'/('previous_current_before_round1_'+stamp+'.json'))
state={'project':'osteosarcoma_CAF','current_phase':'phase58_F1F5_round1_20261004','status':'round1_complete','scope':'original_expanded_8174_F1_F5_full_transcriptome','report':str(O/'delivery/CAF_F1F5_round1_report.pdf'),'results_zip':str(out),'acceptance':str(O/'delivery/FINAL_ACCEPTANCE.json'),'round2_eggshell_CAF':'not_started','CNV_rerun':False,'prior_phase55_57_role':'ancillary_annotation_and_reclassification_not_replacement_for_original_round1','updated_at':qa['at']}
dump(current,state);dump(O/'delivery/CURRENT_PROGRESS.json',state);print(json.dumps(res),flush=True)