from pathlib import Path
import json,hashlib,datetime,pandas as pd,re,shutil
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005';B=P/'phase62_F1F5_internal_20261005'
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(8388608),b''):h.update(z)
 return h.hexdigest()
for d in ['resources','tables','objects','figures','QA','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
libs=json.loads((P/'phase58_F1F5_round1_20261004/resources/libraries_canonical.json').read_text())
specs={
'collagen_biosynthesis':('reactome','REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES','main_matrix_mechanisms'),
'collagen_assembly':('reactome','REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES','main_matrix_mechanisms'),
'collagen_crosslinking':('reactome','REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS','main_matrix_mechanisms'),
'collagen_degradation':('reactome','REACTOME_COLLAGEN_DEGRADATION','main_matrix_mechanisms'),
'ECM_degradation':('reactome','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX','main_matrix_mechanisms'),
'ECM_organization':('reactome','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','main_matrix_mechanisms'),
'collagen_components':('matrisome','NABA_COLLAGENS','main_matrix_mechanisms'),
'contraction':('reactome','REACTOME_SMOOTH_MUSCLE_CONTRACTION','secondary_functional_axes'),
'mineralization':('go_bp','GOBP_BONE_MINERALIZATION','secondary_functional_axes'),
'inflammation':('hallmark','HALLMARK_INFLAMMATORY_RESPONSE','secondary_functional_axes')}
sets={k:sorted(set(libs[lib][term])) for k,(lib,term,family) in specs.items()};roles={k:v[2] for k,v in specs.items()};sets['LOX_family']=['LOX','LOXL1','LOXL2','LOXL3','LOXL4'];roles['LOX_family']='secondary_target_panel';dump(O/'resources/PROGRAMS.json',{'sets':sets,'roles':roles,'definitions':specs,'panel_note':'LOX family is a targeted RNA panel, not enzyme activity or a validated classifier'})
protocol={'date':datetime.datetime.now(datetime.timezone.utc).isoformat(),'objective':'Find cell states and targetable processes relevant to collagen production, assembly, crosslinking and matrix deposition around tumors and treatment resistance; classification is a tool. Eggshell geometry / immune neighborhood is a side branch.','cells':'Freeze all original8174 and33 subgroup IDs and all confidence fields; no new global clustering or CNV.','primary_sc':'Every subgroup vs other original8174 cells, matched source/donor/site; raw and corrected; full-library CP10k denominator. Entire fixed pathway sets, not selected top markers. Other-F background serves mechanism prioritization only, not an identity gate.','sc_sensitivity':'Corrected same contrast minimum5 cells each arm, and corrected Fibro-restricted target/reference minimum1. Main minimum1 is mathematical evaluability. No biological cell exclusion.','patient_effect':'Mean gene log2(CP10k+.1) target-minus-reference, source-specific measurement coverage. Site-equal then independent source-local donor-equal; record gene-specific N and program-specific N. Existing80% pathway coverage and min(3,size) required in a pair; unassayed never zero.','statistics':'Exact two-sided patient signflip N<=15 else19999 Monte Carlo+observed seed20261001 using retained phase58 helper; N<2/allzero no P; N<10 exploratory. Nominal raw target P<.05 hard criterion when tests valid. FDR retained per source/matrix/scope/threshold and main vs secondary program families. Genes form separate family across33 subgroup contrasts.','gene_scope':'Focused mechanism gene follow-up includes union of every fixed pathway gene; broad full-transcriptome discovery and within-F GSEA from phase62 retained with original families, explicitly distinct from new targeted tests.','ranking':'Mechanism-specific evidence ledger with positive nominal patient support, descriptive positive direction, sensitivity contradictions and negative results; no arbitrary composite score or mandatory fixed cluster count.','spatial':'All26 GSE293065 slices,19 patient records,53309 spots; frozen raw counts, author signed cancer/fibroblast components and prior adjacency. Cell program scores mapped as programs, not automatic33 subgroup assignments.','spatial_endpoints':'Primary mainline: program vs local cancer signed component, adjacent cancer component and fibroblast component; raw and composition/UMI/XY-adjusted where registered. Side: local/adjacent CD8 RNA associations. Per-slice effects average event then site then patient. No spot-as-N P. Patient signflip and BH per endpoint x signature_version x main/secondary program family.','spatial_sensitivity':'Remove author component contributor genes from program definitions; record resulting membership and reduced-program status. Minimum3 remaining genes, except original panels of size<3. This diagnostic sensitivity is not the original full program. Primary-site-only patient estimates retained.','geometry':'No accepted independent tumor masks; do not manufacture tumor boundary from tested matrix RNA. Report tumor-neighborhood association; circumferential coverage and physical shell geometry remain conditional. Existing histology images and overlaps retained.','clinical':'Audit integrated107sample discovery and26slice metadata for actual Huvos, necrosis, drug response or dose-response endpoints. Treatment history, recurrence, metastasis, survival status are not resistance labels. Absent fields mean unavailable in current integrated data, not proven absent from all original supplements. No fabricated resistance tests.','delivery':'Full plan, executed mechanism and spatial results, evidence ledger, drug-response availability, eggshell side evidence, conditional next steps; maintain original labels/confidence.'}
dump(O/'resources/PROTOCOL_FROZEN.json',protocol);dump(O/'resources/PREVIOUS_CURRENT.json',json.loads((P/'CAF_NEXT_STAGE_CURRENT.json').read_text()))
paths=[B/'tables/F1F5_internal_cell_labels.tsv.gz',B/'tables/F1F5_subgroup_summary.tsv',P/'phase58_F1F5_round1_20261004/resources/libraries_canonical.json',P/'phase33_client_report_20260926/tables/final_clinical_sample_manifest.tsv',P/'phase36_external_validation_20260927/resources/sample_patient_manifest_corrected.tsv',P/'phase48_clinical_programs_delivery_20261001/delivery/HISTOLOGY_IMAGE_AUDIT.json']
dump(O/'resources/INPUT_MANIFEST.json',[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in paths])
available=[]
for path in paths[3:5]:
 f=pd.read_csv(path,sep='\t',low_memory=False);response=[c for c in f.columns if re.search('huvos|necros|drug_response|chemotherapy_response|resistance|ic50|auc_drug',c,re.I)]
 available.append({'file':str(path),'records':len(f),'datasets':f.dataset.nunique(),'patients_source_local':f.donor_record.nunique() if 'donor_record' in f else f.patient_id.nunique(),'response_fields':' ; '.join(response),'direct_response_available':bool(response),'status':'current_integrated_metadata_has_no_direct_response_endpoint' if not response else 'requires_value_and_endpoint_review'})
 shutil.copy2(path,O/'resources'/path.name)
pd.DataFrame(available).to_csv(O/'tables/clinical_endpoint_availability.tsv',sep='\t',index=False)
for phase,file in [('phase46_spatial_functional_20261001','spatial_patient_statistics.tsv'),('phase50_F3_spatial_decomposition_20261002','statistics.tsv')]:shutil.copy2(P/phase/'tables'/file,O/'tables'/('retained_'+phase.split('_')[0]+'_'+file))
shutil.copy2(paths[-1],O/'resources/HISTOLOGY_IMAGE_AUDIT.json')
g=pd.read_csv(B/'tables/all_GSEA_results.tsv.gz',sep='\t');g=g[g.Term.isin([v[1] for v in specs.values()])];g.to_csv(O/'tables/reused_full_transcriptome_withinF_GSEA.tsv.gz',sep='\t',index=False)
plan='''骨肉瘤基质相关治疗抵抗研究 下一阶段执行计划
日期 2026-10-05

最终目标
定位参与肿瘤周围胶原生成、纤维装配与交联、基质沉积的候选细胞及可干预环节，检验其与药物作用或治疗应答的关系。蛋壳状CAF及局部免疫分布为侧线，不以必须找到环形结构或新的CAF类型作为完成条件。

执行顺序与交付
第一步 全部亚群的机制定位 本轮执行
以原8174细胞和33内部亚群为索引，保留全部置信度。分别检验胶原合成、装配、交联、胶原降解、ECM降解、ECM组织及胶原组分；收缩、矿化、炎症和LOX家族另列。每个亚群与同患者同部位其他候选细胞比较，复核原始/去污染、细胞数敏感性及Fibro限定结果。结合已有全转录组结果形成按机制分列的候选表，不将降解重塑高直接写成净沉积高，不用综合分数掩盖矛盾。

第二步 空间定位 本轮执行可评估部分
使用26张空间切片及19个患者记录，将预先冻结的完整机制程序映射到组织坐标，分析其与肿瘤和成纤维样表达成分、相邻肿瘤成分的关系。以患者为单位汇总；记录组成调整和去除成分模型贡献基因后的敏感性，以及原发部位子集。程序定位与具体亚群定位分别表述，不把共有胶原RNA直接当作某个CAF亚群定位。

第三步 治疗应答关联 本轮完成数据条件核查
逐项核对已整合的107份单细胞样本台账和26张空间样本台账，寻找病理坏死率/Huvos、药物反应和药敏资料。现有治疗史、复发/转移及生存状态不得改造成耐药标签。若无直接终点，交付缺口表并保留来源与样本联结入口；下一步定向核查原始补充资料或接入明确授权的临床表，具备有效患者对比后再冻结并执行治疗应答模型。当前不补造耐药显著性。

第四步 蛋壳侧线 本轮保留并补充可评估证据
复用已验收的三基因基质/CD8空间结果并完整保留其模型依赖。新增完整机制程序与局部/邻近CD8 RNA的患者级关联。26张组织图像存在，但独立肿瘤边界尚未验收，因此本轮不计算或宣称已确认蛋壳环绕程度。后续获得独立组织学边界后，预先登记边缘富集、周向覆盖、连续性和内外分布终点。

第五步 候选收敛与干预验证设计 本轮提出数据驱动优先环节
候选细胞、机制分子和检测标志分开列出。根据患者效应、空间程序证据、来源和质量敏感性决定优先级。实验主终点应分别包括实际胶原沉积/纤维组织、药物进入或作用以及肿瘤细胞药物反应，设置匹配的药物和基质干预对照；实际实验不计为本轮计算已完成内容。

统计与验收
独立患者为N；有效检验按原始P<0.05筛选名义候选，分别保留同类家族FDR。N<10定位探索性；缺重复/全零不造P。空间spot不作独立样本；现有IC为带符号表达成分，不称细胞比例。发现资料与空间来源的全局患者独立性未完全解决，来源不合并作新增独立验证。
本轮完成条件为计划与程序定义冻结、可评估机制和空间模块实际运行、临床及几何缺口明确、候选证据与阴性矛盾完整、标签和统计通过必要验收。仅将本轮计算阶段标记完成，不将真实耐药机制或实验干预标记为已证实。
'''
(O/'delivery/CAF基质耐药目标导向执行计划_20261005.txt').write_text(plan,encoding='utf-8');(O/'FREEZE_COMPLETE').write_text('SUCCESS\n');print('FREEZE_COMPLETE',len(sets),sum(len(v) for v in sets.values()),'membership entries',flush=True)
