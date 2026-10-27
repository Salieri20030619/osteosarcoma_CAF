from pathlib import Path
import json,datetime,hashlib,shutil,zipfile
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase93_expanded_admission_20261016';S=P/'phase92_copykat_validation_20261015'
def dump(f,d):f.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
assert (O/'logs/build.exit_code').read_text(encoding='utf-8').strip()=='0'
build=json.loads((O/'BUILD_COMPLETE.json').read_text(encoding='utf-8'));assert build['registry_cells']==8980
assert all(x['pass'] for x in json.loads((O/'resources/FINAL_INPUT_INTEGRITY.json').read_text(encoding='utf-8')))
for name in ['admitted_expanded_cells.tsv.gz','expanded_registry_8980.tsv.gz']:
 f=O/'tables'/name;c=pd.read_csv(f,sep='\t',encoding='utf-8')
 if 'formal_admission' in c:c=c.drop(columns='formal_admission')
 if 'not_ready_for_formal_admission' in c:c=c.rename(columns={'not_ready_for_formal_admission':'historical_not_ready_for_formal_admission'})
 c.to_csv(f,sep='\t',encoding='utf-8',index=False)
clusters=pd.read_csv(O/'tables/cluster_review.tsv',sep='\t',encoding='utf-8')
clusters['final_scope']='retain_expanded_candidate_with_original_identity_risks'
m=clusters['sample'].eq('NCC__BS3')&clusters.cluster.eq(4)
clusters.loc[m,'final_scope']='retain_expanded_osteochondral_matrix_overlap'
clusters.loc[m,'manual_reason']='Supplemental raw-count review: COL2A1/COL10A1/COL11A1/FN1/IBSP/SPP1 with COL1A1/COL1A2/LUM/DCN; osteochondral-matrix overlap, not confirmed pure fibroblast. Immune/endothelial/myogenic markers not dominant. Descriptive expression, not a new DE test.'
clusters.to_csv(O/'tables/cluster_review.tsv',sep='\t',encoding='utf-8',index=False)
cnv=pd.read_csv(O/'tables/chromosome_profile_cell_diagnostics.tsv.gz',sep='\t',encoding='utf-8');risk=[]
for sm,g in cnv.groupby('sample'):
 refmax=g[g.group.eq('heldout_normal')].max_abs_chromosome_mean.max()
 q=g[g.group.eq('candidate')].copy();q['above_observed_heldout_max']=q.max_abs_chromosome_mean.gt(refmax)
 risk.append({'sample':sm,'candidate_cells':len(q),'above_observed_heldout_max':int(q.above_observed_heldout_max.sum()),'heldout_max':float(refmax),'role':'Descriptive risk only; not an independent malignancy test or new admission threshold.'})
dump(O/'resources/CNV_PROFILE_REVIEW.json',{'utc':now,'review':'Candidate group chromosome medians are lower in maximum absolute amplitude than corresponding sampled native aneuploid context in all seven libraries. Cell-level tails and lineage overlap persist; native diploid is not proof of pure CAF or genomic normality. No new exclusion cutoff applied. Aneuploid context is first100 sorted non-heldout IDs, not a random patient sample.','samples':risk,'HMM':'Any-segment binary rejected intended-normal controls and remains descriptive only; native outputs are preserved under Phase90/91.'})
disposition=[]
for sm in ['NCC__BS2','NCC__BS9','GSE234187__GSM7453693']:
 v=json.loads((S/'resources'/(sm+'_VALIDITY_RECHECK.json')).read_text(encoding='utf-8'))
 disposition.append({'sample':sm,'decision':'not_admitted_in_this_version','reason':'insufficient_evaluable_Mast_reference_9_below_10' if sm.startswith('GSE') else 'reference_classification_failed','checks':v['CopyKAT_reference_checks'],'raw_evidence_preserved':True,'threshold_relaxed':False})
dump(O/'resources/EXCLUDED_LIBRARY_DISPOSITIONS.json',disposition)
copyfiles={
 'native_computation_reconciled.json':S/'COMPUTATION_RECONCILED.json',
 'patient_counts_before_admission.tsv':S/'patient_count_verification/patient_counts.tsv',
 'NCC_provenance.json':P/'phase90_patient_validation_20261015/resources/PATIENT_PROVENANCE_REVIEW.json',
 'NCC_fingerprint_complete.json':P/'phase90_patient_validation_20261015/resources/FINGERPRINT_COMPLETE.json',
 'lung_provenance.json':P/'phase91_lung_patient_validation_20261015/resources/SOURCE_PROVENANCE.json',
 'HMM_role_amendment.json':S/'resources/HMM_ROLE_AMENDMENT.json'}
for n,f in copyfiles.items():shutil.copy2(f,O/'resources'/n)
for f in (S/'resources').glob('*_VALIDITY_RECHECK.json'):shutil.copy2(f,O/'resources'/f.name)
report='''CAF 拓展层扩充与患者数量终核（内部记录，2026-10-16）

结论：本轮拓展层扩充完成。新增 NCC 来源内6位作者标识患者、7份文库、806个拓展候选细胞；与历史8174个细胞组成8980个细胞的版本化登记资源。新增细胞及患者明确按拓展探索层纳入，严格层不作为统一入库门槛。

患者清单（新增拓展细胞 / 去污染支持子集 / 旧严格支持标记子集）
Patient1：41 / 9 / 0
Patient3：53 / 19 / 0
Patient4：277 / 162 / 0；BS4和BS14合计，只计1位患者。
Patient5：56 / 46 / 0
Patient7：351 / 325 / 176
Patient12：28 / 11 / 0
合计：806 / 572 / 176。三个数字是包含关系及证据分层，不是相加的三批细胞。

本轮检查与结果
1. 10份文库、92755个完整provider-called细胞的原生CopyKAT均完成，逐文库退出码为0且完成标记存在。7份文库通过预设正常参考诊断；BS3只有可用免疫参考，保留其参考谱系覆盖限制。主调度器曾因pidfd接口缺失异常；未启动的第二路已单独恢复，最终按逐样本记录核对，无缺失和重复运行计数。
2. 新增806个细胞均来自已人工复核的拓展保留名单，具有完整输入下Fibroblast竞争评分主导、原生diploid分类，无旧/新软件双细胞阳性；明确竞争谱系群不纳入。去污染支持和旧严格层作为附加列保留，不因未进入严格层删除拓展候选。
3. 复核全部新增细胞涉及的群背景。BS3 cluster4补充原始计数表达检查，显示成骨/软骨与基质重叠；保留拓展身份标签。其余身份未定、成骨、收缩/周细胞重叠、CTA及混合谱系风险均保留，不写成确认纯CAF。
4. 对7份纳入文库检查原生CopyKAT染色体均值概览及逐细胞幅度。候选群概览幅度低于各自抽取的原生aneuploid背景，个别细胞尾部信号仍在；不据此宣称基因组正常或另设阈值筛数。HMM任意片段即排除的二分类因正常参考失败退出准入判定，原始区段与调用保留。
5. 806个新细胞与8174个历史细胞ID无重叠；Patient4两文库去重完成。7份新增原始计数对象逐项回读，与对应来源计数完全相同。4个登记基线文件及9份NCC完整输入共13项SHA256检查全部通过。
6. 计数对象按来源保留全部实测基因。COUNT_COLLECTION.json将历史对象、新增7个对象及8980细胞登记表关联起来；未把不同来源未测基因伪填为测得的零。历史33状态保留，新细胞状态待下游重新整合；此次未自动增加或验证状态标签。

本版本不纳入的3份文库
BS2：内皮参考仅8/32 diploid，免疫66/74，未通过参考诊断。
BS9：内皮参考仅2/11 diploid；免疫51/51。失败原因是参考分类，不能据此确认其中候选正常。
肺样本GSE234187：可评估正常参考54/54均diploid；Mast只有9个可评估，低于预设10个，属于覆盖不足而非观察到全部参考异常。保留原始及去污染证据，不将失败标签改成通过，不在本轮扩充人数中计入。

适用范围
本轮新增的是NCC来源内6位患者的拓展分析资源，不声称新增6位已完成跨研究遗传去重的全局独立患者。既有963次NCC与基线矩阵指纹比较未见完全相同计数向量；该检查不能排除同患者不同文库。NCC各样本取材部位仍为Unknown，因此新增骨/肺临床比较有效患者数为0。去污染证据沿用既有预处理保留细胞分支，作为敏感性层，不冒充完整输入下重跑去污染。完整输入用于本轮CNV及双细胞验证；候选身份仍限定于前期人工复核名单，并未宣称对全部新增低质量细胞完成新的CAF发现。
当前数据不支持把8980细胞资源描述为8980个独立确认CAF。该边界与历史拓展层口径一致，不阻止其作为主探索队列使用。本轮没有新增患者层面推断检验，也没有更新旧骨肺统计结果。

完成状态
拓展层名单、患者去重、计数对象、证据分层及完整性验收已完成。未纳入文库已有明确处理决定，不作为本轮等待项。重新整合/分群和重跑下游统计属于后续分析，不是此次扩充入库的未完成验收。
'''
(O/'FINAL_REVIEW_ZH.txt').write_text(report,encoding='utf-8')
dictionary='''Main entry flag: admitted_to_expanded_resource. This means admission to a source-stratified expanded exploratory resource, not independently confirmed CAF identity or eligibility for a clinical Bone/Lung comparison.
historical_* admission flags preserve earlier review status and are not current admission decisions.
New cells have source_review_cluster and new_state_assignment=pending_reintegration. Historical cluster/subgroup/state labels remain historical and are not transferred to new cells.
Inherited Phase88/audit scores, raw detections, corrected lineage flags and decontX values are earlier retained-input evidence. native_copykat_call and full_scDblFinder_class are full-input Phase92 evidence.
The raw count collection retains separate feature universes. Use each object's var_names to determine whether a gene was assayed; missing features are not measured zero.
'''
(O/'DATA_DICTIONARY.txt').write_text(dictionary,encoding='utf-8')
acceptance={'utc':now,'status':'PASS_EXPANDED_RESOURCE_ADMISSION','scope':'NCC source-stratified expanded exploratory resource','baseline_cells':8174,'new_expanded_cells':806,'total_registry_cells':8980,'new_source_resolved_patients':6,'source':'NCC_Zheng2024','added_libraries':7,'ambient_supported_subset':572,'legacy_strict_supported_subset':176,'strict_required_for_admission':False,'new_Bone_Lung_patients':0,'pooled_global_independent_N_established':False,'baseline_unchanged':True,'raw_counts_roundtrip_pass':True,'integrity_checks_passed':13,'visual_QA':'Both count and chromosome figures inspected; labels and evidence layers legible.','remaining_required_work_for_this_admission':[],'next_analysis':'Optional subsequent integration and downstream analysis; new state labels not yet assigned.'}
dump(O/'ACCEPTANCE.json',acceptance)
backup=O/'resources/current_before_publication';backup.mkdir(exist_ok=True)
for name in ['CAF_PATIENT_VALIDATION_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_NEXT_STAGE_CURRENT.json']:
 f=P/name
 if f.exists():
  if not (backup/name).exists():shutil.copy2(f,backup/name)
  d=json.loads(f.read_text(encoding='utf-8'))
 else:d={}
 d.update(phase=O.name,current_phase=O.name,directory=str(O),status='expanded_resource_admission_complete',updated_utc=now,new_source_resolved_patients=6,new_expanded_cells=806,total_expanded_cells=8980,new_Bone_Lung_patients=0,scientific_admission_complete=True,admission_scope=acceptance['scope'],acceptance=str(O/'ACCEPTANCE.json'),report=str(O/'FINAL_REVIEW_ZH.txt'),registry=str(O/'tables/expanded_registry_8980.tsv.gz'),count_collection=str(O/'objects/COUNT_COLLECTION.json'),pooled_global_independent_N_established=False)
 d['new_formal_patients']=None;d['new_formal_patients_note']='Unqualified global patient increment is not established; source-resolved expanded admission is6. No strict/pure-CAF gate.'
 d['next_action']='Use versioned expanded8980 registry and source-aware count collection for subsequent integration; do not reuse old33-state labels on new cells.'
 if 'active_task' in d:d['active_task'].update(status='complete',phase=O.name,directory=str(O),scientific_admission_pending=False,new_source_resolved_patients=6,new_formal_patients=None)
 for key in ['current_state','task_status']: 
  if key in d:d[key]='expanded_resource_admission_complete'
 if 'analysis_now_completed' in d:d['analysis_now_completed']=True
 if 'acceptance_passed' in d:d['acceptance_passed']=True
 if 'execution_status' in d:d['execution_status']['exit_status']='completed_verified'
 if 'process_observation' in d:d['process_observation']='All native jobs complete; source-level expanded admission verified.'
 tmp=f.with_suffix('.phase93.tmp');dump(tmp,d);tmp.replace(f)
dump(P/'CAF_EXPANDED_COHORT_CURRENT.json',dict(acceptance,phase=O.name,registry=str(O/'tables/expanded_registry_8980.tsv.gz'),count_collection=str(O/'objects/COUNT_COLLECTION.json'),report=str(O/'FINAL_REVIEW_ZH.txt')))
# Package new deliverables; old source matrices remain referenced and unmodified.
files=[f for f in O.rglob('*') if f.is_file() and 'current_before_publication' not in f.parts and f.suffix!='.zip']
manifest=[]
for f in files:manifest.append({'path':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
dump(O/'DELIVERY_MANIFEST.json',manifest)
z=O/'CAF_expanded8980_admission_20261016.zip'
with zipfile.ZipFile(z,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
 for f in files+[O/'DELIVERY_MANIFEST.json']:archive.write(f,f.relative_to(O))
with zipfile.ZipFile(z) as archive:assert archive.testzip() is None
dump(O/'DELIVERY.json',{'archive':str(z),'bytes':z.stat().st_size,'sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'zip_integrity_pass':True})
print(json.dumps(acceptance),flush=True);print((O/'DELIVERY.json').read_text(encoding='utf-8'),flush=True)
