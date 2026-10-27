import pathlib,json,datetime,hashlib
import pandas as pd
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase99_six_patient_admission_check_20261020'
url=json.loads((p/'resources/issue2_comments.json').read_text(encoding='utf-8'))[0]['html_url']
m=pd.read_csv(p/'tables/sample_site_information_needed.tsv',sep='\t',encoding='utf-8',keep_default_na=False)
m=m.rename(columns={'patient_id':'legacy_patient_id'})
m['author_corrected_patient_id']=m.legacy_patient_id
m.loc[m['sample'].eq('NCC__BS14'),'author_corrected_patient_id']='NCC__Patient14'
m['author_correction_evidence']=url
m['treatment_corrected']=m['sample'].map({'NCC__BS1':'Pre-chemotherapy','NCC__BS12':'Pre-chemotherapy','NCC__BS14':'Pre-chemotherapy','NCC__BS3':'Post-chemotherapy','NCC__BS4':'Pre-chemotherapy','NCC__BS5':'Post-chemotherapy','NCC__BS7':'Post-chemotherapy'})
m['tumor_status']=m['sample'].map(lambda x:'Recurrent' if x=='NCC__BS7' else 'Primary')
assert m.author_corrected_patient_id.nunique()==7 and m.cells.sum()==806
m.to_csv(p/'tables/author_corrected_crosswalk_pending_anatomy.tsv',sep='\t',encoding='utf-8',index=False)
d=pd.read_csv(p/'tables/new_library_14gene_descriptive_expression.tsv',sep='\t',encoding='utf-8')
d=d.rename(columns={'patient_id':'legacy_patient_id'})
d['author_corrected_patient_id']=d['sample'].map(m.set_index('sample').author_corrected_patient_id)
d.to_csv(p/'tables/author_corrected_library_expression_descriptive_only.tsv',sep='\t',encoding='utf-8',index=False)
report=f"""新增NCC患者纳入骨肺主分析与参考敏感性分析的核查说明
2026年10月20日

结论
目前不能完成“加入后结论是否成立”的验证。已核查7份文库、806个候选细胞；14个目标基因均在每份文库的测量特征中，原始计数对象SHA-256与既有清单一致，逐文库描述性表达已准备。
但实际测序病灶的骨/肺取材部位仍未确认，新增骨肺比较有效患者数仍为0。
此外，作者更正表明，此前登记为6位患者有误：BS14为Patient14；BS4与BS15才是同一患者治疗前后。现有纳入文库没有BS15，因此对应7位来源内患者，而不是6位。该计数不等同于已完成跨研究遗传去重的全局新增独立患者数。
原始主分析及≥20个正常参考细胞敏感性分析的数值保持原版本；未运行新增分组检验，不能把“数值未变”写成“加入新患者后仍然成立”。

问题与待补条件
CAF-ADD-01｜阻断／待补信息
位置：phase93_expanded_admission_20261016/tables/admitted_patients.tsv，全部6条旧患者记录；phase94_NCC_site_verification_20261016/SITE_VERIFICATION.json中的samples。
证据：7份文库取材部位均为Unknown，Bone_Lung_comparison_eligible均为False。今日重查作者仓库，文件树未变；Issue2作者回复仍仅提供Primary/Recurrent，未提供骨/肺解剖分组。
影响：无法构造新增样本的骨肺效应或更新以其为基础的正常参考背景对比。
处理建议：补充BS1、BS3、BS4、BS5、BS7、BS12、BS14逐文库实际测序病灶部位，以及对应临床/作者证据。不能用患者原发肿瘤位置、Primary/Recurrent标签或表达聚类反推本次取材部位。
复核要求：每份文库能追溯到Bone、Lung或Other/Unknown；保留证据来源，不将未知强行二分。

CAF-ADD-02｜高／患者编号必须修订后方可进入新增推断
位置：phase93_expanded_admission_20261016/tables/admitted_patients.tsv中NCC__Patient4行，libraries=NCC__BS14;NCC__BS4；对应新细胞登记表patient_id。
证据：作者在2026-10-15的Issue2回复澄清BS14属于Patient14、为治疗前样本；BS4与BS15才属于Patient4治疗前后配对。
影响：旧表将两位患者错误合并，低估来源内患者数；若沿用合并会错误聚合表达并影响独立样本量。
处理建议：本轮另存带legacy_patient_id及author_corrected_patient_id的工作对照表；原入库基准和原分析结果未覆盖。后续重分析采用有作者依据的新映射。
复核要求：7份现有文库对应7位来源内患者；Patient4仅计BS4的251个候选细胞，Patient14计BS14的26个；总细胞数806不变，BS15不在本次纳入对象中。
作者证据：{url}

后续统计口径
原两个肿瘤来源SCPCP000017与SCPCP000023保留各自患者级主分析。NCC为第三个来源，不能直接把其患者数加入这两个来源。明确取材部位后，先检查NCC两侧是否有独立患者重复，以及与既有来源的患者重叠证据，再决定来源内比较与跨来源汇总。
正常参考细胞阈值沿用每位正常供者至少20个所选参考细胞；该阈值不自动用于筛掉新肿瘤患者。
不能依据期望P值、候选基因表达或现有分群为未知样本指定骨/肺标签。

核查范围
现有7个已纳入计数对象及其哈希、14基因覆盖、旧患者登记、既有部位核查记录、作者最新仓库树、当前Patient metadata、Issue1与Issue2及现有评论。未重新开展全项目QC或新的细胞鉴定；未向作者发送消息。

来源内更正对照
"""
for r in m.itertuples():
 report+=f"{r.sample}：{r.legacy_patient_id} → {r.author_corrected_patient_id}；{r.cells}个细胞；取材部位Unknown。\n"
report+="\n辅助工作表、描述性表达与作者证据均保存在服务器本阶段目录；本报告不是加入新患者后的统计结果。\n"
out=p/'delivery/CAF_new_patients_analysis_prerequisites_20261020.txt';out.write_text(report,encoding='utf-8')
status=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='admission_check_complete_inference_not_run_missing_sequenced_site',old_registered_patients=6,author_corrected_source_patients=7,libraries=7,cells=806,new_valid_Bone_Lung_patients=0,inferential_tests_run=False,baseline_results_unchanged=True,findings=['CAF-ADD-01','CAF-ADD-02'],author_evidence=url,report=str(out),sha256=hashlib.sha256(out.read_bytes()).hexdigest())
(p/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
(root/'CAF_NEW_PATIENT_RECHECK_CURRENT.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(status,ensure_ascii=False,indent=2))
