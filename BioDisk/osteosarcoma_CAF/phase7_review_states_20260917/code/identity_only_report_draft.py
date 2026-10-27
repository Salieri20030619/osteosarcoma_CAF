"""SYNTAX_CHECK: python -m py_compile code/freeze_review_report.py
RUN_COMMAND: python code/freeze_review_report.py
EXPECTED_SUCCESS_MARKER: PAUSED_AFTER_REVIEW_FOR_INTERIM_REPORT
"""
from pathlib import Path
import json,datetime,hashlib
import pandas as pd,numpy as np
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=r/'phase7_review_states_20260917';old=r/'phase6_nonpdx_cnv_20260916';delivery=p/'delivery';delivery.mkdir(exist_ok=True)
assert (p/'identity/COMPLETE').read_text().strip()=='PHASE6_IDENTITY_INTEGRATION_COMPLETE'
q=pd.read_csv(p/'identity/tool_reference_QC.tsv',sep='\t').fillna('');s=pd.read_csv(p/'identity/sample_identity_summary.tsv',sep='\t');before=pd.read_csv(old/'identity/sample_identity_summary.tsv',sep='\t');qc=pd.read_csv(r/'phase5_nonpdx_qc_20260916/results/sample_summary.tsv',sep='\t')
allcells=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t');strict=pd.read_csv(p/'identity/caf_strict_cells.tsv',sep='\t');ext=pd.read_csv(p/'identity/caf_extended_cells.tsv',sep='\t')
assert len(allcells)==125212 and allcells.cell_id.is_unique
assert len(strict)==int(s.caf_strict.sum())==int(allcells.caf_strict.sum())
assert len(ext)==int(s.caf_extended.sum())==int(allcells.caf_extended.sum())
assert set(strict.cell_id)<=set(ext.cell_id)
assert set(pd.read_csv(old/'identity/caf_strict_cells.tsv',sep='\t').cell_id)<=set(strict.cell_id)
checks=[]
for x in s.to_dict('records'):
 sm=x['sample'];d=allcells[allcells['sample']==sm];assert len(d)==int(x['retained'])
 assert int(d.caf_strict.sum())==x['caf_strict'] and int(d.caf_extended.sum())==x['caf_extended']
 for row in q[(q['sample']==sm)&(q.state=='complete')].to_dict('records'):
  metrics=json.loads(row['lineage_checks']);expected=bool(metrics) and all(v['evaluable']>=10 and v['normal_fraction']>=.9 for v in metrics.values()) and float(row['heldout_evaluable'])>=10 and float(row['heldout_normal_fraction'])>=.9
  assert bool(row['reference_check_pass'])==bool(expected),(sm,row['method'])
 if d.caf_strict.any():
  z=d[d.caf_strict];assert z.reference_QC_pass.all() and (z.n_genes>=500).all() and (z.det_Fib_review>=2).all() and (z.lineage_reviewed=='FibLike').all() and (z.cnv_concordance=='both_low').all() and (~z.CTA_context).all() and (z.cluster_high_fraction<.2).all()
 checks.append({'sample':sm,'cell_accounting':'PASS','reference_gate_recalculation':'PASS','strict_criteria':'PASS','strict_cells':int(d.caf_strict.sum()),'extended_cells':int(d.caf_extended.sum())})
# Cross-check the identity lists against an independent recalculation, not only completion flags.
(p/'review/final_identity_audit.json').write_text(json.dumps(checks,indent=2))
compare=before[['sample','reference_QC_pass','caf_strict','caf_extended']].merge(s[['sample','reference_QC_pass','caf_strict','caf_extended','selected_secondary_method']],on='sample',suffixes=('_before','_after'))
compare['strict_added']=compare.caf_strict_after-compare.caf_strict_before;compare['extended_change']=compare.caf_extended_after-compare.caf_extended_before
compare.to_csv(delivery/'复核前后候选集合对照.tsv',sep='\t',index=False)
s.to_csv(delivery/'样本身份与CAF数量.tsv',sep='\t',index=False);q.to_csv(delivery/'工具参考检查完整结果.tsv',sep='\t',index=False)
clinical=pd.read_csv(r/'phase4_download_20260914/clinical_source_audit.tsv',sep='\t',encoding='utf-8-sig').fillna('')
clinical.to_csv(delivery/'临床分组依据与缺失字段.tsv',sep='\t',index=False)
now=datetime.datetime.now(datetime.timezone.utc).astimezone(datetime.timezone(datetime.timedelta(hours=8)))
strict_samples=s[s.caf_strict>0];patients=strict.patient_id.unique().tolist();lung=strict[strict.site=='Lung'];nl=int(lung.patient_id.nunique());nb=int(strict.loc[strict.site!='Lung','patient_id'].nunique());npass=int(s.reference_QC_pass.sum())
raw_discord=int(before.CNV_discordant_or_ambiguous.sum());review_discord=int(s.CNV_discordant_or_ambiguous.sum());pending=int(s.Identity_pending.sum())
failed_samples=['GSE299023__OS.C.21.BX.OCT','GSE299023__OS.C.32.BX','GSE299023__OS.C.33.BX','GSE152048__BC3','GSE152048__BC6']
review_decisions=[]
for sm in failed_samples:
 row=s[s['sample']==sm].iloc[0];methods=q[(q['sample']==sm)&(q.state!='not_run')]
 issues=[]
 for z in methods.to_dict('records'):
  if z['state']!='complete':issues.append(z['method']+'运行失败');continue
  m=json.loads(z['lineage_checks'])
  for role,values in m.items():
   if values['evaluable']<10:issues.append(f"{z['method']} {role}仅{values['evaluable']}个可评估细胞")
   elif values['normal_fraction']<.9:issues.append(f"{z['method']} {role}正常比例{values['normal_fraction']:.1%}")
 decision='通过预设参考门槛，采用'+row.selected_secondary_method+'与inferCNV联合判定' if row.reference_QC_pass else '未通过双方法参考门槛，不纳入严格CAF集合'
 review_decisions.append({'sample':sm,'decision':decision,'strict_CAF':int(row.caf_strict),'extended_CAF':int(row.caf_extended),'issues_retained':'; '.join(issues)})
pd.DataFrame(review_decisions).to_csv(delivery/'五个重点样本复核结论.tsv',sep='\t',index=False)
lines=[]
def add(x=''):lines.append(x)
def table(headers,rows):
 add('| '+' | '.join(headers)+' |');add('| '+' | '.join(['---']*len(headers))+' |')
 for row in rows:add('| '+' | '.join(str(x).replace('|','／').replace('\n',' ') for x in row)+' |')
 add()
add('# 骨肉瘤肺转移相关CAF项目中期报告');add();add('报告日期　'+now.strftime('%Y年%m月%d日'));add('阶段范围　数据接管与质控　细胞谱系复核　CNV交叉鉴定　候选集合冻结');add()
add('## 核心进展');add()
add(f'项目已完成新队列17个组织的逐样本质控、细胞谱系复核及CNV交叉鉴定，并完成5个重点样本的参考证据审查。复核后，新队列保留 **{len(strict):,} 个严格CAF候选**和 **{len(ext):,} 个扩展候选**；扩展集合包含严格集合。严格候选分布于 **{len(strict_samples)} 个组织、{len(patients)} 个队列内已标识患者**，为后续CAF状态分析提供了可追溯的表达对象和细胞清单。');add()
add(f'新队列共有 **{npass}/17 个组织通过双方法参考检查**。此前GSE270231队列的188个严格候选、345个扩展候选保持原结果，作为后续独立队列观察的基础。本次已完成的工作止于身份复核和候选集合冻结；CAF再聚类、成骨／矿化程序比较、差异基因和分选marker分析尚未启动，项目按要求暂停在中期报告节点。');add()
add('## 数据和计数质量已完成核验');add()
add('新队列由GSE299023的6个患者组织及GSE152048的11个患者组织构成，来源核查支持直接患者取材且不属于PDX。29个输入文件全部完成下载，合计1,182,489,788字节，逐文件SHA256核对和gzip完整性读取通过；17个矩阵的维度、非负整数计数、条形码唯一性及基因标识核验通过。GSE152048的多个GSM记录按11个组织矩阵管理，没有按技术记录数增加患者数。');add()
table(['队列','组织或样本','输入条形码','本地保留细胞或核','当前用途'],[['GSE299023','6个组织',int(qc.loc[qc.dataset=='GSE299023','input_cells'].sum()),int(qc.loc[qc.dataset=='GSE299023','retained'].sum()),'新增患者骨组织候选集合'],['GSE152048','11个组织',int(qc.loc[qc.dataset=='GSE152048','input_cells'].sum()),int(qc.loc[qc.dataset=='GSE152048','retained'].sum()),'新增骨与肺病灶候选集合'],['GSE270231','既有7个样本','本次复用已处理对象','28,542','后续跨队列观察']])
add('新队列共160,541个输入条形码；自适应QC后保留131,609个，排除6,397个通过自动阈值检查的双细胞候选，最终保留125,212个细胞并形成263个初始簇。OS.C.32.BX、BC10和BC11的Scrublet自动阈值列入复核，原始分数和预测记录完整保留，没有按固定比例补删细胞。');add()
add('## 联合证据形成逐样本CAF候选集合');add()
add(f'新队列的34项基础CNV任务已全部退出，其中32项成功、2项失败。3个新增SCEVAN核查任务与此前1个SCEVAN补充任务均保留实际运行结果。严格候选数从基础整合的2,047个变为 **{len(strict):,} 个**，扩展候选数从3,167个变为 **{len(ext):,} 个**；变化逐一对应到参考检查、方法选择和细胞清单，未覆盖任何原始工具判定。');add()
rows=[]
for x in s.to_dict('records'):
 rows.append([x['sample'].split('__')[0],x['short'],'肺' if x['site']=='Lung' else '骨',f"{int(x['retained']):,}",'通过' if x['reference_QC_pass'] else '待定',int(x['caf_strict']),int(x['caf_extended'])])
table(['队列','组织编号','取材','保留细胞','双方法参考','严格CAF','扩展CAF'],rows)
add('表中细胞数是身份筛选的描述性结果，不作为独立生物学样本量。严格候选为肿瘤背景下获得成纤维谱系及双方法非恶性证据支持的细胞；扩展集合用于后续敏感性分析。严格集合为零表示当前规则未获得满足全部条件的候选，不等于该组织在生物学上不存在CAF。');add()
add('## 重点复核明确了方法适用范围');add()
add('本轮复核保持每个预设参考谱系至少10个可评估细胞、正常分类比例至少90%的门槛，并要求全部可评估留出参考的正常比例至少90%。补充SCEVAN采用同样本参考，在可行时平衡免疫与内皮数量。方法选择依据参考质量和适用性，没有依据CAF数量挑选结果。');add()
add('复核同时修复了参考汇总逻辑：SCEVAN未用于平衡基线的原始参考细胞仍按冻结时的谱系纳入检查，不能因改标为check_ref而从相应参考检查中消失；联合留出参考的总体正常比例也被明确纳入通过条件。该修复防止只检查被算法选入基线的细胞而高估参考可靠性。');add()
for x in review_decisions:
 add('### '+x['sample'].replace('__',' '));add();add(x['decision']+'。严格候选 '+str(x['strict_CAF'])+' 个，扩展候选 '+str(x['extended_CAF'])+' 个。');add('保留的异常或不足证据：'+(x['issues_retained'] or '本轮实际运行方法未发现低于预设门槛的参考指标')+'。');add()
add(f'复核后的候选集合保留全部125,212个细胞的证据记录，其中双方法不一致或模糊记录为 **{review_discord:,} 个**，身份待定记录为 **{pending:,} 个**。这些记录继续单列，不强行归为CAF或恶性细胞。基础整合的不一致或模糊记录为{raw_discord:,}个；由于部分样本改用通过参考检查的补充方法，两个版本的计数用于审计结果变化，不用于评价方法优劣。');add()
add('## 肺转移组织的可用CAF样本范围已经明确');add()
add(f'GSE152048的BC10和BC17均有明确肺转移取材依据。复核后的严格CAF来自 **{nl} 位肺转移患者**和 **{nb} 位骨病灶患者**；BC10贡献{int(s.loc[s["sample"]=="GSE152048__BC10","caf_strict"].iloc[0])}个严格候选，BC17贡献{int(s.loc[s["sample"]=="GSE152048__BC17","caf_strict"].iloc[0])}个。该分布明确了后续病灶比较的实际数据基础。');add()
add('GSE299023的6个组织保留原始骨取材信息；缺少肺转移史的骨组织未补推为A组或B1组。GSE152048的原发与复发性质、治疗记录分别保存，BC17额外接受抗PD-L1治疗的注释保留在临床上下文表。既有GSE270231严格集合涉及S0058、S0059、SC069和SC072；前两者同属NCHS-009，后两者患者映射未确认，不能据此计算4位独立患者。');add()
add('## 分析方法');add()
add('每个组织独立计算UMI、检出基因数和线粒体比例。检出基因下限为max(200, floor(10^(median(log10(n_genes+1))−3×MAD)−1))，UMI下限为相同形式且最低300；MAD采用未缩放绝对中位差。线粒体比例上限为min(20%, max(5%, 中位数+3×MAD))。高复杂度仅标记，Scrublet预期双细胞率为0.06并固定随机种子。对低细胞数、缺失阈值或异常预测比例保留复核状态。');add()
add('表达按每细胞10,000 UMI归一化并log1p变换，原始整数counts另行保留。使用高变基因、PCA、近邻图和Leiden 0.4／0.8／1.2分辨率进行逐样本聚类，以0.8作为谱系复核入口。30项具名marker复核与竞争谱系表达记录并行保存，区分成纤维样、成骨样、髓系／破骨、内皮、壁细胞、红系、肌源性及待定群体。');add()
add('新队列CNV共输入97,017个细胞，使用完整RNA计数和按稳定Ensembl ID匹配的GRCh38位置。Ensembl 98与三个36,601基因矩阵的ID及符号完全吻合，其余矩阵稳定ID匹配率超过99%；未映射条目逐样本保留。inferCNV使用cutoff 0.1、101基因平滑窗、denoise开启、HMM关闭；CopyKAT使用ngene.chr 5、win.size 25、KS.cut 0.1、LOW.DR 0.05和请求UP.DR 0.1，工具内部实际调整记录在日志中。补充SCEVAN使用官方预处理和分类步骤、ngenes_chr 5、perc_genes 0.1、beta_vega 0.5、4线程，关闭自动成纤维正常签名选择，不强制其余细胞为肿瘤。');add()
add('参考来自同样本已复核免疫／内皮谱系，要求至少500个检出基因、750 UMI、至少3个本谱系标志，本谱系模块UMI占比至少0.5%，竞争间质模块合计不超过0.2%，并控制其他竞争谱系和双细胞分数。保留冻结的参考与留出细胞对应关系。严格CAF要求成纤维谱系、至少2个成纤维模块基因、通过参考检查的两种方法一致支持低CNV／非恶性、无强竞争或癌睾抗原程序、至少500个检出基因，且所在簇联合可评估细胞中的双方法高CNV比例低于20%。');add()
add('本阶段为细胞身份鉴定与质量审查，报告计数、比例和方法一致性，不进行患者组间推断检验。原始失败记录、未评估结果、参考指标和逐细胞判定均进入完整审计表。');add()
add('## 当前证据边界');add()
add('当前数据支持可追溯的CAF候选身份筛选，尚不支持已发现肺转移特异新亚群、成骨／矿化机制或可直接使用的临床分选panel。新队列肺转移严格CAF仅由单一独立患者提供时，后续骨／肺对比只能报告探索性方向和效应量；细胞数不能替代患者重复。缺少随访和肺转移史的骨组织不能承担A与B1比较。');add()
add('表达推断CNV的低信号不等同于DNA二倍体确证，参考检查通过也不证明对所有间质谱系具有相同特异性。低复杂度、参考覆盖不足、工具不一致和环境RNA相关跨谱系表达分别保留。当前进行的是跨谱系表达筛查，没有宣称完成环境RNA校正。对于不能满足既定检查条件的组织，维持严格集合排除，不降低阈值补入。');add()
add('不同队列的建库、取材与治疗背景需在下一阶段分开检查。GSE270231中的重复患者和未确认映射会限制独立验证的有效重复数。当前CAF候选集合的建立，为后续分析提供了事实基础；这些边界决定比较方式和结论强度。');add()
add('## 暂停节点与下一阶段计划');add()
add('本次复核、集合更新和中期报告完成后，项目暂停。已准备的analyze_states.py仅作为后续实施草案，未执行，不计入分析成果。下一阶段拟依次开展严格CAF集合内的状态识别与患者支持检查、myCAF细分、成骨／矿化／钙调节程序与核心基因效应量、扩展集合及QC敏感性分析，以及分队列验证；候选表面marker筛选随后衔接。');add()
add('## 交付与追溯文件');add()
table(['文件','内容'],[['中期报告.md','本报告正文'],['样本身份与CAF数量.tsv','17个组织的最终身份与候选数量'],['五个重点样本复核结论.tsv','5个重点样本的处理结果和保留问题'],['复核前后候选集合对照.tsv','基础整合与本轮复核的逐样本变化'],['工具参考检查完整结果.tsv','各方法、各参考谱系与留出细胞的完整检查'],['临床分组依据与缺失字段.tsv','取材、治疗、肺转移标签及待补字段']])
add('服务器目录：`/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917/`。最终候选清单位于`identity/`，可复用表达对象位于`identity/objects/`；核查结果位于`review/final_identity_audit.json`。阶段6结果完整保留，不覆盖。');add()
add('## 公开数据来源');add()
add('- GSE299023　https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299023');add('- GSE152048　https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152048');add('- GSE270231　https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270231');add('- 新队列来源与临床核查所用原文　https://pmc.ncbi.nlm.nih.gov/articles/PMC12248827/　https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/');add()
(delivery/'中期报告.md').write_text('\n'.join(lines),encoding='utf-8')
summary={'state':'PAUSED_AFTER_REVIEW_FOR_INTERIM_REPORT','paused_at':now.isoformat(),'new_cohort_strict_CAF_candidates':len(strict),'new_cohort_extended_CAF_candidates':len(ext),'reference_pass_samples':npass,'strict_CAF_samples':len(strict_samples),'strict_CAF_dataset_patient_ids':len(patients),'lung_patients':nl,'bone_patients':nb,'remaining_reference_fail_samples':s.loc[~s.reference_QC_pass,'sample'].tolist(),'analysis_started':False,'report':'delivery/中期报告.md','report_format':'Markdown','audit':'review/final_identity_audit.json'}
(p/'STATUS.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(p/'PAUSED_AFTER_REVIEW_FOR_INTERIM_REPORT').write_text(now.isoformat()+'\n')
project=json.loads((r/'PROJECT_STATE.json').read_text());project.update({'updated_at':now.isoformat(),'current_stage':p.name,'dynamic_status_file':p.name+'/STATUS.json','status':'PAUSED_AFTER_REVIEW_FOR_INTERIM_REPORT','next_model_review':'resume downstream CAF analysis only after further user instruction','interim_report':p.name+'/delivery/中期报告.md'})
(r/'PROJECT_STATE.json').write_text(json.dumps(project,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
