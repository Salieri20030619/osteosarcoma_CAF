import pathlib,json,hashlib,shutil
import pandas as pd,numpy as np
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase108_each_source_significant_20261020'
for d in ['code','tables','resources']:(P/d).mkdir(parents=True,exist_ok=True)
f1=R/'phase107_third_NCI_normal_background_20261020/tables/patient_permutation_allgenes_followup.tsv.gz'
f2=R/'phase105_original_rule_NCI_20261020/tables/original_two_source_rule_reproduced.tsv.gz'
a=pd.read_csv(f1,sep='\t');b=pd.read_csv(f2,sep='\t')
z=a.merge(b[['feature','P_17','P_23','FDR_17','FDR_23','mean_lung_minus_bone_17','mean_lung_minus_bone_23']],on='feature',validate='one_to_one')
z['original_each_nominal']=z.P_17.lt(.05)&z.P_23.lt(.05)
z['final_selected']=z.remain_after_both_corrections&z.original_each_nominal
z['final_background_both_FDR']=z.final_selected&z.FDR_pooled.lt(.05)&z.FDR_NCI.lt(.05)
# Reconstruct the combined rule independently from numerical fields.
sign=lambda x:np.where(abs(x)<=1e-12,0,np.sign(x))
s=sign(z.mean_lung_minus_bone)
same=(s!=0)&(s==sign(z.mean_lung_minus_bone_17))&(s==sign(z.mean_lung_minus_bone_23))&(s==sign(z.third_effect))&(s==sign(z.effect_pooled))&(s==sign(z.effect_NCI))
rule=same & z.P.lt(.05)&z.third_P.lt(.05)&z.P_17.lt(.05)&z.P_23.lt(.05)&z.P_pooled.lt(.05)&z.P_NCI.lt(.05)&z.normal_evaluable
assert np.array_equal(rule,z.final_selected)
assert set(z.loc[z.final_selected,'gene'])=={'PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1'}
z.to_csv(P/'tables/all_genes_each_source_gate.tsv.gz',sep='\t',index=False,encoding='utf-8',na_rep='NA')
cols=['feature','gene','P_17','FDR_17','P_23','FDR_23','mean_lung_minus_bone','P','FDR','third_effect','third_P','third_FDR','effect_pooled','P_pooled','FDR_pooled','effect_NCI','P_NCI','FDR_NCI']
z[z.final_selected][cols].to_csv(P/'tables/final_six_genes.tsv',sep='\t',index=False,encoding='utf-8')
steps=[('此前全基因背景筛选',z.remain_after_both_corrections),('增加原两队列各自P<0.05',z.final_selected),('其中两项背景对比均全基因FDR<0.05',z.final_background_both_FDR)]
summary=[dict(step=k,count=int(m.sum()),genes=z.loc[m,'gene'].tolist() if int(m.sum())<20 else []) for k,m in steps]
protocol=dict(user_request='各自显著加上再看',extra_gate='Original SCPCP000017 and SCPCP000023 patient-level permutation P<0.05 each; original uncorrected tumor contrast, not reference-corrected per-source tests',
parent='phase107 whole-gene normal-background followup of phase106 selection',
candidate_restriction=False,all_genes_evaluated=len(z),preserve_all_P_and_FDR=True,FDR_scope='Previously computed full-gene families unchanged, not recomputed within6 selected genes',
all_four_source_and_background_FDR='Not requested; original two sources did not pass genomewide FDR, so no claim of all-comparison FDR significance')
(P/'resources/PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
lines=['加上原两队列各自显著后的筛选结果','2026-10-20 | phase108','',
'在全部20,315个共同可评估肿瘤基因上筛选，没有预先限定原14个。',
'延续phase107：合并比较、NCI对骨比较均原始P<0.05、原两个队列与上述比较同向，且两项扣除正常骨/肺背景后仍同向P<0.05。',
'本次新增：原SCPCP000017与SCPCP000023各自原始P<0.05，沿用最初患者标签置换检验。',
'结果：扣正常背景前3,956个缩至8个；其中7个有正常参考；完成背景筛选后的2,070个缩至6个。',
'最终6个：PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1。',
'其中5个在合并与NCI对骨的背景对比中均全基因FDR<0.05：TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1。PNKD在NCI背景对比FDR=0.071526，仅达到该项名义显著。',
'这里的5个不表示原两个队列的全基因FDR也显著；原两队列采用的是原始P门槛。未在6个结果内重新计算FDR。',
'',
'基因\t原17 P\t原23 P\t合并扣背景P\tNCI扣背景P\t合并扣背景FDR\tNCI扣背景FDR']
for r in z[z.final_selected].itertuples():lines.append('\t'.join([r.gene]+[format(getattr(r,k),'.6g') for k in ['P_17','P_23','P_pooled','P_NCI','FDR_pooled','FDR_NCI']]))
lines+=['','QA：从原始数值独立重建所有P值和方向条件，与逐层筛选结果一致；原P/FDR未改动。各比较共享部分患者及参考，仍属于一致性支持，不作独立验证声明。']
report=P/'each_source_significance_brief_20261020.txt';report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
d=dict(status='complete',phase=str(P),report=str(report),report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),before_normal_after_source_gate=8,after_normal_after_source_gate=6,background_both_FDR=5,selected=z.loc[z.final_selected,'gene'].tolist(),QA='PASS',no_original14_restriction=True)
(P/'FINAL_STATUS.json').write_text(json.dumps(d,indent=2),encoding='utf-8')
(R/'CAF_EACH_SOURCE_SIGNIFICANCE_CURRENT.json').write_text(json.dumps(d,indent=2),encoding='utf-8')
snap=P/'resources/previous_current';snap.mkdir(exist_ok=True);f=R/'CAF_RESUME_STATUS_CURRENT.json'
if not (snap/f.name).exists():shutil.copyfile(f,snap/f.name)
old=json.loads(f.read_text());old.update(authoritative_current_file=str(R/'CAF_EACH_SOURCE_SIGNIFICANCE_CURRENT.json'),current_phase=str(P),status='each_original_source_significance_gate_complete',active_task='Original-source P gates added to full-gene normal-background pipeline;6 retained',latest_analysis_report=str(report))
f.write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(d))
