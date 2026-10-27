from pathlib import Path
import json,os,re
os.environ['MPLBACKEND']='Agg'
import numpy as np,pandas as pd,matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def order(s):return (s.split('_')[0],int(re.findall(r'C(\d+)',s)[0]))
cfg=json.loads((O/'resources/PROGRAMS.json').read_text());s=read(O/'tables/patient_program_statistics.tsv.gz');t=read(O/'tables/posthoc_collagen_structure_statistics.tsv');b=read(O/'tables/subgroup_baseline.tsv').set_index('subgroup_id');g=read(O/'tables/targeted_mechanism_gene_statistics.tsv.gz');sp=read(O/'tables/spatial_patient_statistics.tsv');site=read(O/'tables/spatial_site_statistics.tsv')
assert (O/'MECHANISM_COMPLETE').exists() and (O/'SPATIAL_COMPLETE').exists() and (O/'STRUCTURE_COMPLETE').exists()
allstats=pd.concat([s.assign(analysis='prespecified'),t.assign(analysis='posthoc_structure')],ignore_index=True);ledger=allstats[(allstats.matrix=='corrected')&(allstats.scope=='all')&(allstats.minimum_cells==1)].copy()
def evidence(r):
 if not np.isfinite(r.mean_effect):return '不可评估'
 if pd.notna(r.P) and r.P<.05:return '患者名义正向支持' if r.mean_effect>0 else '患者名义负向'
 if r.N<2:return '描述性正向 无生物学重复检验' if r.mean_effect>0 else '描述性负向 无生物学重复检验'
 return '正向但未达名义阈值' if r.mean_effect>0 else '负向但未达名义阈值'
ledger['goal_evidence']=ledger.apply(evidence,axis=1);ledger['primary_source']=ledger.subgroup_id.map(b.primary_source);ledger['is_primary_source']=ledger.dataset.eq(ledger.primary_source)
for c in ['working_name','functional_interpretation_confidence','lineage_confidence'] if 'lineage_confidence' in b else ['working_name','functional_interpretation_confidence']:ledger[c]=ledger.subgroup_id.map(b[c])
ledger.to_csv(O/'tables/goal_mechanism_evidence_ledger.tsv',sep='\t',index=False,encoding='utf-8-sig')
main=ledger[ledger.is_primary_source].copy();main.to_csv(O/'tables/primary_source_mechanism_ledger.tsv',sep='\t',index=False,encoding='utf-8-sig')
nom=ledger[ledger.P.lt(.05)&ledger.mean_effect.gt(0)];nom.to_csv(O/'tables/nominal_positive_programs_all_sources.tsv',sep='\t',index=False)
# Gene candidates are nominated by actual biological P only, with every P/FDR retained.
genehits=g[(g.matrix=='corrected')&g.P.lt(.05)&g.mean_effect.gt(0)].copy();genehits['program_membership']=[';'.join(k for k,v in cfg['sets'].items() if gene in v) for gene in genehits.gene];genehits['candidate_type']='RNA mechanism candidate; not validated protein gate or intervention target';genehits.to_csv(O/'tables/nominal_positive_mechanism_genes.tsv',sep='\t',index=False)
focus=['F2_C1','F2_C3','F1_originalC10','F3_C1','F4_C8','F4_C10'];sens=allstats[allstats.subgroup_id.isin(focus)&allstats.program.isin(['collagen_biosynthesis','collagen_crosslinking','LOX_family','fibrillar_I_III_V','collagen_IV','cartilage_matrix'])];sens.to_csv(O/'tables/focused_sensitivity_ledger.tsv',sep='\t',index=False)
priorities=[
{'subgroup_id':'F2_C3','priority':'纤维状胶原相关首轮候选','basis':'I/III/V面板探索性正向；交联程序名义正向；基质画像支持','interpretation':'N=7；程序FDR未通过；限定Fibro及每臂5细胞后完整结果保留；先验证胶原定位及实际沉积','mechanism_candidates':'COL1A1/COL1A2为表达方向线索；PCOLCE与LOXL2有名义基因证据；LOX本身未正向支持'},
{'subgroup_id':'F2_C1','priority':'收缩与IV型胶原相关并行候选','basis':'收缩特征明确；COL4A1/COL4A2、PXDN及LOXL2有名义基因证据','interpretation':'区分基底膜/血管相关结构与肿瘤外围纤维状包绕；不直接称蛋壳CAF','mechanism_candidates':'COL4A1;COL4A2;PXDN;LOXL2'},
{'subgroup_id':'F4_C8','priority':'基质组分探索性候选','basis':'I/III/V及软骨相关面板方向为正；原基质组分画像清楚','interpretation':'主要来源有效对比N=1；仅描述性；区分纤维状与软骨样基质','mechanism_candidates':'ASPN;ACAN;COL11A1；其他分子以完整基因表为准'},
{'subgroup_id':'F1_originalC10','priority':'分泌基质探索性候选','basis':'多个机制程序方向为正','interpretation':'主要来源有效N=1，不作患者显著性结论','mechanism_candidates':'保留原SFRP2/PLA2G2A画像，机制分子另按本轮表核查'},
{'subgroup_id':'F3_C1','priority':'保留内部基质画像 暂不列整体交联优先候选','basis':'原F3内LOX/ECM糖蛋白特征成立；本轮对其他8174候选背景的交联/LOX家族效应为负','interpretation':'不同对照问题分别报告；不能用内部标志优势推断跨所有候选的交联优势','mechanism_candidates':'LOX保留为既有内部标志，不视为本轮优先靶点'},
{'subgroup_id':'F4_C10','priority':'成骨矿化侧向候选','basis':'原IBSP/BGLAP矿化相关表达清楚','interpretation':'当前胶原机制并不突出；仅在实际组织矿化参与目标结构时提升研究优先级','mechanism_candidates':'IBSP;BGLAP；不是已确认CAF特异靶点'}]
pr=pd.DataFrame(priorities);pr['original_functional_confidence']=pr.subgroup_id.map(b.functional_interpretation_confidence);pr.to_csv(O/'tables/goal_oriented_candidate_priorities.tsv',sep='\t',index=False,encoding='utf-8-sig')
requirements=pd.DataFrame([
['真实治疗应答','病例/样本ID；治疗前后取材关系；方案及周期；Huvos或坏死率或明确药敏终点','107个单细胞样本和26张空间台账未提供直接终点','定向核查原始补充临床表；若需客户数据，以去标识样本ID联结；不以复发或治疗史代替耐药'],
['独立肿瘤边界','可与切片坐标对应的病理区域标注','26张低/高分辨组织图可用；验收肿瘤mask为0','在独立于待检验胶原RNA的组织学依据上定义肿瘤边缘'],
['胶原沉积及纤维结构','与切片配准的胶原组织学或纤维成像及定量','当前RNA对象没有直接物理沉积/交联测量','分别评估胶原量、排列、连续性和交联；不以RNA分数冒充'],
['干预与药效','匹配候选细胞/基质条件；药物单独、基质干预单独、联合及相应对照','本轮无实际干预数据','联测基质改变、药物进入/作用和肿瘤细胞反应，区分直接细胞毒效应与基质介导效应']
],columns=['endpoint','required_fields','current_status','next_action']);requirements.to_csv(O/'tables/NEXT_DATA_REQUIREMENTS.tsv',sep='\t',index=False,encoding='utf-8-sig')
# Main mechanism landscape, not an arbitrary overall score.
programs=['fibrillar_I_III_V','collagen_IV','cartilage_matrix','collagen_biosynthesis','collagen_assembly','collagen_crosslinking','collagen_degradation','ECM_degradation','LOX_family']
sts=sorted(b.index,key=order);wide=main.pivot(index='subgroup_id',columns='program',values='mean_effect').reindex(index=sts,columns=programs);fig,ax=plt.subplots(figsize=(12,13));im=ax.imshow(wide,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_yticks(range(len(sts)),[z.replace('_original','_') for z in sts]);ax.set_xticks(range(len(programs)),[z.replace('_',' ') for z in programs],rotation=40,ha='right')
for i,st in enumerate(sts):
 for j,p in enumerate(programs):
  r=main[(main.subgroup_id==st)&(main.program==p)]
  if len(r) and r.iloc[0].P<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=12)
fig.colorbar(im,ax=ax,label='Mean paired patient log2(CP10k + 0.1) difference',shrink=.6);ax.set_title('Goal-oriented matrix mechanism landscape | all 33 subgroups',loc='left');fig.text(.02,.015,'* Raw patient P < 0.05; FDR in full ledger. First three columns are posthoc structural panels. N and primary source differ by row.',fontsize=8);fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(O/'figures/subgroup_matrix_mechanisms.png',dpi=170);plt.close(fig)
# Explicit within-patient effects for leading two candidates.
pe=pd.concat([read(O/'tables/patient_program_effects.tsv.gz'),read(O/'tables/posthoc_collagen_structure_patient_effects.tsv')]);show=pe[(pe.dataset=='SCPCP000017')&(pe.matrix=='corrected')&(pe.scope=='all')&(pe.minimum_cells==1)&pe.subgroup_id.isin(['F2_C1','F2_C3'])&pe.program.isin(['fibrillar_I_III_V','collagen_IV','collagen_crosslinking','LOX_family'])].dropna(subset=['effect']);fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True);orderp=['fibrillar_I_III_V','collagen_IV','collagen_crosslinking','LOX_family']
for ax,st in zip(axs,['F2_C1','F2_C3']):
 for j,p in enumerate(orderp):
  v=show[(show.subgroup_id==st)&(show.program==p)].effect.to_numpy();ax.scatter(v,np.full(len(v),j),s=30,alpha=.6);ax.scatter([v.mean()],[j],marker='D',s=65,color='#a84026')
 ax.axvline(0,color='gray',lw=.8);ax.set_yticks(range(4),[z.replace('_',' ') for z in orderp]);ax.set_title(st+' | 7 matched patients');ax.set_xlabel('Patient mean log2(CP10k + 0.1) difference')
fig.suptitle('Collagen structure and crosslinking must be separated',weight='bold');fig.tight_layout();fig.savefig(O/'figures/F2_patient_mechanism_effects.png',dpi=170);plt.close(fig)
# Show major spatial sensitivity rather than only the significant full score.
selprog=['collagen_biosynthesis','collagen_assembly','collagen_crosslinking'];eps=['tumor_local_adjusted','CD8_local_adjusted'];fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True)
for ax,ep in zip(axs,eps):
 for i,p in enumerate(selprog):
  for off,version,col in [(-.1,'full_program','#347c86'),(.1,'minus_author_contributors','#c1874e')]:
   r=sp[(sp.program==p)&(sp.endpoint==ep)&(sp.version==version)].iloc[0];ax.scatter(r.mean_effect,i+off,s=65,color=col,label=version if i==0 else None)
 ax.axvline(0,color='gray',lw=.8);ax.set_yticks(range(3),[x.replace('_',' ') for x in selprog]);ax.set_title(ep.replace('_',' '));ax.set_xlabel('Patient mean partial rank correlation');ax.legend(fontsize=8,frameon=False)
fig.suptitle('Spatial associations depend on component-contributor genes',weight='bold');fig.tight_layout();fig.savefig(O/'figures/spatial_overlap_sensitivity.png',dpi=170);plt.close(fig)
def fmt(x):return 'NA' if pd.isna(x) else f'{x:.6g}'
lines=['骨肉瘤基质相关治疗抵抗研究 首轮实施结果','2026-10-05 phase63','',
'研究结论','本轮已完成目标导向计划冻结、33个内部亚群的机制比较、26张空间切片/19个患者记录的程序定位分析，以及直接治疗应答和独立肿瘤边界的数据条件核查。原8174个细胞、F1–F5、内部编号和全部原置信度保持不变。',
'当前最接近纤维状胶原问题的首轮候选为F2-C3；F2-C1作为收缩及IV型胶原相关并行候选。现有信号属于小样本探索性证据，不将任何亚群标记为已确认耐药执行者。F4-C8、F1-C10保留描述性候选；F3-C1的原内部标志与本轮更广背景下效应分别解释。','',
'主要结果','F2-C3的I/III/V型胶原事后分解面板：N=7，平均效应+1.796，6/7患者正向，P=0.03125，FDR=0.16544；完整交联程序效应+0.673，P=0.046875，FDR=0.20833。均为名义信号，未通过对应FDR阈值。',
'F2-C1的IV型胶原面板：N=7，平均效应+1.356，7/7正向，P=0.015625，FDR=0.11719；I/III/V面板效应+1.134，P=0.0625。完整交联程序效应+0.712，P=0.03125，FDR=0.14583。该结果要求区分基底膜相关成分与外围纤维状胶原。',
'交联集合包含COL4A1/COL4A2和PXDN等。基因层面，F2-C3的PCOLCE效应+1.474（P=0.015625，FDR=0.33453）、LOXL2效应+2.041（P=0.046875，FDR=0.55078）；F2-C1的PXDN效应+3.182、LOXL2效应+2.188（均P=0.015625，FDR=0.33453）。它们是RNA机制候选，不是已经验证的干预靶点。LOX本身在这两个对比中均未显示正向支持。',
'原始与去污染后完整交联方向相近。每臂至少5个细胞时F2-C1/C3均仅N=3；Fibro限定后交联效应分别为+0.408（N=6，P=0.1875）和+0.416（N=7，P=0.15625），均未达到名义阈值。完整敏感性表保留所有P及FDR。',
'F3-C1在原F3内部具有LOX/特定基质组分优势，但本轮对同患者其他8174候选背景的交联效应为-0.86、LOX家族为-0.76（N=2）；当前不据原名称将其列为整体交联首选。F4-C8纤维状及软骨相关面板均正向，但主要来源N=1，仅作描述。','',
'空间主线','完整胶原生成程序与肿瘤成分的组成调整后患者平均相关系数为+0.16078，N=19，16/19正向，P=0.000103，FDR=0.000561；装配程序为+0.15875，18/19正向，P=0.000160，FDR=0.000561。相邻肿瘤成分关联也有支持。',
'该关联不能作为独立空间验证：剔除作者成分模型贡献基因后，胶原生成与局部肿瘤成分的平均相关为-0.00957（P=0.3863），装配为+0.02035（P=0.3437）；原发肿瘤子集的胶原生成关联为+0.05409（N=9，P=0.4414）。完整交联程序与局部肿瘤成分没有名义支持（P=0.2319）。主表不隐藏这些结果。','',
'蛋壳侧线','胶原生成程序与局部CD8 RNA在组成调整后呈弱负相关（平均-0.03984，N=19，15/19负向，P=0.000172，FDR=0.001202）；未经调整P=0.5834，去除成分贡献基因后P=0.6185，原发肿瘤子集P=0.3711。因此保留模型依赖的空间关联，不能称真实免疫排斥或保护。',
'既有COL1A1/LOX/MMP14三基因结果及其原始FDR原样保留。本轮未确认蛋壳环绕结构：组织图像存在，但独立肿瘤边界、物理胶原沉积和纤维组织没有直接测量。','',
'临床终点与下一步','当前107份单细胞样本及26张空间样本整合台账没有直接Huvos/坏死率/药敏终点字段。治疗史、复发、转移和生存状态均未替代耐药标签；本轮没有生成耐药P值。',
'下一步优先补齐能够区分肿瘤周围纤维状胶原与基底膜/血管结构的独立组织学定位，并定向核查原始临床补充表中的真实治疗应答及样本联结。F2-C3与F2-C1分别作为纤维状胶原和收缩/IV型胶原候选推进；若原亚群内部无法定位相应细胞，再局部细分，不先重做全面分类。',
'实验设计应分别测定胶原沉积/纤维组织、药物进入或作用、肿瘤细胞药物反应，并设置药物单独、基质干预单独、联合和匹配对照。机制分子、检测标志与实际可用干预工具分别核实。','',
'方法','冻结7个完整主机制程序、3个完整辅助程序和LOX家族5基因参考面板，机制基因并集658个。使用原始与DecontX后完整文库归一化，计算同来源/患者/部位的目标亚群减其他8174候选细胞的基因log2(CP10k+0.1)差，程序为可测基因效应均值，部位等权后患者等权。每臂至少1细胞为主分析，至少5细胞及Fibro限定为敏感性。来源分开分析，缺测不记零，保留来源特征覆盖。',
'患者双侧符号翻转检验在N<=15精确计算，较大N使用19999次Monte Carlo加观察值。基因与程序分开定义检验家族；程序按来源/矩阵/范围/细胞阈值/主辅角色分家族；原始P<0.05用于名义筛选，FDR完整报告。N<2或全零不计算P，N<10定位探索性。三类结构面板在主结果后追加，独立记录事后检验家族，原主检验不改动。',
'空间复用原始26切片、53309 spots、19患者记录及既有邻接图；每切片计算程序与成分的秩相关/偏秩相关，事件、部位、患者逐级等权。所有推断N为患者，不以spot作N。作者成分为带符号IC值，不解释为细胞比例。原发及其他部位子集分别报告；各版本/终点/机制家族内BH校正。基因重叠剔除为诊断性敏感性，保留定义变化及覆盖，不能与完整程序混称。',
'完整全转录组发现沿用phase62，对应GSEA及原检验家族另表保留；本轮658基因检验是目标机制跟进，不冒称新一次全转录组发现。','',
'证据范围','目前结果定位到RNA程序、来源内患者效应和模型依赖空间关联。独立肿瘤边界、物理胶原与交联结构、真实药效或干预仍需数据。跨队列患者对应未完全解决，不宣称新增独立患者重复验证。所有正负、不可评估及关键敏感性结果保留。',
'胶原IV/PXDN的基底膜相关生物学依据用于解释RNA候选，不替代本项目组织定位。参考：Reactome R-HSA-2243919 https://reactome.org/content/detail/R-HSA-2243919 ；Peroxidasin forms sulfilimine chemical bonds using hypohalous acids in tissue genesis https://pubmed.ncbi.nlm.nih.gov/22842973/ 。骨肉瘤ECM与药物反应的模型依据：Cortini等2023 https://pmc.ncbi.nlm.nih.gov/articles/PMC9954237/ 。','',
'交付索引','计划：CAF基质耐药目标导向执行计划_20261005.txt。候选：goal_oriented_candidate_priorities.tsv；全证据：goal_mechanism_evidence_ledger.tsv。患者基因、程序、结构分解及空间完整结果均在tables目录。空间地图26张在figures/spatial_maps；每张均为RNA叠加，无虚构肿瘤边界。下一步必要数据见NEXT_DATA_REQUIREMENTS.tsv。']
(O/'delivery/CAF基质耐药目标导向首轮结果_20261005.txt').write_text('\n'.join(lines),encoding='utf-8')
facts={'cells':8174,'subgroups':33,'mechanism_genes':658,'programs':11,'posthoc_panels':3,'spatial_slices':26,'spatial_patients':19,'spatial_spots':53309,'direct_response_endpoint_available':False,'accepted_tumor_masks':0,'goal_achieved':'first goal-oriented computational round only','true_resistance_mechanism_validated':False,'eggshell_confirmed':False}
(O/'delivery/REPORT_FACTS.json').write_text(json.dumps(facts,indent=2));(O/'INTEGRATION_COMPLETE').write_text('SUCCESS\n');print('INTEGRATION_COMPLETE',len(ledger),len(genehits),flush=True)
