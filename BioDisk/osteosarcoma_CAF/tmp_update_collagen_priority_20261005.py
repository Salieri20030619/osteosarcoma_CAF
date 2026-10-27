from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'priority_collagen_matrix_20261005';O.mkdir(exist_ok=True)
S=B/'phase66_doublet_repair_DSC1_20261005'
sets=json.loads((S/'resources/PROGRAMS.json').read_text())['sets']
roles=['collagen_biosynthesis','collagen_assembly','collagen_crosslinking','collagen_components','ECM_organization']
genes=set().union(*(set(sets[k]) for k in roles))
g=pd.read_csv(S/'tables/gene_before_after_comparison.tsv.gz',sep='\t')
summary=pd.read_csv(B/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv',sep='\t').set_index('subgroup_id')
z=g[g.gene.isin(genes)].copy();z['primary_source']=z.subgroup_id.map(summary.primary_source);z['is_primary_source']=z.dataset.eq(z.primary_source)
z['program_membership']=z.gene.map(lambda x:';'.join(k for k in roles if x in sets[k]))
z['also_in_degradation_program']=z.gene.map(lambda x:x in sets['collagen_degradation'] or x in sets['ECM_degradation'])
z['evidence_scope']='patient_site_pseudobulk_target_vs_other_candidate_cells_within_source'
z['up_nominal_all']=z.mean_effect_all.gt(0)&z.P_all.lt(.05)
z['up_nominal_excluded']=z.mean_effect_excluded.gt(0)&z.P_excluded.lt(.05)
z.to_csv(O/'collagen_matrix_gene_evidence_all.tsv.gz',sep='\t',index=False)
q=z[z.matrix.eq('corrected')&z.is_primary_source&z.up_nominal_excluded].copy()
q['evidence_level']=q.FDR_excluded.lt(.05).map({True:'FDR_significant',False:'nominal_exploratory'})
q=q.sort_values(['subgroup_id','mean_effect_excluded'],ascending=[True,False])
q.to_csv(O/'collagen_matrix_upregulated_primary.tsv',sep='\t',index=False)
focus=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL5A3','LOX','LOXL1','LOXL2','LOXL3','LOXL4','PCOLCE','BMP1','PLOD2','P4HA1','SERPINH1','FN1','DCN','LUM','VCAN','COL4A1','COL4A2','PXDN']
watch=z[z.matrix.eq('corrected')&z.is_primary_source&z.gene.isin(focus)]
watch.to_csv(O/'focus_genes_all_primary_subgroups_including_nonsignificant.tsv',sep='\t',index=False)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
plan={'at':now,'user_instruction':'客户表示先看胶原纤维联合和基质沉淀的基因上调，耐药抑制如果一时证明不了先搁置','current_priority':'胶原纤维形成 装配与交联 以及基质沉积相关基因上调','interpretation_of_customer_wording':'将胶原纤维联合落实为装配与交联两个相关环节；基质沉淀在报告中规范表述为基质沉积。','active_work':['逐亚群整理胶原组成 合成加工 装配与交联 基质组织相关基因的正向效应','报告效应量 原始P FDR 有效患者数及同向患者比例','结合原始与去污染 双联体排除及细胞来源证据排列候选优先级','将胶原合成/装配/交联与降解重塑分别呈现，基因上调不等同实际净沉积'],
'deferred_work':['耐药机制验证','抑制相关机制验证'],'deferred_policy':'暂缓，非取消；当前胶原及基质主线的候选筛选和交付不以证明耐药或抑制为前提。历史阴性及已完成结果保留。','statistics':'原始P<0.05及正效应筛选名义上调；FDR完整保留，不作唯一门槛；有效独立患者少于10例按探索性呈现。无有效重复不补造P。','current_gene_evidence':'复用Phase66完整基因统计；主要来源由既有亚群细胞数预定；主表与排除26个模型阳性后的敏感性证据并列保存；本次未重算模型或缩小FDR检验家族。','comparison':'目标亚群与同来源 同患者 同部位其他候选细胞比较；各部位等权后按患者等权汇总log2(CP10k+0.1)差；不等同治疗前后差异。','screening_summary':q.groupby('subgroup_id').size().to_dict(),'all_gene_evidence':str(O/'collagen_matrix_gene_evidence_all.tsv.gz'),'nominal_upregulated_primary':str(O/'collagen_matrix_upregulated_primary.tsv'),'negative_and_unevaluable_focus_genes':str(O/'focus_genes_all_primary_subgroups_including_nonsignificant.tsv'),'input_sha256':hashlib.sha256((S/'tables/gene_before_after_comparison.tsv.gz').read_bytes()).hexdigest()}
(O/'CUSTOMER_PRIORITY_UPDATE.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
text=['CAF项目客户优先级调整及已有基因上调证据','2026-10-05','',
'当前主线为胶原纤维形成、装配与交联，以及基质沉积相关基因上调。耐药及抑制相关机制验证暂缓，不作为当前阶段完成条件。','',
'已完成的证据整理',
'依据Phase66现有统计，在去污染、主要来源、排除26个模型阳性细胞的敏感性视图中，以正效应及原始P<0.05提取胶原和基质组织相关名义上调。原始主表、原始表达、全部来源、非显著及不可评估结果均保留，未重新计算或改变FDR家族。',
'主要来源筛选得到：'+str(plan['screening_summary'])+'。这些是基因与亚群对应关系条数，不是独立样本量。','',
'F2-C3优先核查示例']
for gene in ['COL3A1','COL5A1','PCOLCE','LOXL2','FN1','VCAN']:
 r=q[(q.subgroup_id=='F2_C3')&(q.gene==gene)].iloc[0]
 text.append(f"{gene}：N={int(r.N_excluded)}，log2(CP10k+0.1)平均差={r.mean_effect_excluded:.6f}，同向{int(r.positive_patients_excluded)}/{int(r.N_excluded)}，P={r.P_excluded:.6f}，FDR={r.FDR_excluded:.6f}；名义显著，探索性。")
text+=['','上调的比较对象为同患者同部位其他候选细胞，不是治疗前后对照。胶原和基质组织集合同时包含结构、酶及调控成员，不能将全部入选者称为沉积促进因子。降解重塑证据单列；实际胶原沉积与纤维结构作为后续独立测量终点。','',
'后续顺序：先形成基因级证据与表达图，再按亚群、来源及患者一致性收敛验证候选；耐药和抑制机制留待后续重启。']
(O/'客户主线调整与胶原基质上调证据_20261005.txt').write_text('\n'.join(text)+'\n')
backup=O/'previous_progress';backup.mkdir(exist_ok=True)
for name in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
 p=B/name
 if not (backup/name).exists():shutil.copy2(p,backup/name)
 d=json.loads(p.read_text());d['previous_user_priority']=d.get('user_priority')
 d.update({'user_priority':plan['current_priority'],'customer_priority_update':str(O/'CUSTOMER_PRIORITY_UPDATE.json'),'active_research_focus':plan['active_work'],'deferred_objectives':plan['deferred_work'],'deferred_objectives_policy':plan['deferred_policy'],'next_action':'优先呈现胶原纤维形成 装配/交联及基质沉积相关基因上调，按效应 原始P FDR和患者一致性整理；耐药及抑制相关机制暂缓，不作为当前阶段门槛。','current_collagen_gene_evidence':str(O/'collagen_matrix_upregulated_primary.tsv'),'priority_updated_at':now})
 p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
assert (q.mean_effect_excluded>0).all() and (q.P_excluded<.05).all()
(O/'ACCEPTANCE.json').write_text(json.dumps({'passed':True,'selection_rows':len(q),'unique_genes':q.gene.nunique(),'all_sources_retained_in_complete_table':True,'raw_and_corrected_retained':True,'FDR_unchanged':True,'progress_entries_updated':3,'new_model_fit':False},indent=2))
print(json.dumps({'passed':True,'nominal_rows':len(q),'unique_genes':q.gene.nunique(),'by_subgroup':plan['screening_summary'],'F2_C3_examples':q[(q.subgroup_id=='F2_C3')&q.gene.isin(['COL3A1','COL5A1','PCOLCE','LOXL2'])][['gene','N_excluded','mean_effect_excluded','P_excluded','FDR_excluded']].to_dict('records')},ensure_ascii=False))
