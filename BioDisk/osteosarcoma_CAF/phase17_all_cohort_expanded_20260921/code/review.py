from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase17_all_cohort_expanded_20260921'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
q=json.loads((O/'delivery/QA.json').read_text());s=read('sample_layer_counts.tsv');c=read('expanded_cluster_support.tsv');a=read('strict_expanded_direction_agreement.tsv');labels=read('expanded_cell_labels.tsv.gz')
assert len(s)==52 and s.QC_cells.sum()==299483 and s.strict.sum()==4099
assert len(labels)==q['expanded_cells'] and labels.cell_id.is_unique
assert labels['strict'].sum()==4099 and labels.supported.sum()==q['supported_cells']
supported=read('supported_cell_labels.tsv.gz');strict=read('strict_cell_labels.tsv.gz')
overlap=[];supported_summary=[]
for name,other in [('supported',supported),('strict',strict)]:
 joined=labels[['cell_id','cluster']].merge(other[['cell_id','cluster']],on='cell_id',validate='one_to_one',suffixes=('_expanded','_refit'))
 for (cl,refit),g in joined.groupby(['cluster_expanded','cluster_refit']):
  overlap.append({'layer':name,'expanded_cluster':cl,'refit_cluster':refit,'shared_cells':len(g),'fraction_of_expanded_cluster_remaining_in_layer':len(g)/sum(joined.cluster_expanded==cl)})
pd.DataFrame(overlap).to_csv(O/'tables/refit_cluster_correspondence.tsv',sep='\t',index=False)
for cl,g in supported.groupby('cluster'):
 counts=g.groupby('donor_record').size()
 supported_summary.append({'cluster':cl,'cells':len(g),'studies':g.dataset.nunique(),'donor_records_ge10':int((counts>=10).sum()),'strict_cells':int(g.strict.sum()),'myCAF_raw':int(g.myCAF_customer_positive.sum()),'mural_fraction':float(g.mural_flag.mean())})
pd.DataFrame(supported_summary).to_csv(O/'tables/supported_cluster_support.tsv',sep='\t',index=False)
programs=[z for z in supported if z.endswith('_score')];se=[]
for cl in sorted(supported.cluster.unique()):
 for donor,g in supported.groupby('donor_record'):
  for layer,gg in [('supported',g),('strict',g[g.strict])]:
   inside=gg[gg.cluster==cl];outside=gg[gg.cluster!=cl]
   if len(inside)<10 or len(outside)<10:continue
   for program in programs:se.append({'cluster':cl,'donor_record':donor,'layer':layer,'program':program,'delta_mean_logCP10k':float(inside[program].mean()-outside[program].mean())})
se=pd.DataFrame(se);se.to_csv(O/'tables/supported_frozen_cluster_effects.tsv',sep='\t',index=False)
sa=[]
for (cl,program),g in se.groupby(['cluster','program']):
 pv=g.pivot(index='donor_record',columns='layer',values='delta_mean_logCP10k')
 if 'strict' not in pv:continue
 both=pv[['supported','strict']].dropna()
 sa.append({'cluster':cl,'program':program,'paired_donor_records':len(both),'same_direction':int((np.sign(both.supported)==np.sign(both.strict)).sum()),'supported_median_effect':both.supported.median(),'strict_median_effect':both.strict.median()})
sa=pd.DataFrame(sa);sa.to_csv(O/'tables/supported_strict_direction_agreement.tsv',sep='\t',index=False)
obj=ad.read_h5ad(O/'objects/supported_integrated.h5ad',backed='r');xy=obj.obsm['X_umap'].copy();assert obj.obs_names.tolist()==supported.cell_id.tolist();obj.file.close()
fig,axs=plt.subplots(1,3,figsize=(17,5))
for ax,col in zip(axs,['cluster','dataset','myCAF_customer_positive']):
 for name in sorted(supported[col].astype(str).unique()):
  mask=supported[col].astype(str).eq(name).to_numpy();ax.scatter(xy[mask,0],xy[mask,1],s=1,alpha=.6,label=name,rasterized=True)
 ax.set_title(col);ax.set_xticks([]);ax.set_yticks([]);ax.legend(markerscale=5,fontsize=6,loc='best')
fig.tight_layout();fig.savefig(O/'figures/supported_joint_umap.png',dpi=180);plt.close(fig)
surf=read('surface_screen_all_QC.tsv.gz')
sd=surf.groupby(['donor_record','site','markers'])[['QC_cells','target_cells','gate_cells','target_in_gate','mural_in_gate','CNV_high_in_gate']].sum().reset_index()
sd['coverage']=sd.target_in_gate/sd.target_cells.replace(0,np.nan);sd['known_target_fraction']=sd.target_in_gate/sd.gate_cells.replace(0,np.nan)
sd['high_fraction']=sd.CNV_high_in_gate/sd.gate_cells.replace(0,np.nan);sd['mural_fraction']=sd.mural_in_gate/sd.gate_cells.replace(0,np.nan)
sd.to_csv(O/'tables/surface_donor_descriptions.tsv.gz',sep='\t',index=False)
rr=[]
for marker,g in sd[sd.target_cells>=10].groupby('markers'):
 rr.append({'markers':marker,'donor_site_records':len(g),'median_coverage':g.coverage.median(),'median_known_target_fraction':g.known_target_fraction.median(),'median_high_fraction':g.high_fraction.median(),'median_mural_fraction':g.mural_fraction.median(),'rank_score':(g.coverage*g.known_target_fraction.fillna(0)).median()})
rank=pd.DataFrame(rr).sort_values('rank_score',ascending=False);rank.to_csv(O/'tables/surface_exploratory_ranking.tsv',sep='\t',index=False)
cohort=s.groupby('dataset')[['strict','supported','expanded','supported_myCAF_raw','expanded_myCAF_raw']].sum().reset_index()
cohort.to_csv(O/'tables/dataset_layer_counts.tsv',sep='\t',index=False)
lines=['# 全队列拓展层探索结果','',
f'纳入全部52条样本记录的299,483个QC后细胞。拓展候选共{q["expanded_cells"]:,}个，来自{q["expanded_samples"]}条样本记录；其中至少一种CNV低信号、另一种无高信号且参考检查通过的支持层有{q["supported_cells"]:,}个。原严格层4,099个全部保留。',
f'拓展层中{q["expanded_raw_myCAF"]:,}个细胞符合客户6基因面板的原始表达myCAF条件；严格层按同一原始表达规则为{q["strict_raw_myCAF"]:,}个。此口径不要求原始和去污染后同时阳性，因此不能直接与此前1,926个双重表达条件候选比较增量。','',
'## 全队列联合分析已形成探索性表达群','',
f'使用完整来源矩阵筛选候选后，取{q["shared_genes"]:,}个共同实测基因作联合建模。Harmony校正研究/技术批次，固定Leiden分辨率0.5，得到{q["clusters"]}个表达群。另在支持层和严格层独立重拟合相同流程。客户面板分数使用各来源完整基因空间，避免共同基因交集造成的缺失误判。','',
'|群|细胞|研究数|≥10细胞供者记录|严格层细胞|原始面板myCAF|血管壁标志比例|CNV冲突比例|',
'|---|---:|---:|---:|---:|---:|---:|---:|']
for r in c.to_dict('records'):lines.append(f'|{r["cluster"]}|{r["cells"]}|{r["datasets"]}|{r["donor_records_ge10"]}|{r["strict_cells"]}|{r["raw_myCAF_cells"]}|{r["mural_fraction"]:.1%}|{r["CNV_conflict_fraction"]:.1%}|')
lines+=['','### 无高信号冲突的支持层独立聚类','', '|群|细胞|研究数|≥10细胞供者记录|严格层细胞|原始面板myCAF|血管壁标志比例|','|---|---:|---:|---:|---:|---:|---:|']
for r in supported_summary:lines.append(f'|{r["cluster"]}|{r["cells"]}|{r["studies"]}|{r["donor_records_ge10"]}|{r["strict_cells"]}|{r["myCAF_raw"]}|{r["mural_fraction"]:.1%}|')
lines+=['','支持层固定群标签的myCAF程序对照：','']
for r in sa[(sa.program=='myCAF_customer_score')&(sa.paired_donor_records>=2)].to_dict('records'):
 lines.append(f'- 支持层群{r["cluster"]}：{r["paired_donor_records"]}条可比较供者记录，{r["same_direction"]}条与严格细胞子集同方向；支持层/严格子集的群内减群外中位效应为{r["supported_median_effect"]:.3f}/{r["strict_median_effect"]:.3f}。')
lines+=['','## 八个来源均进入支持层','', '|来源|严格层|支持层|广义拓展层|支持层原始面板myCAF|','|---|---:|---:|---:|---:|']
for r in cohort.to_dict('records'):lines.append(f'|{r["dataset"]}|{r["strict"]}|{r["supported"]}|{r["expanded"]}|{int(r["supported_myCAF_raw"])}|')
lines+=['','## 严格层与拓展层的方向对照','',
'在拓展层冻结群标签后，分别在拓展、支持及严格细胞集合中计算同供者记录的群内减群外面板平均logCP10k差值。群内和群外各至少10个细胞用于本次探索性描述；10个细胞不是生物学样本量。完整效应与可比较记录数保存在within_donor_layer_effects.tsv和strict_expanded_direction_agreement.tsv。','',
'|群|程序|可比较供者记录|同方向记录|拓展层中位效应|严格层中位效应|','|---|---|---:|---:|---:|---:|']
aa=a[(a.comparison=='strict')&(a.evaluable_donor_records>=2)]
for r in aa.to_dict('records'):lines.append(f'|{r["cluster"]}|{r["program"]}|{r["evaluable_donor_records"]}|{r["same_direction_records"]}|{r["expanded_median_effect"]:.3f}|{r["comparison_median_effect"]:.3f}|')
if len(aa)==0:lines.append('|无达到至少2条供者记录的严格层对照||||||')
lines+=['','## 拓展纳入明显增加了可分析范围，也保留了身份风险信息','',
f'拓展层包含{q["CNV_conflict_cells"]:,}个CNV一低一高细胞、{q["CNV_no_low_cells"]:,}个没有低信号支持的未定细胞，以及{q["mural_flag_cells"]:,}个血管壁程序阳性细胞；这些标记可重叠。它们参与广义探索池，但不能因聚到同一群就视为已确认CAF。支持层对照用于检查探索结果是否依赖这些细胞。',
'原始Fib≥2，或Fib≥1且ECM≥3，作为拓展间质表达入口；最低基因数300。较强的免疫/内皮/上皮/肌源/破骨程序及已标记双细胞仍排除。成骨、血管壁、去污染不足及原簇风险改为记录属性；双方法均高、或仅高且没有低信号支持的细胞不进入本次联合池。原严格层不改写。','',
'## 表面组合已按拓展目标重新筛查','',
'16个既有候选及两两组合按原始RNA≥1筛查；排除门为PTPRC、PECAM1、EPCAM未检出。分母为有拓展候选样本的全部QC细胞，未包含零拓展候选样本，不能写为全52样本通用性能。目标定义为拓展层原始面板myCAF。表中门内已知目标比例属于RNA模拟结果。排序表同时保留血管壁和CNV高信号比例。','']
for r in rank.head(5).to_dict('records'):lines.append(f'- {r["markers"]}：{r["donor_site_records"]}条供者/部位记录；覆盖中位数{r["median_coverage"]:.1%}，门内已知目标比例中位数{r["median_known_target_fraction"]:.1%}。')
lines+=['','## 解释与下一步','',
'本轮实现了全队列宽入口联合探索及严格层对照，没有把有限的新增ScPCA子集代替全项目。各群先按表达程序与供者支持评估，随后优先复核有跨记录支持、严格层方向保留且CNV/血管壁风险较低的群。',
'当前群和myCAF均为操作性候选，未建立转移相关新亚型或最终流式panel。程序定义与表达验证使用部分相同基因，因此方向一致属于内部稳定性描述，并非独立验证。新增跨来源患者对应仍未完成核实，不将供者记录数直接作为独立患者数，不进行细胞级疾病差异检验。未将非肺样本统一解释为无转移A组；已有临床标签保留于逐细胞表，未进行未知临床状态的强制填补。',
'本轮未重跑原始QC/CNV、未修改旧候选集合或客户原始数据；完整参数、逐细胞分层、完整矩阵对象与图件均保存在本阶段目录。']
(O/'delivery/全队列拓展探索结果.md').write_text('\n'.join(lines)+'\n')
q['review_counts_pass']=True;q['surface_scope']='all QC cells only in samples contributing expanded candidates';(O/'delivery/QA.json').write_text(json.dumps(q,indent=2))
h={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['tables','figures','delivery','code'] for p in (O/folder).glob('*') if p.is_file() and p.name!='SHA256.json'}
(O/'delivery/SHA256.json').write_text(json.dumps(h,indent=2))
(O/'REVIEW_COMPLETE').write_text('DESCRIPTIVE_RESULTS_REVIEWED\n')
statefile=P/'PROJECT_STATE.json'
snapshot=O/'PROJECT_STATE_before_phase17.json'
if not snapshot.exists():snapshot.write_bytes(statefile.read_bytes())
state=json.loads(statefile.read_text());state.update(updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),current_stage=O.name,current_status='EXPANDED_EXPLORATION_COMPLETE',status='EXPANDED_EXPLORATION_COMPLETE',dynamic_status_file=str(O/'STATUS.json'),latest_report=str(O/'delivery/全队列拓展探索结果.md'),next_action='Review supported and broad expression states and per-donor direction agreement; exploratory candidates, no validated metastasis-specific subtype')
state['phase17_all_cohort_expanded_20260921']=q
statefile.write_text(json.dumps(state,ensure_ascii=False,indent=2))
print(json.dumps(q),flush=True)
