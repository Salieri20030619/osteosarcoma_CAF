from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012';B=P/'phase82_allCAF_crossdataset_markers_20261012'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
ret=read('sample_retention_waterfall.tsv');directions=read('all14_markers_all_sources_direction_evidence.tsv');summary=read('MARKER_DIRECTION_SUMMARY.tsv');rec=read('SCPCP000018_late_recovery_direction_sensitivity.tsv');refs=read('designated_reference_CNV_calibration_summary.tsv')
assert len(ret)==30 and ret.gate_mismatches.sum()==0 and ret.cell_ids_unique.all() and ret.final_ids_subset.all()
parts=['fibroblast_known_doublet','fibroblast_concordant_alteration','fibroblast_single_alteration','fibroblast_method_conflict','fibroblast_unassessed','stage27_CAF_candidates']
assert (ret[parts].sum(axis=1)==ret.provisional_fibroblast_cells).all()
assert (rec.original_direction==rec.recovery_direction).all()
frozen=pd.read_csv(B/'tables/all_candidate_source_layer_design_statistics.tsv',sep='\t');main=pd.read_csv(B/'tables/primary_same_layer_crossdataset_markers.tsv',sep='\t')
probe=directions[directions.effect.notna()].merge(frozen[(frozen.sensitivity=='all_cells')&(frozen.design=='independent')],on=['layer','dataset','symbol'],suffixes=('_new','_old'),validate='one_to_one')
assert np.allclose(probe.effect,probe.mean_lung_minus_bone)
assert np.allclose(probe.P_new,probe.P_old,equal_nan=True) and np.allclose(probe.FDR_new,probe.FDR_old,equal_nan=True)
ref17=refs[refs['sample'].eq('GSE152048__BC17')].iloc[0];assert ref17.designated_reference_cells==471 and ref17.infercnv_alteration_called==471
aggregate=ret.groupby('dataset')[['input_cells','provisional_fibroblast_cells']+parts+['final_candidate_cells','after_stage27_removed']].sum();save(aggregate.reset_index(),'retention_by_dataset.tsv')
site=pd.read_csv(O/'tables/official_SCPCP000018_metadata.tsv',sep='\t',keep_default_na=False)
site=site[site.scpca_sample_id.isin(['SCPCS000527','SCPCS000528','SCPCS000530'])][['scpca_sample_id','participant_id','diagnosis','subdiagnosis','disease_timing','tissue_location','tissue_ontology_term_id']]
site['decision']='retain_Unknown_or_mixed';site['reason']='Region label does not identify sampled bone versus soft-tissue compartment; no precise primary evidence obtained';save(site,'ANATOMICAL_SITE_REVIEW_DECISIONS.tsv')
sources=['SCPCP000017','SCPCP000023','SCPCP000018','GSE152048','GSE318840']
summary['any_opposite']=summary.opposite_sources.gt(0)
summary=summary.sort_values(['any_opposite','same_sources','symbol'],ascending=[True,False,True]).reset_index(drop=True)
summary['direction_readout']=np.where(summary.opposite_sources.eq(0),'no_opposite_among_observed_effects','discordant_sources_present')
save(summary,'DIRECTION_REVIEW_PRIORITY.tsv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
codes={'opposite':-1,'zero':0,'same':1,'not_evaluable':2};grid=[];text=[]
for r in summary.itertuples():
    vals=[];labels=[]
    for ds in sources:
        q=directions[(directions.layer=='expanded')&(directions.symbol==r.symbol)&(directions.dataset==ds)].iloc[0]
        vals.append(codes[q.direction]);labels.append(f'{q.effect:+.2f}'+('*' if q.P<.05 else '') if pd.notna(q.effect) else 'NA')
    grid.append(vals);text.append(labels)
fig,ax=plt.subplots(figsize=(12,8.5));cmap=ListedColormap(['#b85050','#ececec','#3e8c82','#d9e1ee']);norm=BoundaryNorm([-1.5,-.5,.5,1.5,2.5],4)
ax.imshow(grid,cmap=cmap,norm=norm,aspect='auto')
for i,row in enumerate(text):
    for j,val in enumerate(row):ax.text(j,i,val,ha='center',va='center',color='white' if grid[i][j] in [-1,1] else '#222',fontsize=10)
ax.set_yticks(range(len(summary)),summary.symbol);ax.set_xticks(range(len(sources)),['S17\nBone12 / Lung3','S23\nBone4 / Lung3','S18\nBone1 / Lung4','GSE152048\nBone2 / Lung1','GSE318840\nBone1 / Lung2']);ax.set_title('All14 candidate markers: agreement across five source datasets',fontsize=15,pad=15)
ax.legend(handles=[Patch(color='#3e8c82',label='Same direction as original pair'),Patch(color='#b85050',label='Opposite direction'),Patch(color='#ececec',label='Zero effect'),Patch(color='#d9e1ee',label='Not evaluable / unmapped')],loc='upper center',bbox_to_anchor=(.5,-.10),ncol=2,frameon=False)
fig.text(.015,.02,'Values: lung minus bone mean log2(CP10k + 0.1). * original raw P < 0.05. S18/GEO directions have insufficient biological replication.\nEBF1 reference direction is bone-higher (negative effect); all other reference directions are lung-higher. Zero or NA is not positive support.',fontsize=9)
fig.tight_layout(rect=[0,.13,1,.98]);fig.savefig(O/'figures/01_five_source_marker_directions.png',dpi=180);plt.close(fig)

fig,ax=plt.subplots(figsize=(12,8.5));y=np.arange(len(summary));rr=rec.set_index('symbol').loc[summary.symbol]
ax.scatter(rr.original_effect,y-.12,color='#326b9b',s=38,label='Frozen candidate pool (1 bone /4 lung)');ax.scatter(rr.recovery_effect,y+.12,color='#ce6840',s=38,label='+384 provisional recovered cells, same donors')
for i,r in enumerate(rr.itertuples()):ax.plot([r.original_effect,r.recovery_effect],[i-.12,i+.12],color='#777',lw=.7)
ax.axvline(0,color='#555',lw=.8);ax.set_yticks(y,summary.symbol);ax.invert_yaxis();ax.set_xlabel('SCPCP000018 mean lung minus bone log2(CP10k + 0.1)');ax.set_title('Late CNV recovery does not change marker directions',fontsize=15,pad=15);ax.legend(frameon=False,loc='lower right')
fig.text(.015,.015,'Recovered384cells satisfy the original provisional entry rule; final cluster/lineage review remains unresolved.\nThey belong to a lung donor already represented by447cells. No additional independent patient and no new inferential P.',fontsize=9)
fig.tight_layout(rect=[0,.07,1,.97]);fig.savefig(O/'figures/02_late_recovery_direction_sensitivity.png',dpi=180);plt.close(fig)

protocol=dict(user_request='检查上游问题；如果解决不了，在方向上看是否一致',new_inferential_P_values=0,site_reassignment='None without tissue-compartment evidence',formal_candidate_pool='unchanged; source data and earlier outputs preserved',late_recovery='Separate versioned384-cell provisional recovery; descriptive patient-level sensitivity only',direction='Same/opposite/zero/not evaluable retained; original P/FDR unchanged; one-sided cohorts never forced into two-group contrast',sources='Public ScPCA metadata, source processed objects, original scripts and CNV outputs, published paper and attempted supplementary retrieval')
(O/'resources/PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
source_note='''Public source verification
ScPCA project SCPCP000018 https://scpca.alexslemonade.org/projects/SCPCP000018
Associated primary article: Al Shihabi et al., The landscape of drug sensitivity and resistance in sarcoma, Cell Stem Cell (2024), DOI10.1016/j.stem.2024.08.010
https://pmc.ncbi.nlm.nih.gov/articles/PMC12318355/
The accessible article describes SARC0105 as primary osteosarcoma. This does not independently identify the sampled bone compartment within Arm.
Original ScPCA single-cell_metadata.tsv and processed objects were inspected; Arm/Leg/Thigh and tissue ontology region labels agree.
Supplement PDF retrieval did not succeed: NCBI returned an HTML preparation page; EuropePMC API returned500 on server and local TLS failed. Do not claim the supplementary clinical table was read.
GSE152048 original study https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/ and GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152048
'''
(O/'resources/PUBLIC_SOURCE_REVIEW.txt').write_text(source_note,encoding='utf-8')
report='''CAF上游保留问题核查与跨来源marker方向复核
2026年10月12日

主要结论
跨来源方向复核改变了候选优先级。NEBL在4个来源中同向、GSE152048为零差，当前没有观察到跨来源反向；GRID2、RGS2-AS1、AL096854.1、PEF1-AS1在3个来源同向，其余为零差或不可评估。RTKN、TAF10、TOR1AIP1在原2个来源名义显著，但另外3个可看方向的来源均反向，不宜仅凭双数据集P<0.05将其作为广泛适用marker优先项。
上游核查未发现细胞编号或既定入口规则的实现错配。确有384个后续CNV成功后可恢复的暂定候选未进入冻结池；单独加入恢复敏感性分析后，来源18的14个marker方向全部不变，也没有增加独立患者数。Arm等部位仍缺乏足够组织定位证据，暂未强行改组。

数据与方法
核查GSE152048、SCPCP000018和GSE318840共30份样本，逐细胞重建成纤维谱系、CNV与既有双细胞入口，并核对最终候选ID。检查提供方临床TSV、原processed对象、CNV原始调用及后补计算的退出码/完成QA。方向复核沿用Phase81患者级骨肺效应和Phase82规范特征映射，支持层与扩展层分开；主表呈现原14个扩展层候选在全部9个来源的状态，其中5个来源有两侧方向可看。
保留原始P、FDR、有效N、零效应和缺测原因，不把零或NA当同向支持，不把无生物学重复的来源补造P。恢复敏感性将同一患者两份肺样本的原始计数及完整文库相加后归一化，再按患者等权计算方向；不改变原细胞标签或正式候选池。

上游筛选减少的具体位置
'''
report+=aggregate.to_string()+'\n'
report+='''
GSE152048的BC17原有4032个输入细胞，操作性成纤维样表达3194个。其中3185个为两种CNV方法均有异常、8个方法冲突、1个已标记双细胞，因此原入口候选为0。不能把“未进入当前候选池”写成“该肺病灶生物学上没有CAF”。GSE152048其他样本也主要受CNV异常/冲突及后续竞争谱系复核影响，完整样本瀑布表保留。
入口规则重建与保存标签逐细胞一致，全部最终细胞ID都是原候选子集，30份样本未发现ID重复或静默丢失。但被指定的正常参考也常被inferCNV叫到非中性区段：BC17为471/471。来源18 Femur/Calcaneus样本同样存在参考异常调用。该现象提示参考注释或CNV校准风险，不能把所有排除细胞直接判为恶性，也不能未经验证就自动补为CAF。
既有Phase34在另两份样本进行免疫/内皮及淋巴参考重算，仍出现大量参考非中性区段；这不能直接修复BC17。当前没有可信的单行解析错误或独立正常基准，故未通过随意改变CNV门槛恢复样本。

SCPCP000018解剖部位复核
'''
report+=site.to_string(index=False)+'\n'
report+='''
SARC0105/SCPCS000528在提供方原始metadata和processed对象均为Arm，诊断Osteosarcoma，disease_timing=Primary；关联原始论文也描述其为原发骨肉瘤。但“原发”不能替代准确取材组织定位，当前未获得可区分骨组织与周围软组织的直接记录。SARC0103为Leg、SARC0129为Thigh（后者标为Metastasis）；同样保留区域字段，不推断骨取材。原论文补充文件未成功取得，未声称已完成补充临床表核查。
因此保留原Unknown_or_mixed标签。SCPCP000018明确骨部位候选仍仅1名患者；扩大范围需要原始病理或具体组织定位，而不是根据希望获得显著性改组。

后补CNV证据的恢复版本
SCPCS000525的后续CopyKAT重跑于2026-09-27完成，native_exit_code=0，41706个细胞验证通过。将这一更新CNV证据与原谱系/双细胞记录对应，384个细胞满足原暂定入口规则。已输出版本化的逐细胞恢复表及证据来源。
该样本和已纳入447个细胞的SCPCS000526属于同一肺侧患者SARC0088，补入敏感性后该患者合计831个候选记录，来源18仍是骨1/肺4。384个细胞没有经历正式最终聚类/谱系验收，故作为暂定恢复敏感性，不冒充已确认CAF或新增独立患者。14个marker的方向均未改变。

逐marker的跨来源方向
'''
report+=summary[['symbol','expected_direction','same_sources','opposite_sources','zero_sources','same_source_list','opposite_source_list','zero_source_list']].to_string(index=False)+'\n'
report+='''
方向较一致的候选
NEBL：来源17、23、18和GSE318840均肺侧较高，GSE152048为零差。GRID2：三个ScPCA来源肺侧较高，两个GEO来源为零差。RGS2-AS1及PEF1-AS1：三个ScPCA来源肺侧较高，GSE152048为零差，GSE318840特征不可评估。AL096854.1：三个ScPCA来源同向，两个GEO来源没有对应可评估特征。
这些额外来源只提供方向支持，仍只有来源17和23具备同层名义显著复现。方向排序不新增筛选硬门槛，14个原候选均完整保留。PEF1-AS1在来源18的效应很小（原+0.0358、恢复后+0.0296），应按实际效应大小理解。

关键反向结果
NDUFS7、PNKD：4个来源同向，但来源18反向。EBF1：4个来源骨侧较高，来源18肺侧较高，且来源17两名配对患者均与原独立组方向相反。BEST1、ZNF526：3个同向、2个反向。CALHM5：3个同向、GSE152048反向、GSE318840为零。RTKN、TAF10、TOR1AIP1：2个同向、另外3个来源反向。
上一轮提到RTKN/BEST1/RGS2-AS1在同来源支持层还有名义信号；该证据不等于跨来源广泛一致。本轮保留这些原始事实，同时将外部方向矛盾纳入候选解释。

讨论
扩展候选池包含身份未完全解决和成骨样表达细胞，上游CNV对正常参考的表现也有风险。NEBL等跨来源方向较稳的候选仍不能直接称为已验证CAF特异marker：NEBL在来源23 Fibro限定中无表达差，在来源17配对患者也未一致同向。来源18、GSE152048及GSE318840分别只有骨1/肺4、骨2/肺1、骨1/肺2名候选患者，方向结果不得升格为新的显著复现。原两数据集的14个候选仍均未通过相应FDR阈值。

文件
sample_retention_waterfall.tsv为30份样本逐步保留；ANATOMICAL_SITE_REVIEW_DECISIONS.tsv为部位审查；SCPCS000525_provisional_recovered_cells.tsv.gz为384个恢复记录；SCPCP000018_late_recovery_direction_sensitivity.tsv为恢复前后效应；all14_markers_all_sources_direction_evidence.tsv为全来源完整效应/P/FDR/N；DIRECTION_REVIEW_PRIORITY.tsv为方向摘要。图1展示五来源全部14marker，图2展示后补CNV恢复敏感性。既有正式结果未覆盖。

公开资料
https://scpca.alexslemonade.org/projects/SCPCP000018
https://pmc.ncbi.nlm.nih.gov/articles/PMC12318355/
https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/
'''
(O/'delivery/CAF上游核查与跨来源方向复核_20261012.txt').write_text(report,encoding='utf-8')
qa=dict(status='numeric_checks_passed_visual_review_pending',retention_samples=30,gate_mismatch_cells=0,retention_partition_verified=True,final_cell_ids_verified=True,source_effect_P_FDR_records_verified=len(probe),source_effect_P_FDR_unchanged=True,late_recovery_cells=384,late_recovery_changes_no_marker_sign=True,independent_patients_added=0,site_labels_changed=False,new_inferential_P_values=0,visual_review_pending=True,figures=2)
(O/'QA/NUMERIC_ACCEPTANCE.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps(qa));print('AGGREGATE\n'+aggregate.to_string())
