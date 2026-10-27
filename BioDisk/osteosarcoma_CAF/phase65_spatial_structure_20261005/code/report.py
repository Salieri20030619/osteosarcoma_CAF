from pathlib import Path
import pandas as pd,json,hashlib,datetime,zipfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005'
assert json.loads((O/'QA/ACCEPTANCE.json').read_text())['passed']
st=pd.read_csv(O/'tables/spatial_patient_statistics.tsv',sep='\t');ds=pd.read_csv(O/'tables/discovery_structure_statistics.tsv',sep='\t');lo=pd.read_csv(O/'tables/discovery_patient_leave_one_out.tsv',sep='\t');dp=pd.read_csv(O/'tables/discovery_patient_effects.tsv',sep='\t');con=pd.read_csv(O/'tables/source_preference_statistics.tsv',sep='\t')
def spatial(scope,p,c,v='full',m='depth_xy'):
 r=st[(st.cohort_scope==scope)&(st.program==p)&(st.context==c)&(st.version==v)&(st.model==m)].iloc[0]
 return f"N={int(r.N)}，患者平均偏秩相关r={r.mean_effect:+.3f}，{int(r.positive_patients)}/{int(r.N)}例正方向，P={r.P:.8f}，FDR={r.FDR:.8f}"
def disc(state,program,version,scope='all',minimum=1,matrix='corrected'):
 return ds[(ds.dataset=='SCPCP000017')&(ds.subgroup_id==state)&(ds.program==program)&(ds.version==version)&(ds.scope==scope)&(ds.minimum_cells==minimum)&(ds.matrix==matrix)].iloc[0]
def desc(r):return f"N={int(r.N)}，平均效应={r.mean_effect:+.3f}，{int(r.positive_patients)}/{int(r.N)}例正方向，P={r.P:.6f}，FDR={r.FDR:.6f}"
report="""骨肉瘤CAF：胶原结构、患者稳健性与空间来源验证
2026-10-05｜phase65｜完整分析报告

一、关键结论与项目选择

本轮支持继续优先研究F2-C3相关纤维胶原生成，但不把它直接升级为已验证的耐药CAF。F2-C3去除COL1A1/COL1A2后仍有7例患者一致的剩余纤维胶原优势，逐一去掉任何患者后平均方向也保持为正。

Ⅳ型胶原需要单独处理：F2-C1的优势集中于COL4A1/COL4A2；组织空间中Ⅳ型胶原与血管内皮及血管壁程序的关联，比与间质程序的关联更强。这提示后续须优先排查血管／基底膜相关解释，不能把COL4A1高表达直接当成蛋壳状CAF证据。

原发灶纤维胶原同时关联间质与成骨程序，没有证据确认间质关联强于成骨关联。因此，F4及骨样基质来源必须继续保留在主线比较中，不能把所有胶原都归给CAF。

二、Methods

1. 数据与预先固定的比较
沿用8174个原始细胞、F1–F5及33个内部亚群；逐细胞原标签和各类置信度完整保留。采用phase63患者/部位配对框架，同时分析原始与校正计数、全8174内部参考与Fibro限制参考，以及每侧至少1或5个细胞的条件。
空间数据为GSE293065共26张切片、19个患者记录、53309个spot。原发灶分析含9个患者记录。MR4与PT4各有1个零库spot，保留条目，规范化分数设为缺失，不进入相关系数计算。
在查看新增结果前冻结3个胶原结构面板、4个不重叠来源背景面板及敏感性方案。来源面板是用于比较RNA背景的分析者定义，不是经验证的细胞分类器。

2. 面板
纤维胶原：COL1A1、COL1A2、COL3A1、COL5A1、COL5A2、COL5A3。
Ⅳ型胶原：COL4A1–COL4A6。
软骨相关基质：COL2A1、ACAN、COL9A1、COL9A2、COL9A3、COL11A1、COL11A2。
间质背景：DCN、LUM、PDGFRA、FBLN1、DPT、PI16。
血管壁背景：RGS5、CSPG4、NOTCH3、KCNJ8、ABCC9、MCAM。
内皮背景：PECAM1、CDH5、VWF、KDR、ESAM、EMCN。
成骨背景：RUNX2、SP7、ALPL、BGLAP、IBSP、IFITM5。
三项结构敏感性分别剔除COL1A1/COL1A2、COL4A1/COL4A2、COL2A1/ACAN；不按新结果调整基因。

3. 计分与患者统计
单细胞候选分析以同来源、同患者、同部位伪bulk的全库CP10k加0.1后取log2差，再对面板基因等权平均；部位等权后形成患者效应。
空间分数为面板基因log1p(全库CP10k)的均值，要求至少3个且覆盖≥80%基因。计算每切片原始秩相关、调整测序深度与x/y坐标的偏秩相关，以及进一步调整其余3个背景程序的敏感性模型。重复切片在患者/部位/事件内平均，事件在患者/部位内平均，再对患者的不同部位等权平均。所有P值以独立患者效应为单位，不使用spot级P。
对患者效应作双侧精确符号翻转检验；N<2或全零不计算P。原始P<0.05用于名义候选筛选。空间12项结构-背景关联在各版本/模型/患者范围内单独BH校正；9项预设关联差值构成独立家族。发现集在来源/计数版本/参考范围/细胞阈值/面板版本内对全部亚群结构检验统一BH。绝对效应<1e-12作为数值零处理。
逐一删除患者只评价平均效应方向及范围，不对多轮删除结果重新挑选最小P值。原发灶N=9及发现集小样本均定位为探索性。

三、Results

3.1 F2-C3纤维胶原特征不依赖COL1A1/COL1A2或单一患者

"""
for v in ['full','drop_dominant_pair']:report+=('完整六基因面板：' if v=='full' else '去除COL1A1/COL1A2：')+desc(disc('F2_C3','fibrillar_I_III_V',v))+'。\n'
rr=lo[(lo.dataset=='SCPCP000017')&(lo.subgroup_id=='F2_C3')&(lo.program=='fibrillar_I_III_V')&(lo.version=='drop_dominant_pair')&(lo.matrix=='corrected')&(lo.scope=='all')&(lo.minimum_cells==1)].iloc[0]
report+=f"去除主导基因对后，逐一去掉任一患者的平均效应范围为{rr.loo_min:+.3f}至{rr.loo_max:+.3f}，全部保持正方向。该结果为小样本名义显著证据，未达到发现集相应家族FDR<0.05。\n"
for scope,minimum in [('Fibro',1),('all',5)]:
 for v in ['full','drop_dominant_pair']:
  r=disc('F2_C3','fibrillar_I_III_V',v,scope,minimum)
  report+=f"敏感性（{scope}参考、每侧≥{minimum}细胞、{v}）："+desc(r)+'。\n'
report+="""
3.2 F2-C1的Ⅳ型胶原优势集中于COL4A1/COL4A2，功能轴仍有重叠

"""
for v in ['full','drop_dominant_pair']:report+=('完整Ⅳ型胶原面板：' if v=='full' else '去除COL4A1/COL4A2：')+desc(disc('F2_C1','collagen_IV',v))+'。\n'
report+='F2-C1的纤维胶原剩余基因面板（去除COL1A1/COL1A2）也呈正效应：'+desc(disc('F2_C1','fibrillar_I_III_V','drop_dominant_pair'))+'。\n'
report+="因此采用“COL4A1/COL4A2相关结构支路”和“纤维胶原生成支路”的可重叠描述，避免把F2-C1与F2-C3当成互斥功能类别。去除特定胶原链后不升高，表示链组成不同，不等于原始COL4A1/COL4A2信号无效。\n\n3.3 纤维胶原与间质背景相关，但原发灶仍存在成骨来源混合\n\n"
report+='全部19个患者记录，纤维胶原—间质：'+spatial('all_patients','fibrillar_I_III_V','fibro_context')+'。\n'
report+='剔除COL1A1/COL1A2后：'+spatial('all_patients','fibrillar_I_III_V','fibro_context','drop_dominant_pair')+'。\n'
report+='进一步调整其他背景程序后：'+spatial('all_patients','fibrillar_I_III_V','fibro_context','full','mutual_contexts')+'。\n'
report+='原发灶9例的探索性结果，纤维胶原—间质：'+spatial('primary_only','fibrillar_I_III_V','fibro_context')+'。\n'
report+='原发灶纤维胶原—成骨：'+spatial('primary_only','fibrillar_I_III_V','osteogenic_context')+'。\n'
r=con[(con.cohort_scope=='primary_only')&(con.version=='full')&(con.program=='fibrillar_I_III_V')&(con.model=='depth_xy')&(con.contrast=='fibro_context_minus_osteogenic_context')].iloc[0]
report+=f"原发灶“间质关联减成骨关联”的患者平均差为{r.mean_effect:+.3f}，P={r.P:.6f}，FDR={r.FDR:.6f}。当前数据不支持在原发灶将纤维胶原解释为间质专属来源。整体队列的间质偏向不能替代原发灶的这一判断。\n"
report+='原发灶去除COL1A1/COL1A2后，剩余纤维胶原—间质：'+spatial('primary_only','fibrillar_I_III_V','fibro_context','drop_dominant_pair')+'；属于名义显著探索性信号，未通过该敏感性家族FDR<0.05。\n'
report+="\n3.4 Ⅳ型胶原与血管背景的关联稳定，不能直接解释为蛋壳\n\n"
for c in ['endothelial_context','mural_context','fibro_context']:report+='全部患者，Ⅳ型胶原—'+c+'：'+spatial('all_patients','collagen_IV',c)+'。\n'
for c in ['endothelial_context','mural_context']:report+='原发灶探索性结果，Ⅳ型胶原—'+c+'：'+spatial('primary_only','collagen_IV',c)+'。\n'
report+='全部患者去除COL4A1/COL4A2后，Ⅳ型剩余基因—内皮：'+spatial('all_patients','collagen_IV','endothelial_context','drop_dominant_pair')+'。\n'
report+='对应原发灶结果：'+spatial('primary_only','collagen_IV','endothelial_context','drop_dominant_pair')+'。部分原发切片剩余程序为常量，相关系数不可评估，实际有效患者减为7例，未填造P值。\n'
report+="\n3.5 软骨相关基质更接近成骨背景，继续保留F4来源比较\n\n"
report+='全部患者，软骨基质—成骨：'+spatial('all_patients','cartilage_matrix','osteogenic_context')+'。\n'
report+='原发灶探索性结果：'+spatial('primary_only','cartilage_matrix','osteogenic_context')+'。\n'
report+="该结果支持后续组织学区分骨样/软骨基质与间质纤维胶原，不能据此将空间区域映射为某个具体F4亚群或判定良恶性。\n\n四、Discussion\n\n"
report+="本轮验证加强了胶原结构和RNA来源背景的区分，未新增耐药因果证据。phase64两套真实应答队列的阴性主结果继续保留；不能用本轮空间显著性抵消临床关联缺失。\n"
report+="背景面板与胶原面板无直接基因重叠，减少了共享基因造成的直接循环，但spot仍是混合组织，多个程序共存不能证明由同一种细胞产生，也不能定位到F2-C1或F2-C3。成骨程序可来自成骨样肿瘤细胞及其他成骨来源，不等于正常成骨细胞；背景标志基因也不是绝对特异。\n"
report+="原发灶仅9个患者记录，发现集主候选仅7例，须保留探索性定位；更严格细胞数门槛降低有效患者数。COL4A3–A6在部分原发灶缺乏可评估变化，链组成敏感性不应解释为否认Ⅳ型胶原本身存在。\n"
report+="现有26张图像均为lowres预览，没有可直接采用的高分辨率病理图或已接受独立肿瘤轮廓。RNA与H&E叠加只显示位置，不等于胶原蛋白沉积、交联强度、闭合包裹或药物排斥。本次从本机重试BS3原论文公开Code Ocean入口仍返回HTTP403，未获得其原始空间数据，没有绕过访问限制。\n\n五、下一步进入独立组织验证\n\n"
report+="主线优先F2-C3相关纤维胶原，组织学比较必须同时包含血管相关COL4A1/COL4A2与成骨/软骨基质。当前可执行计算已完成，继续确认物理屏障需补充：可判读的高分辨率H&E及独立病理轮廓、实际胶原纤维/蛋白沉积图层、以及后续受控药物进入与效应实验。\n"
report+="已准备26张匿名H&E预览和26份空GeoJSON轮廓模板，匿名图像与spot坐标配准已检查。病理复核包不含RNA分数，空模板均标明未标注，当前接受的独立肿瘤轮廓仍为0。低分辨率图不足以判断时必须标记不确定并获取原始高分辨率资料，不能强行补画。独立轮廓通过后再做边缘距离、周界覆盖、厚度和连续性；这些几何结果仍需胶原特异染色/纤维读数支撑。\n"
report+="干预阶段应以改善药物进入或肿瘤细胞处理反应为功能终点，同时设置基质和直接细胞作用的区分对照。现阶段不根据RNA富集直接决定临床用药或把候选酶认定为已验证靶点。蛋壳仍作为侧枝，不预设其存在。\n\n六、验收与交付\n\n"
report+="独立复核重新计算全部面板分数与切片相关系数、患者聚合、各检验家族FDR及逐一删除患者结果；重现phase63全部1032行完整结构效应；核对8174个细胞原标签和置信度；独立抽核精确符号翻转P值。全部通过。全26张切片图均交付，未只选择外观支持主线的切片。\n"
report+="关键文件：tables/spatial_patient_statistics.tsv、source_preference_statistics.tsv、discovery_structure_statistics.tsv、两类patient_leave_one_out.tsv；全部原始P/FDR、正负结果及实际N均保留。病理匿名编码对应表单独保存在resources/PATHOLOGY_PRIVATE_KEY.tsv，不在独立病理复核ZIP内。\n"
report+="程序背景参考（本轮面板为分析者预设，非原文已验证分类器）：\nhttps://www.nature.com/articles/s41467-020-17740-1\nhttps://pmc.ncbi.nlm.nih.gov/articles/PMC8436943/\n空间数据：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293065\n原论文：https://doi.org/10.1038/s41467-026-74603-x\nBS3资料入口：https://codeocean.com/capsule/9535428/tree/v1\n"
name='CAF胶原结构与空间来源验证报告_20261005.txt';(O/'delivery'/name).write_text(report,encoding='utf-8')
# Patient level visualization with paired full / reduced panels; original donor identifiers stay in statistical tables.
fig,axs=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
for ax,state,program in zip(axs,['F2_C3','F2_C1'],['fibrillar_I_III_V','collagen_IV']):
 q=dp[(dp.dataset=='SCPCP000017')&(dp.subgroup_id==state)&(dp.program==program)&(dp.matrix=='corrected')&(dp.scope=='all')&(dp.minimum_cells==1)]
 w=q.pivot(index='donor_record',columns='version',values='effect')
 for _,r in w.iterrows():ax.plot([0,1],[r['full'],r.drop_dominant_pair],'-o',alpha=.7,lw=1.1,ms=4)
 ax.axhline(0,color='grey',ls='--',lw=.8);ax.set_xticks([0,1],['Full panel','Remove dominant pair']);ax.set_ylabel('Patient paired mean log2(CP10k + 0.1) difference');ax.set_title(state+' | '+program+'\n7 patients; exploratory')
fig.savefig(O/'figures/discovery_dominant_gene_sensitivity.png',dpi=170);plt.close(fig)
dec=pd.DataFrame([
['F2_C3','纤维胶原机制候选优先','剔除COL1A1/2后7/7正方向，LOPO同向；小样本名义信号','未建立低应答/因果耐药关联','优先验证胶原实际沉积及其细胞来源'],
['F2_C1','COL4A1/2及收缩比较支路，可与纤维胶原重叠','COL4优势由A1/A2驱动；空间IV与血管背景相关','不能解释为已证实蛋壳或确定F2_C1空间位置','对照血管/基底膜与间质纤维胶原'],
['F4_osteogenic','主线来源鉴别对照','原发纤维胶原同时关联成骨；软骨基质与成骨关联','未判定空间信号良恶性或具体F4亚群','保留骨样/软骨基质病理分区'],
['eggshell','独立病理验证侧枝','26图已整理；0独立轮廓；均lowres','尚无包裹几何与药物屏障验证','补高分辨率切片与独立轮廓后才做几何']
],columns=['candidate','role','support','unresolved','next'])
dec.to_csv(O/'tables/candidate_decisions.tsv',sep='\t',index=False)
(O/'resources/BS3_ACCESS_RECHECK.json').write_text(json.dumps({'url':'https://codeocean.com/capsule/9535428/tree/v1','route':'local Windows curl','HTTP_status':403,'raw_BS3_obtained':False},indent=2))
print('REPORT_WRITTEN',name,'characters',len(report),flush=True)
