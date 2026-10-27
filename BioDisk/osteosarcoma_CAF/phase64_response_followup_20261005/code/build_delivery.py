from pathlib import Path
import pandas as pd,json,datetime,hashlib,shutil,zipfile,sys,importlib.metadata
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
assert json.loads((O/'QA/ACCEPTANCE.json').read_text())['passed']
assert json.loads((O/'FINAL_RUN_STATUS.json').read_text())['success']
st=pd.read_csv(O/'tables/program_response_statistics.tsv',sep='\t');ge=pd.read_csv(O/'tables/fixed_program_gene_response_statistics.tsv',sep='\t');cv=pd.read_csv(O/'tables/program_coverage.tsv',sep='\t')
def row(c,p,m='unadjusted_Welch'):return st[(st.cohort==c)&(st.program==p)&(st.model==m)].iloc[0]
def line(c,p,m='unadjusted_Welch'):
 r=row(c,p,m);return f"{c}：差值 {r.effect:+.3f}，95% CI [{r.ci_low:+.3f}, {r.ci_high:+.3f}]，P={r.p:.6f}"
report="""骨肉瘤 CAF：真实化疗应答验证与下一步机制研究
2026-10-05｜phase64｜内部完整分析报告

一、结论与项目决策

两套治疗前活检队列均未显示整体Ⅰ/Ⅲ/Ⅴ型纤维胶原程序与较差病理化疗应答显著相关。主分析共纳入83例有明确应答标签的患者，按队列分别计算，未合并表达矩阵。当前数据不支持把该程序直接作为通用耐药标志，也不支持将F2-C3命名为已经确认的“耐药CAF”。

上一阶段确立的F2-C3产纤维胶原候选特征仍是机制研究的入口。现阶段研究重点应转向“何种细胞在何处沉积何种胶原，以及改变这些基质是否改善药物作用”，而非继续依赖更换CAF分类名称。F2-C1保留为收缩／Ⅳ型胶原比较支路；F4成骨相关状态用于区分成骨基质与肿瘤周围间质沉积。

空间侧枝核查发现：GSE293065的公开作者代码仍以表达成分定义癌区，没有识别到可直接采用的独立病理肿瘤轮廓。BS3提供单病例的COL4A1高表达坏死前沿文献线索，尚不能确立普遍存在的蛋壳结构。

二、Methods

1. 数据与患者单位
GSE39055为37例独立诊断活检mRNA芯片样本，GPL14951平台；按坏死率≥90%分为高应答，<90%分为低应答。23例低应答、13例高应答纳入主分析；1例仅记载“>80% necrosis”，保留为分组不确定病例，另分别归入两组进行敏感性分析。
GSE21257为53例治疗前活检，GPL10295平台；Huvos I/II为低应答、III/IV为高应答。29例低应答、18例高应答纳入主分析；6例Unknown保留记录。患者为统计单位，未用细胞或spot数替代样本量。

2. 基因程序与探针注释
在查看关联结果前冻结主要终点：COL1A1、COL1A2、COL3A1、COL5A1、COL5A2、COL5A3六基因程序。其余13个程序沿用phase63，三项结构面板在phase63属事后分解，在本轮作为固定外部假设使用。
表达使用作者提交的lumi VST与标准化结果，不再取对数。采用官方GPL完整注释，优先依据HGNC Entrez与去版本号RefSeq识别基因；标识冲突者排除，其他探针回退至确定的当前／历史符号映射。同一基因的多个确定探针取中位数。对各队列全部提交样本计算基因z分数，再等权平均成程序分数；覆盖率≥80%且至少3个基因的程序进入分析。
两队列14个程序均满足覆盖要求；主要六基因、Ⅳ型胶原六基因均完整覆盖。GPL10295缺少可确定映射的PXDN，涉及该基因的程序按实际覆盖计算，未补造表达值。旧符号PRG2对应的NM_024888实际为PLPPR3，已按转录本标识处理。完整映射、缺失基因和修订记录归档。

3. 统计分析
主要效应为“低应答减高应答”的程序均值差，报告Welch双侧P值、95%置信区间和Hedges g；主要程序另作99,999次标签置换。调整模型为程序分数~应答组+年龄+性别，GSE21257另加诊断时转移状态，采用HC3稳健协方差及t推断。未来转移事件未用作调整协变量。
每队列单一主要终点不要求FDR；13个次要程序在各队列、各模型内单独BH校正，同时提供全部14个程序的补充BH结果。固定程序基因分别构成各队列基因检验家族。原始P<0.05为名义筛选标准，FDR不作为候选进入的硬门槛。标签未知者不插补，不将生存或推算IC50替代病理应答。

三、Results

3.1 整体纤维胶原程序未获得跨队列低应答关联支持

"""
for c in ['GSE39055','GSE21257']:
 r=row(c,'fibrillar_I_III_V');a=row(c,'fibrillar_I_III_V','adjusted_HC3');pm=row(c,'fibrillar_I_III_V','permutation_99999')
 report+=line(c,'fibrillar_I_III_V')+f"，Hedges g={r.hedges_g:+.3f}，补充全部14程序FDR={r.fdr_all14:.6f}。\n"
 report+=f"调整后差值 {a.effect:+.3f}，95% CI [{a.ci_low:+.3f}, {a.ci_high:+.3f}]，P={a.p:.6f}；置换P={pm.p:.6f}。\n"
report+="""
GSE39055中将唯一模糊标签分别归入高、低应答组，差值分别为+0.366和+0.090，P分别为0.198和0.721。结果未出现名义显著关联。两套队列效应方向不一致，均有置信区间跨零。

3.2 胶原交联及Ⅳ型胶原程序未显示低应答组升高的统计证据

"""
for p in ['collagen_crosslinking','collagen_IV']:
 report+=p+'\n'
 for c in ['GSE39055','GSE21257']:
  r=row(c,p);report+=line(c,p)+f"，次要程序家族FDR={r.fdr_secondary13:.6f}。\n"
report+="""
两队列全部14个程序在未调整和调整模型中均未达到原始P<0.05。最接近阈值的未调整结果为GSE21257矿化程序：差值−0.089，P=0.051777，FDR=0.673098，不列为名义显著候选。所有方向与统计值均保留在完整表中。

3.3 固定候选基因未形成跨队列重复的应答信号

GSE39055评估654个固定程序基因，其中13个原始P<0.05；GSE21257评估649个，其中24个原始P<0.05。两队列没有同时达到原始P<0.05的共同基因，两个基因家族均无FDR<0.05结果。单队列名义信号完整保留，不因FDR未通过而删除；它们属于单队列探索性线索。
"""
for gene in ['PCOLCE','LOXL2','LOX','PXDN']:
 report+=gene+'：'
 parts=[]
 for c in ['GSE39055','GSE21257']:
  d=ge[(ge.cohort==c)&(ge.gene==gene)]
  if len(d):
   r=d.iloc[0];parts.append(f"{c}效应{r.effect:+.3f} VST，P={r.p:.6f}，FDR={r.fdr_gene_family:.6f}")
  else:parts.append(f'{c}无确定映射，未检验')
 report+='；'.join(parts)+'。\n'
report+="""
上述主要靶点未获得本轮病理低应答关联支持。RNA定位或内部亚群差异不能因此直接升级为患者耐药靶点证据。

3.4 蛋壳状CAF仍是空间验证支路

已检查GSE293065原临床补充表、Source Data工作簿及作者Supplementary_Code.Rmd。代码以Cancer_sweep>0.4选取癌相关spot；该值来自表达IC成分。Seurat的boundaries/centroids字段记录坐标结构，不能视为病理画出的肿瘤边界。26张切片、19个患者记录中，本轮接受的独立病理轮廓为0；未新增环状包裹比例、厚度或药物排斥的验证性统计。
Zheng等2024年论文报告BS3在3次化疗后坏死率<50%，COL4A1高表达cluster6位于坏死前沿。这是Ⅳ型胶原空间支路的单病例文献线索，本轮未获得其原始空间数据并重分析。修正后的公开入口为Code Ocean capsule9535428，本次访问返回HTTP403。Zenodo14594702是另一组儿童骨肉瘤肺转移空间数据，不将其误当BS3，也未下载其多GB对象。

四、Discussion

本轮否定的是“当前bulk队列已经支持通用高胶原低应答关联”的说法，不能据此否定既有F2-C3内部产胶原特征。bulk组织包含肿瘤、成骨、血管及间质来源，无法直接定位F2-C3或识别边缘胶原结构。RNA分数也不等于实际胶原沉积量、成熟交联或药物通透性。
两队列使用不同芯片与处理流程，效应按队列解释；部分基因覆盖不同。GSE39055为FFPE资料且样本量36，GSE21257按未来转移结局选择入组，并非无选择的临床连续队列。现有区间仍允许一定方向与大小的效应，阴性结果不是已证明生物学效应为零。此次分析没有建立预测模型、治疗因果关系或临床可用的干预方案。
旧符号映射修正在最终统计前完成；中间单队列结果被正式结果替代。代码、方法冻结文件、注释修订、失败下载记录和最终验收均归档。

五、按最终目标安排下一步

第一优先：定位实际沉积的胶原。
对既有组织材料建立独立H&E肿瘤、坏死、骨样基质与血管轮廓，标注时不参考本轮胶原RNA分数。注册到空间坐标后，区分Ⅰ/Ⅲ/Ⅴ型纤维胶原、Ⅳ型血管／基底膜胶原和骨样基质。F2-C3作为纤维胶原来源候选，F2-C1作为收缩／Ⅳ型胶原比较，F4相关状态作为成骨来源比较；保留原F1–F5、33内部亚群及其置信度。
建议优先获得胶原蛋白／纤维组织学读数，并与原有RNA候选对应。仅H&E与RNA共定位不足以量化成熟交联。

第二优先：确认影响药物作用的基质性质。
若独立组织定位显示肿瘤边缘纤维胶原沉积，再评估纤维排列、密度、成熟交联与基质力学；并在同一实验中测量药物进入肿瘤区域的分布及细胞反应。主要功能终点应直接对应药物作用改善，不以CAF数量下降或胶原RNA下降代替。
功能验证采用肿瘤细胞与候选间质细胞／其沉积基质的对照体系，设置化疗、基质相关干预及联合处理，记录独立生物学重复；进行必要的基质或靶点恢复对照以区分基质效应和直接细胞毒作用。PCOLCE、LOXL2、PXDN目前为候选，具体先后顺序取决于蛋白定位与酶／交联读数，不能仅按本轮阴性的bulk关联硬排“耐药靶点”。

侧枝：蛋壳几何。
仅在独立肿瘤轮廓及胶原沉积图层具备后，计算周界覆盖率、内外侧梯度、连续性和厚度；以患者为推断单位，单病例只作描述。COL4A1坏死前沿、血管基底膜与Ⅰ/Ⅲ/Ⅴ胶原环应分别检验。若未形成连续边缘结构，保留局灶基质沉积解释，不强行命名蛋壳。

决策门槛：
（1）内部亚群RNA特征支持——可作来源候选；
（2）独立组织定位与物理基质读数支持——可作局部屏障机制候选；
（3）受控干预改善药物进入或抗肿瘤反应，并有机制对照——方可提升为功能性耐药机制证据。
本轮完成真实临床应答关联及空间资料核查；病理轮廓、实际胶原沉积／交联读数和干预实验仍未完成。

六、给客户的简短回复建议

“我们已把分析重点转向胶原基质与化疗反应，并完成两套治疗前活检队列的真实病理应答核查。目前整体胶原表达量没有显示稳定的低应答关联，因此不能直接把某个CAF亚群命名为耐药CAF。接下来保留已有内部亚群作为细胞来源线索，重点验证胶原沉积的位置、结构及改变基质后药物作用是否改善。蛋壳状结构作为空间支路继续检查，但不预设它一定存在。”

七、交付与来源

主统计：tables/program_response_statistics.tsv
患者标签与分数：tables/patient_metadata.tsv、patient_program_scores.tsv
完整基因与跨队列证据：tables/fixed_program_gene_response_statistics.tsv、cross_cohort_gene_evidence.tsv
覆盖及探针映射：tables/program_coverage.tsv、*_probe_mapping.tsv
主要结果图：figures/fibrillar_response.png / .pdf
冻结方案与修订：resources/CLINICAL_PROTOCOL_FROZEN.json、ANNOTATION_QA_AMENDMENT.json
空间核查：resources/SPATIAL_FOLLOWUP_CONCLUSIONS.json、SPATIAL_CODE_AUDIT_EXCERPTS.json
正式运行与验收：FINAL_RUN_STATUS.json、QA/ACCEPTANCE.json

GEO GSE39055：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39055
GEO GSE21257：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21257
原始mRNA研究：https://doi.org/10.1186/gm406
BS3研究：https://doi.org/10.1186/s13045-024-01598-7
GSE293065论文：https://doi.org/10.1038/s41467-026-74603-x
作者空间代码：https://zenodo.org/records/19369637
"""
name='CAF真实化疗应答验证与下一步计划_20261005.txt'
(O/'delivery'/name).write_text(report,encoding='utf-8')
# Compact decision table, separate the established cellular phenotype from unproven resistance.
decisions=pd.DataFrame([
['F2_C3','纤维胶原来源候选','phase63患者配对探索性支持','两队列整体程序未支持低应答关联','先做独立组织定位与沉积读数','机制候选；不是已验证耐药CAF'],
['F2_C1','收缩/Ⅳ型胶原比较','phase63患者配对探索性支持','Ⅳ型胶原程序未支持低应答关联','区分血管基底膜、坏死前沿及间质边缘','比较支路'],
['F4相关成骨状态','成骨/骨样基质比较','保留原内部亚群与置信度','未获临床耐药升级依据','区分骨样基质与间质纤维胶原','来源鉴别支路'],
['eggshell_CAF','环状包裹几何假设','当前0独立病理轮廓；BS3仅文献单病例','未建立药物排斥/因果耐药证据','待独立轮廓及胶原蛋白图层后检验','侧枝；不预设成立']
],columns=['candidate','role','cellular_or_spatial_evidence','response_evidence','next_action','confidence_statement'])
decisions.to_csv(O/'tables/updated_candidate_decisions.tsv',sep='\t',index=False)
for f in ['program_response_statistics.tsv','updated_candidate_decisions.tsv']:
 shutil.copy2(O/'tables'/f,O/'delivery'/f)
(O/'resources/SOFTWARE_VERSIONS.json').write_text(json.dumps({p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib']},indent=2))
(O/'QA/FIGURE_VISUAL_QA.json').write_text(json.dumps({'passed':True,'file':'figures/fibrillar_response.png','verified':'Patient points, group Ns, separate cohort axes, labels and P values legible; no clipping.'},indent=2))
# Hash exactly the retained source and output files; incomplete downloads/interim files stay remote only.
manifest=[]
for folder in ['code','sources','resources','tables','objects','figures','QA','logs']:
 for p in sorted((O/folder).rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and 'incomplete_attempt' not in p.name and not p.name.startswith('probe_'):
   manifest.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for p in [O/'FINAL_RUN_STATUS.json',O/'COMPUTATION_COMPLETE.json',O/'delivery'/name]:
 manifest.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'resources/DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
zip_path=O/'delivery/CAF_response_validation_20261005.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for rec in manifest:z.write(O/rec['path'],arcname='phase64_response_followup_20261005/'+rec['path'])
 z.write(O/'resources/DELIVERY_MANIFEST.json',arcname='phase64_response_followup_20261005/resources/DELIVERY_MANIFEST.json')
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
receipt={'zip':str(zip_path),'bytes':zip_path.stat().st_size,'sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest(),'files':len(manifest)+1,'report':str(O/'delivery'/name),'local_verified':False}
(O/'delivery/DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
print(json.dumps(receipt,ensure_ascii=False),flush=True)
