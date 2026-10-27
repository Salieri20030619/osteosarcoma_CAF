from pathlib import Path
import hashlib
import json
import shutil
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT = Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
OUT = PROJECT/'phase80_scheme1_anatomical_20261012'
SRC = PROJECT/'phase74_cohort_only_rescreen_20261009'
FEAS = PROJECT/'task_anatomical_cohort_feasibility_20261012_1115'
for name in ['tables','figures','delivery','QA','resources']:
    (OUT/name).mkdir(parents=True, exist_ok=True)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return pd.read_csv(p, sep='\t')
def save(df,name): df.to_csv(OUT/'tables'/name,sep='\t',index=False)
def dump(obj,p): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
now = datetime.now(timezone.utc).isoformat()
protocol = dict(confirmed_at=now,user_instruction='按照方案一进行',
    scope='F1/F3/F4 and existing internal states, supported/expanded separately',
    primary_groups=['All explicitly bone-site lesions','All explicitly lung-site lesions'],
    metastatic_history='Not a grouping or eligibility requirement; preserved annotation',
    other_or_mixed_sites='Retain separately; do not force into Bone or Lung',
    unit='Source-specific independent biological patients; patient-site aggregation',
    repeated_sampling='Same patient with both sites: paired comparison, removed from corresponding independent contrast',
    statistics='Reuse accepted Phase73 strict_v2 and Phase74 full-family P/FDR unchanged',
    gate='Valid target raw P<0.05 nominates exploratory results; no FDR gate for inclusion',
    no_new_cell_cutoffs=True,original_labels_preserved=True,new_inferential_tests=0,
    endpoints=['Captured proportions','Within-class/state gene expression','Patient-level functional programs'],
    enrichment='ORA modeled gene-list enrichment kept distinct from patient functional tests',
    flow_marker_stage='deferred; no marker feasibility or antibody filter',
    third_cohort_stage='not active under approved scheme1')
dump(protocol, OUT/'resources/PROTOCOL.json')
inputs = list((SRC/'tables').glob('*')) + [FEAS/'testable_class_state_bone_lung.tsv',FEAS/'all_class_state_bone_lung_feasibility.tsv',FEAS/'clinical_sample_manifest.tsv']
input_hashes = {str(p):digest(p) for p in inputs if p.is_file()}
dump(input_hashes,OUT/'resources/INPUT_SHA256.json')
for p in (SRC/'tables').glob('*'):
    if p.is_file(): shutil.copy2(p,OUT/'tables'/p.name)
for name in ['testable_class_state_bone_lung.tsv','all_class_state_bone_lung_feasibility.tsv','clinical_sample_manifest.tsv']:
    shutil.copy2(FEAS/name,OUT/'tables'/name)

genes = read(SRC/'tables/cohort_nominal_genes_ALL.tsv.gz')
programs = read(SRC/'tables/cohort_nominal_programs_ALL.tsv')
allprog = read(SRC/'tables/cohort_all_program_statistics.tsv.gz')
props = read(SRC/'tables/cohort_all_proportion_statistics.tsv')
cross = read(SRC/'tables/cohort_direction_consistency_summary.tsv')
customer = read(SRC/'tables/customer_program_cross_design_evidence.tsv')
other = read(SRC/'tables/gene_other_source_evidence.tsv.gz')
eligible = read(FEAS/'testable_class_state_bone_lung.tsv')
assert len(eligible)==10 and len(genes)==2013 and len(programs)==783
assert genes.P.lt(.05).all() and programs.P.lt(.05).all()
assert genes.FDR.lt(.05).sum()==0 and programs.FDR.lt(.05).sum()==0
assert eligible.nominal_up.sum()+eligible.nominal_down.sum()==len(genes)
assert set(genes.context_id) <= set(eligible.context_id)
assert not props.P.lt(.05).any()

# Include every eligible contrast, including null results, without evidence-based exclusion.
coverage=[]
for _,r in eligible.iterrows():
    for site in ['bone','lung']:
        items = str(r[f'{site}_cells_by_patient']).split(';')
        for item in items:
            donor,count=item.rsplit(':',1)
            coverage.append(dict(context_id=r.context_id,layer=r.layer,dataset=r.dataset,target=r.target,design=r.design,site=site,donor_record=donor,cells=int(count)))
coverage=pd.DataFrame(coverage)
save(coverage,'scheme1_patient_cell_coverage.tsv')
for (ctx,site),d in coverage.groupby(['context_id','site']):
    n=int(eligible.loc[eligible.context_id.eq(ctx),f'{site}_N'].iloc[0])
    assert d.donor_record.nunique()==n
for _,r in eligible.iterrows():
    d=coverage[coverage.context_id.eq(r.context_id)]
    b=set(d[d.site.eq('bone')].donor_record);l=set(d[d.site.eq('lung')].donor_record)
    assert b==l if r.design=='paired' else not b&l

program_counts=programs.groupby('context_id').size()
eligible['nominal_program_records']=eligible.context_id.map(program_counts).fillna(0).astype(int)
eligible['paired_direction_available_genes']=eligible.context_id.map(cross.set_index('context_id').paired_effect_available).fillna(0).astype(int)
eligible['paired_same_direction_genes']=eligible.context_id.map(cross.set_index('context_id').paired_mean_same_direction).fillna(0).astype(int)
eligible['paired_opposite_direction_genes']=eligible.context_id.map(cross.set_index('context_id').paired_mean_opposite_direction).fillna(0).astype(int)
save(eligible,'scheme1_all_eligible_comparisons.tsv')
lib=programs.copy();lib['library']=lib.program.str.split('|').str[0]
save(lib.groupby(['context_id','layer','dataset','target','library']).agg(nominal_programs=('program','size'),min_P=('P','min'),min_FDR=('FDR','min')).reset_index(),'scheme1_program_library_summary.tsv')

other['opposite_direction']=other.primary_effect*other.other_effect<0
other['zero_effect']=other.other_effect.eq(0)
other['valid_other_P']=other.other_P.notna()
other['other_nominal']=other.other_P.lt(.05)
osum=other.groupby(['primary_context','other_dataset']).agg(records=('feature','size'),same_direction=('same_direction','sum'),opposite_direction=('opposite_direction','sum'),zero_effect=('zero_effect','sum'),valid_other_P=('valid_other_P','sum'),other_nominal=('other_nominal','sum')).reset_index()
save(osum,'scheme1_other_source_consistency.tsv')

# Focus panels are selected by customer programs and matrix relevance; full tables remain intact.
focus=programs[(programs.program.str.startswith('customer|')) | (programs.target.eq('F4_C9') & programs.program.isin(['matrisome|NABA_COLLAGENS','matrisome|NABA_BASEMENT_MEMBRANES','matrisome|NABA_CORE_MATRISOME']))].copy()
save(focus,'scheme1_highlighted_programs.tsv')
font={'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False}
plt.rcParams.update(font)
labels=[]
for _,r in eligible.iterrows():
    prefix='S' if r.layer=='supported' else 'E'
    labels.append(f"{prefix}-{r.target.replace('original','')} | {r.dataset.replace('SCPCP0000','S')} | {'paired' if r.design=='paired' else 'independent'}")
y=np.arange(len(eligible))
fig,axes=plt.subplots(1,2,figsize=(14,7),gridspec_kw={'width_ratios':[1,1.35]})
axes[0].barh(y-.17,eligible.bone_N,height=.32,color='#326b9b',label='Bone')
axes[0].barh(y+.17,eligible.lung_N,height=.32,color='#ce6840',label='Lung')
axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].set_xlabel('Biological patients per comparison');axes[0].legend(frameon=False)
axes[1].barh(y,eligible.nominal_up,color='#ce6840',label='Higher in lung')
axes[1].barh(y,-eligible.nominal_down,color='#326b9b',label='Lower in lung')
axes[1].axvline(0,color='#555',lw=.7);axes[1].set_yticks(y,[]);axes[1].invert_yaxis();axes[1].set_xlim(-230,1000)
for i,(_,r) in enumerate(eligible.iterrows()): axes[1].text(max(r.nominal_up,0)+15,i,f"{int(r.nominal_up)} / {int(r.nominal_down)}",va='center',fontsize=9)
axes[1].set_xlabel('Nominal gene records (raw P < 0.05)');axes[1].legend(frameon=False,loc='lower right')
fig.suptitle('Scheme 1 | all eligible class/state expression comparisons',fontsize=15,x=.54)
fig.text(.02,.035,'S: supported layer; E: expanded layer. Paired N is the same patients at both sites. Zero results retained.\nNo listed gene contrast reaches FDR < 0.05. Classes/layers and parent/internal populations overlap; counts are not independent discoveries.',fontsize=9)
fig.tight_layout(rect=[0,.095,1,.94]);fig.savefig(OUT/'figures/01_all_eligible_comparisons.png',dpi=180);plt.close(fig)

focus=focus.sort_values(['target','program'])
fig,ax=plt.subplots(figsize=(13,6))
f_labels=[]
for i,(_,r) in enumerate(focus.iterrows()):
    f_labels.append(f"E-{r.target.replace('original','')} | {r.program.split('|')[1].replace('NABA_','').replace('_customer','').replace('_',' ')}")
    ax.plot([r.LOPO_effect_min,r.LOPO_effect_max],[i,i],color='#326b9b',lw=3,alpha=.55)
    ax.scatter(r.mean_lung_minus_bone,i,color='#326b9b',s=55,zorder=3)
    if r.paired_followup_N>0:
        ax.scatter(r.paired_followup_mean_effect,i,marker='D',color='#ce6840',s=55,zorder=3)
    ax.text(.59,i,f"P={r.P:.4f}  FDR={r.FDR:.3f}",va='center',fontsize=9)
ax.axvline(0,color='#555',lw=.8);ax.set_yticks(range(len(focus)),f_labels);ax.invert_yaxis();ax.set_xlim(-.24,.99)
ax.set_xlabel('Mean program effect: lung minus bone (log1p CP10k score)')
ax.set_title('Functional signals and paired-patient direction',fontsize=15,pad=18)
fig.text(.02,.03,'Blue: independent comparison; horizontal segment: leave-one-patient-out effect range (not a confidence interval).\nOrange diamond: paired mean (F1/C4 N=1; F4 N=2). F4-C9 lung coverage: 1 and 3 cells; no paired state comparison.\nDisplayed programs are customer panels plus selected matrix themes; all 783 nominal program records are retained in the accompanying tables.',fontsize=9)
fig.tight_layout(rect=[0,.14,1,.98]);fig.savefig(OUT/'figures/02_program_effects_and_discordance.png',dpi=180);plt.close(fig)

fig,ax=plt.subplots(figsize=(12,5.5))
for i,(_,r) in enumerate(eligible.iterrows()):
    d=coverage[coverage.context_id.eq(r.context_id)]
    for j,(site,col) in enumerate([('bone','#326b9b'),('lung','#ce6840')]):
        z=d[d.site.eq(site)]
        offset=np.linspace(-.08,.08,len(z)) if len(z)>1 else np.array([0])
        ax.scatter(z.cells,np.full(len(z),i)+(j-.5)*.22+offset,color=col,s=28,alpha=.75,label=site.title() if i==0 else None)
ax.set_xscale('log');ax.set_xlim(.8,max(coverage.cells)*1.4);ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Cells per biological patient at each site (log scale)');ax.set_title('Patient representation underlying each expression comparison',fontsize=14);ax.legend(frameon=False)
fig.tight_layout();fig.savefig(OUT/'figures/03_patient_cell_coverage.png',dpi=180);plt.close(fig)

report = '''骨肉瘤 CAF 解剖病灶部位差异分析
方案一正式结果 2026年10月12日

研究结论
按实际取材部位比较全部明确骨病灶与肺病灶，F1/F3/F4及固定内部状态中可评估的表达比较已完成。主要差异线索位于扩展F1、支持F3、扩展F4及内部F1-C4、F4-C9；现有数据未显示预设捕获比例的名义显著变化。功能差异需同时结合配对患者、其他来源及细胞覆盖解释。

数据与方法
沿用支持层1224、扩展层8174细胞及原始F类和内部标签。两层分别分析，不作为独立验证队列。骨、肺按已核实的解剖部位分组，不以是否发生肺转移、A/B1/B2或转移史完整性限制入组。其他、混合和不明确部位单列。转移史、治疗、谱系/CNV及双细胞记录保留。
本阶段分析F1/F3/F4主类及全部已有内部状态。按来源分别，以独立患者为统计单位；同一患者同部位先聚合，有两部位的患者进入配对比较，并从该端点非配对比较中移除。基因表达使用原始计数患者部位伪合并、完整文库归一化的log2(CP10k+0.1)；功能评分采用完整冻结面板的平均log1p(CP10k)。所有效应定义为肺减骨。
比例分别使用全部来源捕获细胞、相应候选层、内部群父类作为分母。保留零分子；零分母不可评估。缺有效重复或全零/无变化不补造P。沿用已验收的患者置换或符号翻转检验、同类家族FDR。有效原始P<0.05作为名义候选标准，FDR作为证据报告。ORA单独表示基因列表富集模型统计。
本次明确采用用户确认的方案一，复用Phase73 strict_v2与Phase74完整计算，新增统一结果表、患者覆盖图和跨设计证据图；不重复运行相同检验，不改变原始P、FDR或检验家族。

全部可评估比较
共10项表达比较具备有效重复，另外29项仅作方向描述；其余登记组合一侧或两侧群体缺失。10项不是10个独立验证研究。完整名义基因2013条、程序783条，按群体/对比计记录，不是唯一基因数。各项见下表。
'''
tab=eligible[['layer','dataset','target','design','bone_N','lung_N','nominal_up','nominal_down','nominal_program_records']]
report += tab.to_string(index=False)+'\n'
report += '''
扩展F1及内部C4具有抗原呈递相关名义差异
SCPCP000017扩展F1独立比较为骨8、肺3名患者，153条名义上调和9条下调基因。apCAF抗原呈递程序平均效应+0.3343，P=0.024242，FDR=0.393467。唯一配对患者效应-0.1196，与独立组相反。内部C4为152条上调、10条下调；抗原呈递程序效应+0.3636，P=0.018182，FDR=1，唯一配对患者效应-0.1106。父类与C4共享141条名义基因且同向，属于嵌套群体信号。支持层F1的2对2比较没有名义差异基因。

支持F3显示广泛表达差异但功能主题尚不集中
SCPCP000017支持F3骨6、肺2名患者，758条名义上调基因、0条下调，37条名义程序均来自Reactome。客户预设面板、Hallmark和Matrisome患者级程序未出现名义差异。保持全转录组观察，不强行归结为胶原沉积或免疫抑制。扩展F3为骨2、肺1名患者，保留方向而不作推断；两层成员不同，不能相互替代。

扩展F4的独立组差异与配对方向存在分歧
SCPCP000017扩展F4骨11、肺2名患者，633条名义上调、161条下调。炎症样程序平均效应+0.3147，P=0.012821，FDR=0.392392。两名配对患者该程序均为负向，平均-0.1421。794条名义基因在配对平均效应中同向257条、反向496条、零41条。其他来源SCPCP000023的3对3比较没有名义基因，主来源794条候选在那里同向436条、反向350条、零8条。保留独立组信号及矛盾，当前不能概括为普遍的肺侧F4功能增强。

内部F4-C9保留基质功能候选线索
SCPCP000017 F4-C9骨5、肺2名患者，123条名义上调、14条下调。胶原程序平均效应+0.4820，基底膜+0.4112，核心基质+0.2271；均P=0.047619、FDR=1。137条名义基因仅6条与父F4名义列表重合，保留该内部状态自身结果。肺侧患者分别仅1与3个细胞，没有同状态有效配对表达比较。

阴性与配对结果完整保留
SCPCP000017 F4及F4-C1均有2对患者，可计算精确配对效应与P，当前无名义基因。SCPCP000023的F4-C3有2对2独立患者，也没有名义基因。相应最小可达P受离散样本量限制，不据阴性否定已观察到的变化方向。原始全细胞口径有71条有效比例P，均未达到P<0.05；群内表达差异不能代替群体比例变化。

讨论与证据边界
全部名义基因及患者级名义程序均未通过相应FDR阈值，应称名义显著或探索性信号。F1/F4配对与非配对方向的分歧会影响跨患者概括，已在结果和图中直接保留。主要来源肺侧F1有3名患者，其中2名仅2与1个细胞；F4肺侧两名患者分别1与642细胞；F4-C9肺侧共4细胞。保留这些对象而不增加入口硬阈值，同时据实际覆盖解释。
患者数、来源数和细胞数分别记录；不能把重叠层、父子类或跨来源患者记录简单合并为独立重复。部位关联不等于转移发生的因果作用，程序表达不等于蛋白功能、物理基质屏障或已验证新CAF类型。新结果作为10月6日功能分类与胶原基质报告的部位比较补充，不覆盖既有报告。

交付与进度
方案一计算结果整理、全可评估对比登记与图表已完成。当前没有后台分析作业；流式标志筛选保持后置，未按表面蛋白或抗体可用性筛除队列差异。所有名义上下调、完整程序/比例/ORA、跨患者与来源方向及不可评估原因随统计表保存。完整未筛选单基因检验仍位于Phase73 strict_v2，来源和哈希见INPUT_MANIFEST及既有档案。
图1 全部可评估表达比较与名义基因；图2 功能效应及配对反向证据；图3 各比较患者细胞覆盖。图2逐患者剔除范围是描述性效应范围，不是置信区间。
'''
(OUT/'delivery/骨肉瘤CAF_骨肺病灶对比_方案一结果_20261012.txt').write_text(report,encoding='utf-8')
shutil.copy2(PROJECT/'phase73_F134_bone_lung_20261009/strict_v2/resources/INPUT_MANIFEST.json',OUT/'resources/PHASE73_INPUT_MANIFEST.json')
assert all(digest(Path(p))==h for p,h in input_hashes.items())
qa=dict(status='computed_checks_passed_visual_review_pending',new_inferential_tests=0,source_inputs_unchanged=True,nominal_gene_records=len(genes),nominal_program_records=len(programs),gene_FDR_significant=int(genes.FDR.lt(.05).sum()),program_FDR_significant=int(programs.FDR.lt(.05).sum()),eligible_expression_comparisons=len(eligible),patient_counts_match=True,paired_and_independent_patient_sets_verified=True,all_negative_eligible_comparisons_retained=True,flow_inputs_used=False,visual_review_pending=True,figures=3,completed_at=datetime.now(timezone.utc).isoformat())
dump(qa,OUT/'QA/COMPUTED_CHECKS.json')
(OUT/'ANALYSIS_FINISHED').write_text('Result assembly completed; final acceptance requires visual review.\n')
print(json.dumps(qa,ensure_ascii=False))
