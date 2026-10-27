from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,html,zipfile,datetime,itertools,shutil
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase120_patient_supported_CAF_review_20261026';I=R/'phase119_Bone_Lung_functional_subgroups_20261024'
def rd(n):return pd.read_csv(O/'tables'/n,sep='\t')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
g=rd('source_candidate_gate.tsv');d=rd('functional_audit.tsv');e=rd('patient_effects.tsv.gz');sp=rd('node_donor_support.tsv');m=rd('cell_identity_audit.tsv.gz');nodes=pd.read_csv(I/'tables/hierarchical_nodes.tsv',sep='\t');rv=pd.read_csv(I/'tables/target_collagen_candidate_review.tsv',sep='\t');lit=json.loads((O/'resources/LITERATURE_CROSSCHECK.json').read_text(encoding='utf-8'));qa=json.loads((O/'AUDIT_COMPLETE.json').read_text());run=json.loads((O/'RUN_RECEIPT.json').read_text());assert run['exit_code']==0
assert len(nodes)==65 and nodes[nodes.level.eq('parent')].fit_cells.sum()==6423
assert sp.groupby('node_id').cells.sum().reindex(nodes.node_id).to_numpy().tolist()==nodes.fit_cells.tolist()
assert g.qualified.sum()==0
focus=['Bone_07','Bone_07_s1','Lung_04','Lung_04_s2']
recheck=e[e.node_id.isin(focus)&e['mode'].eq('main_min5')&e.endpoint.isin(['program__mCAF_matrix','audit__FNDC1_paper','audit__collagen_processing_nonCOL'])];maxerr=0;tests=0
for key,a in recheck.groupby(['node_id','source','mode','endpoint']):
 v=a.effect.to_numpy();n=len(v);r=d[(d.node_id==key[0])&(d.source==key[1])&(d['mode']==key[2])&(d.endpoint==key[3])].iloc[0];assert n==r.N
 if n<2:assert pd.isna(r.P);continue
 vals=[abs(sum(s*x for s,x in zip(signs,v))) for signs in itertools.product([-1,1],repeat=n)];p=sum(t>=abs(v.sum())-1e-10 for t in vals)/len(vals);maxerr=max(maxerr,abs(p-r.P));tests+=1
assert maxerr<1e-12
from scipy.stats import false_discovery_control
ferr=0
for _,a in d[d.P.notna()].groupby(['site','source','level','mode','family']):ferr=max(ferr,float(np.max(abs(false_discovery_control(a.P.to_numpy())-a.FDR))))
assert ferr<1e-12
rows=[]
for n in nodes.itertuples():
 choices=g[g.node_id.eq(n.node_id)].copy();actual=sp[sp.node_id.eq(n.node_id)].groupby('source').cells.sum();choices['source_cells']=choices.source.map(actual).fillna(0)
 a=choices[choices.source_cells.gt(0)].sort_values(['N_min5','N_identity5','source_cells','source'],ascending=[False,False,False,True]).iloc[0]
 b=d[d.node_id.eq(n.node_id)&d.source.eq(a.source)&d['mode'].eq('main_min1')];N1=int(b.N.max())
 reason=[]
 if a.N_min5<6:reason.append('有效患者不足6')
 if pd.isna(a.nominal_endpoints):reason.append('无满足人数/方向/名义P的目标程序')
 if not a.processing_support:reason.append('非胶原加工面板支持未达规则')
 if a.N_identity5<6:reason.append('身份敏感性有效患者不足6')
 if a.source_patient_dominance>.5:reason.append('单患者细胞占比>50%')
 if not a.stable:reason.append('种子稳定性<0.8')
 tier='优先保留的来源内基质线索，未达最终门槛' if n.node_id=='Lung_04' else '骨侧基质线索，患者数不足' if n.node_id=='Bone_07' else '细分探索线索，非最终候选' if n.node_id in focus else '本轮未进入最终ECM候选'
 rows.append(dict(父群或子群=n.node_id,层级=n.level,细胞=n.fit_cells,覆盖记录=n.source_patients,评估来源=a.source,最少1细胞有效患者=N1,最少5细胞有效患者=int(a.N_min5),最少10细胞有效患者=int(a.N_min10),身份敏感性有效患者=int(a.N_identity5),结论=tier,未达标原因='；'.join(reason)))
overview=pd.DataFrame(rows);overview.to_csv(O/'tables/node_decisions_ZH.tsv',sep='\t',index=False,encoding='utf-8')
focus_summary=overview[overview['父群或子群'].isin(focus)]
endpoints=['program__collagen_biosynthesis','program__collagen_crosslinking','program__fibrillar_I_III_V','program__mCAF_matrix','audit__collagen_processing_nonCOL','audit__FNDC1_paper','audit__MYH11_paper','audit__CD36_paper','audit__C7_paper','audit__osteogenic_identity']
detail=d[d.node_id.isin(focus)&d['mode'].isin(['main_min5','main_min10','identity_min5'])&d.N.gt(0)&d.endpoint.isin(endpoints)][['node_id','source','mode','endpoint','N','mean_effect','median_effect','positive_patients','P','FDR','LOO_min_mean']]
parents=overview[overview['层级'].eq('parent')]
cnv=pd.read_csv(I/'tables/node_inherited_CNV_records.tsv',sep='\t');cnv=cnv[cnv.node_id.isin(focus)]
hold=pd.read_csv(I/'tables/holdout_statistics.tsv',sep='\t');hold=hold[hold.node_id.isin(focus)&hold.endpoint.eq('program__mCAF_matrix')&hold.background.isin(['within_site','within_parent'])][['node_id','source','background','N','mean_effect','positive_patients','P','FDR']]
source_roster=m[~m.doublet].groupby(['source','site']).agg(细胞=('cell_id','size'),来源内患者=('donor','nunique'),样本记录=('sample','nunique')).reset_index();source_roster.to_csv(O/'tables/source_patient_roster_counts.tsv',sep='\t',index=False)
assert not m.groupby(['source','sample']).donor.nunique().gt(1).any();assert not m.groupby(['source','donor']).site.nunique().gt(1).any()
lung=d[d.node_id.eq('Lung_04')&d.source.eq('NCI')&d['mode'].eq('main_min5')&d.endpoint.eq('audit__FNDC1_paper')].iloc[0];assert lung.N==7 and lung.positive_patients==7
for f,h in json.loads((O/'INPUT_HASHES.json').read_text()).items():assert sha(f)==h
extra=[R/'phase102_NCI_author_CAF_20261020/tables/author_CAF_cell_registry.tsv.gz',I/'tables/node_inherited_CNV_records.tsv',R/'phase118_reference_based_broad_classes_20261023/references/SOURCE_MANIFEST.json',R/'phase118_reference_based_broad_classes_20261023/references/four_class_report_extracted.txt',R/'phase118_reference_based_broad_classes_20261023/references/CD36_paper_extracted.txt']
dump('resources/REFERENCE_PROVENANCE.json',{str(f):sha(f) for f in extra})
def tab(x):return '<div class="scroll">'+x.to_html(index=False,border=0,float_format=lambda x:f'{x:.6g}',na_rep='不可评估')+'</div>'
intro='已完成固定6,633细胞队列中29个父群及36个子群的患者级核查。按预先冻结的患者支持和身份敏感性规则，最终计算候选为0；Lung_04是当前证据最强的来源内基质线索，Bone_07保留为骨侧线索。阴性结论指未达到最终准入规则，不表示没有基质功能。'
threshold='最低6位来源内去重、可配对患者；每位目标群及对照各至少5个非软件双细胞阳性细胞；至少75%同向，并至少一个目标基质/胶原程序原始P<0.05。推荐10位及以上作为更高支持等级。对照为同患者同病灶同来源的其余候选细胞；子群主要与同父群其余细胞比较。6是统计可评估性的最低操作门槛，不是充分功效或文献通用标准：双侧精确符号翻转最小P为2/2^N，N=5时0.0625，N=6时0.03125。5细胞是沿用既有敏感性分析的低限，本轮同时检查10细胞；它不是高精度保证。'
rules='额外要求：合格来源内单患者占该群细胞不超过50%；既有种子Jaccard至少0.8；去掉任一患者后目标效应均值仍为正；不含COL家族的胶原加工面板至少75%同向且均值及留一均值为正，避免共享胶原基因重复计证据。身份敏感性层仍需6位有效患者及75%同向，不强求重复一次名义显著。若另一个来源有至少3位min5有效患者且目标程序75%同向，可标记来源重复支持，但不能称外部独立验证。所有数值规则均为本项目操作规则，不是领域共识阈值。'
findings='Lung_04主分析min5为7位，交联、纤维型胶原、基质面板均7/7同向；FNDC1论文六基因面板也7/7同向，效应+0.317835，P=0.015625，FDR=0.06875。非COL加工面板6/7同向，效应+0.228204，P=0.03125，FDR=0.105769。这些是名义显著的探索性信号。其MYH11、CD36和C7论文面板没有相应整体增强，故更接近基质/成骨重叠方向，不能称CD36型或文献亚群完全复现。'
attrition='Lung_04身份敏感性降至5位：NCI_PATIENT_6546目标细胞34→24、对照16→4；NCI_PATIENT_7978目标14→13、对照5→3。因此这两位主要因合格对照不足退出配对，并非目标群消失。其余5位的基质、纤维型胶原、FNDC1及非COL加工面板仍5/5升高；P=0.0625符合5位时的精确检验分辨率，不应解读为功能被推翻。身份门槛未达标与实质污染证据应分开表述。'
bone='Bone_07主分析min5仅4位，交联和基质均4/4同向；身份敏感性3位。Bone_07_s1同父群对照min5为3位，身份敏感性2位；其非COL加工面板主分析均值为负，因此不能仅凭共享胶原面板认定广泛胶原生产增强。Lung_04_s2同父群min5仅1位。保留方向线索，但不升格为最终亚群。Bone_07_s1实际覆盖5位来源内患者记录，之前对话中的4位是原min1同父群有效配对N，不应写成总覆盖人数。'
identity='身份敏感性使用原成纤维支持层；NCI进一步限定作者CopyKAT diploid，并排除PTPRC与LST1/TYROBP/CD3D/CD3E任一共表达或PECAM1与VWF/CDH5共表达的细胞。所有规则仅用于敏感性切片，未删除原细胞或改身份。共表达可能受环境RNA影响，不能当作已证实污染；RGS5单基因不作为删除依据。diploid或single_method_no_alteration也不是正常基因组或CAF身份的证明。'
literature_boundary='用户四分类报告自身将iCAF/myCAF/mCAF/osteogenic定义为整合框架，不是一篇论文确立的统一四分法。论文FNDC1型是ECM参照，MYH11偏血管相关，不能简单等同myCAF。Liu2021的C2_CAF与本项目F2_C3只是不同编号；F2_C3仍需依赖原项目字典及细胞对应，不能从名称推断一致。LOX研究复用了早期GSE152048/GSE162454；多篇文章不能累计为独立患者支持。文献中的68.35%也不作为本队列的目标比例。'
methods='保留Phase119全部细胞ID和层级，主分析6,423细胞（骨2,788，肺3,635），210个原软件双细胞阳性保留在总名单。未重聚类、未重做总体骨肺差异、未画科学图。程序采用既有log1p(CP10k)均分及原秩评分，新增论文面板在同一未整合表达空间评分。来源内同患者目标均值减对照均值，每患者等权；精确双侧符号翻转（N≤16），更大N用65,536次蒙特卡洛。N<2或全零不造P。BH按病灶、来源、父子层级、分析层和终点家族划分；完整P/FDR保留。新增audit面板作为独立的复核面板家族，含全部11个预设面板，未只对显著项校正。'
limits='这是发现结果上的事后核查，门槛在本轮计算前冻结，但并非原研究前瞻预注册。来源内donor合并样本记录，同一患者多样本不计作多患者；不同来源不直接累加为全球独立N。既有留出结果没有建立独立重复支持。RNA支持的是表达程序，不能证明胶原分泌量、真实纤维沉积、矿化、耐药或转移机制。没有为使某群达标降低阈值，也没有把所有患者特异群合并以凑人数。'
actions='当前应保留Lung_04作为最优先后续验证对象，优先补足清晰身份的同患者对照或独立患者/来源重复；Bone_07需增加有效患者覆盖。两个小子群不作为本轮最终候选。对现有数据的本轮核查、文献交叉和筛选已完成；外部生物学验证为后续研究建议，不是后台仍在执行的任务。'
css='body{font:16px/1.75 "Microsoft YaHei",sans-serif;color:#243543;margin:24px auto;max-width:1260px;padding:24px}h1{font-size:28px}h2{font-size:21px;border-bottom:1px solid #bdcbd1;margin-top:30px}.lead{background:#e9f3f3;padding:18px;border-left:5px solid #287780}.scroll{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%}td,th{border:1px solid #dce3e6;padding:8px;vertical-align:top}th{background:#e9f1f5}a{color:#126880}code{overflow-wrap:anywhere}'
body='<h1>CAF候选患者支持与文献交叉核查</h1><p>内部研究核查｜Phase120｜2026-10-26｜科学图件暂停</p><div class="lead">'+intro+'</div>'
sections=[('1. 推荐并冻结的最终候选门槛',threshold+'\n'+rules),('2. 重点候选的患者支持',findings+'\n'+bone),('3. 严格身份层为何少了两位患者',attrition+'\n'+identity),('4. 文献对应与证据边界',literature_boundary),('5. Methods',methods),('6. 适用范围与后续研究',limits+'\n'+actions)]
for title,txt in sections:
 body+='<h2>'+title+'</h2>'+''.join('<p>'+html.escape(p)+'</p>' for p in txt.split('\n'))
 if title.startswith('2.'):body+=tab(focus_summary)
 if title.startswith('4.'):
  for a in lit:body+='<p><b>'+a['id']+'：<a href="'+a['url']+'">'+html.escape(a['title'])+'</a></b><br>'+html.escape(a['evidence']+'。'+a['finding']+' '+a['boundary'])+'</p>'
body+='<h2>7. 所有父群的准入核查</h2>'+tab(parents)+'<h2>8. 完整数值与过程记录</h2><p>效应为程序均分之差，不是倍数变化。下表N为有效配对患者数；来源内总覆盖数另列，父子细胞不能相加。</p>'
for title,df in [('重点功能与文献面板（含10细胞及身份敏感性）',detail),('全部65节点准入结论',overview),('重点父群/子群的原CNV记录',cnv),('既有留出患者结果',hold),('主分析来源及患者记录',source_roster),('重点群逐患者细胞与身份计数',sp[sp.node_id.isin(focus)])]:body+='<details><summary>'+title+'</summary>'+tab(df)+'</details>'
body+='<p>验收：15,792项既有功能对比重新计算一致；重点精确P独立枚举复核'+str(tests)+'项；全部新BH家族复核通过；输入哈希未变。原始数据、代码、全统计表及过程日志保存在服务器 <code>'+str(O)+'</code>。本报告未对外发送。</p>'
report=O/'delivery/CAF_patient_support_literature_review.html'
if report.exists():
 backup=O/'resources/report_v1_before_source_tie_fix';backup.mkdir(exist_ok=True)
 for f in (O/'delivery').glob('*'):
  if f.is_file() and not (backup/f.name).exists():shutil.copyfile(f,backup/f.name)
 dump('resources/REPORT_SOURCE_TIE_FIX.json',dict(reason='When all min5 N are zero, select an actually represented source and use source cell counts as tie-breaker; avoids showing nonexistent NCI for Lung_11',statistics_changed=False))
report.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CAF患者支持与文献核查</title><style>'+css+'</style>'+body+'</html>',encoding='utf-8')
txt='CAF患者支持与文献交叉核查\n2026-10-26 内部研究核查\n\n'+intro+'\n\n'+'\n\n'.join(t+'\n'+v for t,v in sections)+'\n\n所有父群及子群\n'+overview.to_string(index=False)+'\n\n重点完整统计\n'+detail.to_string(index=False)+'\n\n文献\n'+'\n\n'.join(a['id']+' '+a['title']+'\n'+a['url']+'\n'+a['finding']+'\n'+a['boundary'] for a in lit)
(O/'delivery/CAF_patient_support_literature_review.txt').write_text(txt,encoding='utf-8')
for f in (O/'delivery').glob('*'):
 if f.suffix in ['.txt','.html']:assert '\ufffd' not in f.read_text(encoding='utf-8')
dump('ACCEPTANCE.json',dict(status='PASS_ANALYSIS',cells=6633,main_cells=6423,nodes=65,parents=29,children=36,original_input_hashes_unchanged=True,existing_tests_recomputed=qa['recomputed_comparisons'],independent_exact_tests=tests,max_exact_error=maxerr,max_BH_error=ferr,final_candidates=0,priority_followup='Lung_04',scientific_figures_created=False,remaining=['local delivery and visual QA']))
dump('STATUS.json',dict(status='analysis_complete_delivery_pending',priority='Lung_04',final_candidates=0,report=str(report)))
print('REPORT_READY',sha(report),report.stat().st_size,flush=True)
