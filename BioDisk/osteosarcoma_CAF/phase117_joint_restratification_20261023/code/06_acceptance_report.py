from caf117_common import *
import itertools,html,zipfile
from scipy.stats import false_discovery_control
event('06_acceptance_report','started')
for k,n in [('01','FREEZE_COMPLETE.json'),('02','STRUCTURE_COMPLETE.json'),('03','CORRESPONDENCE_COMPLETE.json'),('04','FUNCTION_COMPLETE.json'),('05','HOLDOUT_COMPLETE.json')]:
 rec=json.loads((O/'receipts'/(k+'.json')).read_text());assert rec['exit_code']==0 and (O/n).exists()
X,m,genes=load();lab=read('joint_cell_labels.tsv.gz');assert len(lab)==10510 and lab.cell_id.is_unique and lab.joint_cluster.notna().all();assert lab.cell_id.equals(m.cell_id)
old=pd.read_csv(R/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');om=m[m.cohort.eq('baseline')].set_index('cell_id');assert len(om)==8174 and (om.original_state==old.loc[om.index].subgroup_id).all()
counts=read('old33_by_joint_counts.tsv');assert int(counts.iloc[:,1:].to_numpy().sum())==8148
split=read('holdout_donor_split.tsv');assert split.groupby('donor_key').validation_split.nunique().max()==1
primary=read('joint_function_statistics.tsv');hold=read('holdout_function_statistics.tsv');maxbh=0.
for df in [primary,hold]:
 assert df.loc[df.N.lt(2),'P'].isna().all();assert df.loc[df.test.isin(['coverage_missing','all_zero_not_tested']),'P'].isna().all()
 for _,g in df.groupby(['source','arm','family']):
  g=g[g.P.notna()].sort_values('P');n=len(g)
  if not n:continue
  bh=np.minimum.accumulate((g.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1);er=float(np.max(abs(bh-g.FDR.to_numpy())));assert er<1e-10;maxbh=max(maxbh,er)
checked=0;maxerr=0.;maxmc=0.;signs={}
effects=read('joint_function_patient_effects.tsv.gz')
test=primary[primary.source.isin(['NCC','NCI'])&primary.arm.isin(['singlet','singlet_min5'])&primary.family.isin(['functional_program','original_signature','target_gene'])]
for key,g in test.groupby(['source','arm','joint_cluster']):
 src,arm,cl=key;ee=effects[(effects.source==src)&(effects.arm==arm)&(effects.joint_cluster==cl)];wide=ee.pivot(index='donor',columns='endpoint',values='effect');n=len(wide)
 if n>=2 and n not in signs:signs[n]=np.array(list(itertools.product((-1.,1.),repeat=n)))
 for _,row in g.iterrows():
  if n==0:assert row.N==0;continue
  v=wide[row.endpoint].to_numpy();assert n==row.N
  if np.isfinite(v).all():assert abs(v.mean()-row.mean_effect)<1e-10;assert (v>1e-12).sum()==row.positive_patients
  if pd.isna(row.P):continue
  exact=float(np.mean(abs(signs[n]@v)>=abs(v.sum())-1e-10));error=abs(exact-row.P)
  if row.test=='exact_signflip':assert error<1e-10;maxerr=max(maxerr,error)
  else:assert error<max(.001,6*np.sqrt(exact*(1-exact)/65536));maxmc=max(maxmc,error)
  checked+=1
hashes=json.loads((O/'resources/INPUT_HASHES.json').read_text())
for h in hashes:assert sha(h['path'])==h['sha256']
qa=dict(status='PASS_COMPUTATION_AND_EVIDENCE_SCOPE',cells=10510,fit_cells=10274,baseline_labels_unchanged=8174,transition_counts_reconciled=8148,input_hashes_unchanged=len(hashes),independent_new_function_tests=checked,max_exact_P_error=maxerr,max_MC_difference=maxmc,max_BH_error=maxbh,donor_holdout_disjoint=True,all_required_native_exit_codes_zero=True,figure_queue='paused_by_user',identity_claim='Computational joint structure and functional candidates, not confirmed CAF identity or deposition mechanism')
dump('FINAL_ACCEPTANCE.json',qa)
structure=json.loads((O/'STRUCTURE_COMPLETE.json').read_text());pilot=json.loads((O/'PILOT_ACCEPTANCE.json').read_text());holdstatus=json.loads((O/'HOLDOUT_COMPLETE.json').read_text());routes=read('old33_joint_correspondence.tsv');clustersummary=read('joint_cluster_origins_and_stability.tsv');route_counts=routes.route.value_counts().to_dict();c3=routes[routes.old_state.eq('F2_C3')].iloc[0];focus_cluster=c3.leading_joint_cluster
core=['fibrillar_I_III_V','collagen_biosynthesis','collagen_crosslinking','myCAF_customer','F2C3_collagen_support','Osteogenic_customer','signature__F2_C3'];fc=primary[primary.joint_cluster.eq(focus_cluster)&primary.source.isin(['NCC','NCI','SCPCP000017','SCPCP000023'])&primary.endpoint.isin(core)&primary.arm.eq('singlet')];fs=primary[primary.joint_cluster.eq(focus_cluster)&primary.source.isin(['NCC','NCI'])&primary.endpoint.isin(core)&primary.arm.eq('singlet_min5')];hc=hold[hold.joint_cluster.eq(focus_cluster)&hold.endpoint.isin(core)]
css='body{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#203440;background:#f3f6f8;margin:0}main{max-width:1200px;background:white;margin:28px auto;padding:36px}h1{font-size:29px}h2{color:#255367;border-bottom:2px solid #dce5e9;margin-top:34px}.lead{background:#eaf4f4;border-left:5px solid #287d86;padding:18px}.note{font-size:14px;color:#526974}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:16px 0}th{background:#254e60;color:white;padding:8px;text-align:left}td{padding:8px;border-bottom:1px solid #dbe4e8;vertical-align:top}tr:nth-child(even){background:#f3f7f9}details{margin:15px 0}code{overflow-wrap:anywhere}a{color:#166f8e}input{padding:9px;width:65%}button{padding:9px;margin:6px}'
def tab(d):return '<div class="scroll">'+d.to_html(index=False,border=0,escape=True,float_format=lambda v:f'{v:.6g}',na_rep='不可评估')+'</div>'
def esc(x):return html.escape(str(x))
cols=['source','joint_cluster','endpoint','cells','N','mean_effect','positive_patients','P','FDR']
mainroutes=routes[['old_state','working_name','original_all_cells','original_fit_cells','leading_joint_cluster','leading_fraction','second_joint_cluster','second_fraction','route','patient_balanced_leading_cluster','patient_balanced_leading_fraction']]
body=f'''<h1>CAF扩展层共同重分层与原33亚群对应</h1><p class="note">内部报告｜Phase117｜2026-10-23｜先冻结方案再顺序实施；科学图件队列保持暂停。</p><div class="lead">原8,174个细胞与NCC 806、NCI 1,530个新增细胞组成10,510个固定扩展资源。主结构拟合使用10,274个非软件双细胞阳性细胞，得到{structure['clusters']}个共同计算群。原33亚群通过真实原细胞ID追踪到共同新群，形成{route_counts.get('mostly_retained',0)}个主要保留、{route_counts.get('split',0)}个拆分、{route_counts.get('dispersed',0)}个分散对应；这些是预设比例规则下的结构描述，不等同于生物学身份确认。</div>
<h2>1. 方案与输入基准</h2><p>方案在执行前写入服务器PLAN_ZH.txt，并登记哈希与PROTOCOL_V1.json。原分类、原报告及原始表达未覆盖；新增数据不先通过F1-F5标签门槛。原26个、NCI210个软件双细胞阳性记录保留在总表，并进入全细胞敏感性。NCC采用修订后的7位来源内患者、NCI采用既已接受的17位分组；不将跨研究记录直接当作已验证独立N。</p>{tab(read('frozen_source_counts.tsv'))}
<h2>2. 先比较表示，再共同重分层</h2><p>比较未整合PCA、按来源Harmony theta=1、theta=2。主表示选择为 <b>{esc(structure['selected_representation'])}</b>。以原高/中置信主要谱系邻域保持作为操作性门槛，同时记录同谱系跨来源邻域比例；没有以最大来源混合或想要的群数选参数。</p>{tab(pd.DataFrame(pilot['decisions']))}{tab(read('pilot_representation_metrics.tsv'))}
<p>主分辨率0.5，敏感性0.3、0.8及另一随机种子；同一群数不是验收要求。未整合、全细胞、去除原标记/核心胶原收缩面板均保留为敏感性。所有主群采用J编号，不强制继承F名。原细胞单独计算的控制分群与共同新群ARI={structure['baseline_to_joint_ARI']:.3f}，用于描述加入新数据后的结构变化。</p>{tab(clustersummary)}
<h2>3. 原33亚群的真实细胞去向</h2><p>主要保留定义为原亚群至少70%细胞进入同一新群；若不满足且至少两个新群各承载20%，记为拆分，其余记为分散。表中同时给出第二去向与患者等权结果；“汇合”指同一新群的旧细胞中至少两个原亚群各占20%，不自动改写生物学命名。新来源富集不等于新型CAF。</p>{tab(mainroutes)}
<h2>4. F2-C3及基质相关亚群的重点对应</h2><p>原F2-C3共{int(c3.original_all_cells)}个细胞，主拟合层{int(c3.original_fit_cells)}个，其中{c3.leading_fraction:.1%}进入{esc(focus_cluster)}，第二去向为{esc(c3.second_joint_cluster)}（{c3.second_fraction:.1%}），按预设规则为{esc(c3.route)}。患者等权的首要去向为{esc(c3.patient_balanced_leading_cluster)}，比例{c3.patient_balanced_leading_fraction:.1%}。本轮优先根据这些原细胞去向判断对应，再查看新细胞贡献和功能，避免单纯按最近表达标签归类。</p>
{tab(clustersummary[clustersummary.joint_cluster.eq(focus_cluster)])}
<p>以下展示F2-C3首要对应群在不同来源内的患者级功能差异。Δ为目标新群与同来源、同患者/部位其他候选细胞的平均log1p(CP10k)评分差，经部位、患者等权汇总；不是倍数变化。N为可评估患者，cells不是检验样本量。小于10例按探索性解读；完整正、负、未通过及不可评估结果保留。</p>{tab(fc[cols])}
<h3>最低细胞数敏感性</h3><p>下表要求每患者/部位两侧均至少5个非软件阳性细胞。人数下降或不显著按实际结果保留，不用主层结论代替。</p>{tab(fs[cols])}
<h2>5. 留出患者的发现与验证</h2><p>按来源分层留出20%患者，小于3位患者的来源不切分。发现集独立选HVG、拟合PCA及整合，独立分群，再用训练细胞交集对应共同J群。留出细胞只投射到发现表示，不参与发现特征选择和拟合。本轮发现细胞{holdstatus['discovery_cells']}个、留出细胞{holdstatus['heldout_cells']}个；留出来源内人数为{esc(holdstatus['heldout_source_patient_records'])}。</p><p>Harmony查询投射采用同来源30近邻学习位移插值，再作15近邻标签预测；这是明确记录的近似投射，不是Harmony原生查询API。与全数据标签的一致率仅作结构诊断；患者级功能方向在留出分组内另外计算。该验证不等于独立外部实验。</p>{tab(read('heldout_structural_agreement.tsv'))}
<p>F2-C3首要对应J群的留出患者结果：</p>{tab(hc[cols])}
<h2>6. Methods与统计口径</h2><p>共用测量基因空间为13,637个唯一符号基因，表达为完整文库log1p(CP10k)。按来源感知选择目标3,000 HVG，排除明确技术/性别特征；主分析保留胶原、收缩等生物学信号。尺度上限10、PCA最多50维、15近邻、Leiden主分辨率0.5、固定种子20261023，线程4。原8174控制使用同一冻结HVG列表和参数；留出发现集重新在发现细胞选择HVG。来源校正不使用临床部位标签。</p><p>功能检验沿用原程序、原通路基因成员及33亚群标记组合，用未整合表达计算。按来源、患者、部位内目标新群对其余候选细胞比较，先部位等权再患者等权；N≤16精确双侧符号翻转，否则固定种子65536次蒙特卡洛加观测值。功能程序、归档通路、原亚群签名、目标单基因分别在来源×分析层内BH校正。原始P&lt;0.05为名义显著，FDR&lt;0.05才称校正后显著；N&lt;2、全零或覆盖不足不造P。保留全细胞、非阳性、最低5细胞三层。</p>
<h2>7. 解释边界与候选推进</h2><p>共同计算群的来源/患者主导程度、重复分群稳定性、旧谱系组成和去标记敏感性均单列。功能与标记参与主分群，因此全数据统计是计算分组条件下的证据；留出患者用于补充检验，不能用细胞数替代生物学重复。混合谱系、成骨或壁细胞重叠不静默删除，也不直接命名为纯CAF。NCC实际部位和跨研究去重仍沿用既有边界，NCI单核平台与来源效应不能完全分离。本轮不作新增骨肺临床比较；RNA程序增强不是实际纤维沉积或功能机制确认。</p><p>后续应从原亚群细胞去向明确、同时获得新增患者贡献且功能方向一致的J群中确定重点。原亚群若拆分，分别评价各分支；若汇合，保留组成和竞争解释；缺少旧对应的新来源富集群保留候选，不直接宣布新亚型。</p>
<h2>8. 验收与过程记录</h2><p>主清单10,510细胞守恒，原8,174标签与输入哈希未变；旧亚群主拟合对应矩阵总和8,148。新增主功能检验独立复核{checked}项，最大精确P误差{maxerr:.3g}，所有BH家族重算最大误差{maxbh:.3g}。步骤01–05原生退出码和完成标记均通过。过程记录包括events.jsonl、逐步日志、代码哈希、每次尝试receipt及STATUS；本报告通过仅表示本轮计算与交付范围完成。</p><p><a href="CAF_joint_full_tables.html">完整对应、功能、稳定性及留出结果附表</a></p><p class="note">服务器目录：<code>{O}</code>。PLAN_ZH.txt、PROTOCOL_V1.json、receipts/、logs/、tables/、objects/、FINAL_ACCEPTANCE.json均保留；下载目录仅交付报告。科学图件继续暂停。</p>'''
report='<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>CAF共同重分层与原亚群对应</title><style>'+css+'</style><main>'+body+'</main></html>';rp=O/'delivery/CAF_joint_restratification_report.html';rp.write_text(report,encoding='utf-8')
frames={'原33亚群对应':routes,'共同新群来源和稳定性':clustersummary,'全部患者级功能统计':primary,'留出患者功能统计':hold,'原亚群细胞去向':read('old33_by_joint_counts.tsv'),'原亚群患者等权去向':read('old33_patient_balanced_fraction.tsv')}
payload={k:json.loads(v.to_json(orient='records',double_precision=15)) for k,v in frames.items()}
js='''const datasets=PAYLOAD;const sel=document.getElementById('set'),q=document.getElementById('q');let off=0;for(const k of Object.keys(datasets)){let o=document.createElement('option');o.textContent=k;sel.append(o)}function draw(){let data=datasets[sel.value],terms=q.value.toLowerCase().trim().split(/\\s+/).filter(Boolean),rs=data.filter(r=>terms.every(t=>JSON.stringify(r).toLowerCase().includes(t))),cols=Object.keys(data[0]||{});document.getElementById('count').textContent=`匹配${rs.length}行；显示${off+1}至${Math.min(off+60,rs.length)}`;let t=document.getElementById('table');t.replaceChildren();let tr=document.createElement('tr');for(let c of cols){let th=document.createElement('th');th.textContent=c;tr.append(th)}t.append(tr);for(let r of rs.slice(off,off+60)){tr=document.createElement('tr');for(let c of cols){let td=document.createElement('td');td.textContent=r[c]===null?'不可评估':r[c];tr.append(td)}t.append(tr)}return rs.length}sel.onchange=q.oninput=()=>{off=0;draw()};document.getElementById('next').onclick=()=>{if(off+60<draw())off+=60;draw()};document.getElementById('prev').onclick=()=>{off=Math.max(0,off-60);draw()};draw();'''.replace('PAYLOAD',json.dumps(payload,ensure_ascii=False,allow_nan=False))
ap=O/'delivery/CAF_joint_full_tables.html';ap.write_text('<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>CAF共同重分层完整附表</title><style>'+css+'</style><main><h1>CAF共同重分层：完整统计与对应表</h1><p><a href="CAF_joint_restratification_report.html">返回主报告</a>。all为全部候选，singlet为非软件双细胞阳性，singlet_min5为患者/部位每侧至少5细胞。N为可评估患者；Δ为log1pCP10k评分差，P与FDR分别保留。不同来源不合并为独立N。原亚群对应比例不是身份概率。</p><select id="set"></select><p><input id="q" placeholder="搜索，如 NCI J00 collagen"><button id="prev">上一页</button><button id="next">下一页</button></p><p id="count"></p><div class="scroll"><table id="table"></table></div><script>'+js+'</script></main></html>',encoding='utf-8')
for f in [rp,ap]:assert '\ufffd' not in f.read_text(encoding='utf-8')
archive=O/'delivery/CAF_joint_restratification_Phase117.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in [rp,ap]:z.write(f,f.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
dump('FINAL_STATUS.json',dict(status='analysis_complete_local_delivery_pending',phase=O.name,completed_at=now(),cells=10510,fit_cells=10274,clusters=structure['clusters'],original33_routes=route_counts,F2C3_leading_cluster=focus_cluster,F2C3_leading_fraction=float(c3.leading_fraction),server_archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,report=str(rp),figure_queue_status='paused_by_user',remaining_required_work=['local delivery, render and receipt']))
event('06_acceptance_report','completed',clusters=structure['clusters'],archive_sha256=sha(archive))
