from pathlib import Path
import json,html,hashlib,zipfile,datetime,importlib.metadata
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase115_incremental_function_20261023';D=O/'delivery'
def load(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
stats=load('functional_statistics.tsv');meta=load('cell_mapping.tsv.gz');cl=load('query_cluster_membership.tsv.gz');stability=load('query_cluster_stability.tsv');diag=load('heldout_label_diagnostics.tsv');coverage=load('program_coverage.tsv');qa=json.loads((O/'ACCEPTANCE.json').read_text());assert qa['status']=='PASS'
main=stats[stats.cohort.eq('added')&stats.arm.eq('singlet')];f2=main[main.level.eq('class')&main.label.eq('F2')];extra=stats[stats.cohort.eq('added')&stats.arm.eq('singlet_min5')&stats.level.eq('class')&stats.label.eq('F2')]
endpoints=['fibrillar_I_III_V','collagen_crosslinking','myCAF_customer','F2C3_collagen_support','Osteogenic_customer']
cols=['source','label','endpoint','cells','N','mean_effect','positive_patients','P','FDR']
names={'fibrillar_I_III_V':'I/III/V型纤维胶原','collagen_crosslinking':'胶原交联','myCAF_customer':'原六基因myCAF收缩面板','F2C3_collagen_support':'F2-C3五基因支持面板','Osteogenic_customer':'成骨面板','inflammation':'炎症程序','iCAF_customer':'原iCAF面板','ECM_degradation':'ECM降解','mineralization':'矿化程序'}
rename={'source':'来源','arm':'分析层','label':'标签/计算群','endpoint':'程序或基因','cells':'目标细胞数','N':'可评估患者N','mean_effect':'平均效应Δ','positive_patients':'同向升高患者数','P':'原始P','FDR':'BH FDR','level':'层级','source_patients':'来源内患者记录'}
def table(df,columns=None,translate=True):
 q=df[columns].copy() if columns else df.copy()
 if translate:
  if 'endpoint' in q:q['endpoint']=q.endpoint.map(lambda x:names.get(x,x))
  q=q.rename(columns=rename)
 return '<div class="scroll">'+q.to_html(index=False,escape=True,border=0,float_format=lambda x:f'{x:.6g}',na_rep='不可评估')+'</div>'
counts=meta[meta.cohort.eq('added')].groupby(['source','class_compatible']).agg(候选细胞=('cell_id','size'),来源内患者记录=('donor','nunique'),软件双细胞=('doublet','sum')).reset_index()
class_summary=pd.DataFrame([
 ['F1','炎症/外膜样','NCI炎症程序在13/16位患者升高；NCC未复现同类名义信号。','部分功能支持；标签迁移表现偏低。'],
 ['F2','收缩—活化基质','NCC与NCI均见纤维胶原、交联与收缩相关增强，最低细胞数敏感性保留方向。','本轮主要功能增量；保留F2相容候选定位。'],
 ['F3','基质重塑','独立校准未接收F3，无法建立可靠新增F3组作比较。','不可评估；不等于新队列缺乏基质重塑。'],
 ['F4','成骨/基质重叠','NCI F4相容组的成骨面板当前未达名义显著；NCC可评估患者较少。','未形成一致的新功能支持。'],
 ['F5','混合/表达未定','NCI保留收缩/成骨等混合名义信号，炎症程序反向。','不升级为单一功能CAF亚型。']],columns=['原类别','原功能定位','新增证据','增量判断'])
focus=main[main.level.eq('cluster')&main.label.isin(['NCI_Q10','NCI_Q0','NCC_Q4','NCI_Q1'])&main.endpoint.isin(endpoints)]
cluster_sensitive=stats[stats.level.eq('cluster')&stats.label.eq('NCI_Q10')&stats.arm.eq('singlet_min5')&stats.endpoint.isin(endpoints)]
state=diag[diag.label.eq('F2_C3')].iloc[0];newstats=stats[stats.cohort.eq('added')]
css='''body{font:16px/1.7 "Microsoft YaHei","Noto Sans CJK SC",sans-serif;color:#20313c;background:#f4f7f9;margin:0}main{max-width:1180px;margin:30px auto;padding:38px 44px;background:white}h1{font-size:29px;color:#173e52}h2{margin-top:35px;border-bottom:2px solid #dbe5eb;padding-bottom:8px;color:#204e61}h3{font-size:18px}p{margin:12px 0}.lead{background:#eaf4f5;border-left:5px solid #257c85;padding:17px}.note{color:#566874;font-size:14px}.scroll{overflow:auto;margin:20px 0}table{border-collapse:collapse;font-size:13px;line-height:1.5;width:100%}th{background:#234d60;color:white;text-align:left;padding:9px;white-space:nowrap}td{border-bottom:1px solid #dae3e7;padding:8px;vertical-align:top}tr:nth-child(even){background:#f3f7f9}code{overflow-wrap:anywhere}a{color:#136d8b}.negative{background:#fff6e8;padding:14px}li{margin:6px 0}@media print{body{background:white}main{margin:0;padding:0;max-width:none}table{font-size:10px}h2{break-after:avoid}.scroll{overflow:visible}}'''
body=f'''<h1>CAF F1–F5 功能分类：新增队列增量评估</h1>
<p class="note">内部分析报告｜Phase115｜2026-10-23｜图件队列已按用户要求暂停；本轮未生成科研图件。</p>
<div class="lead"><b>主要结论：</b>新增NCC与NCI数据对“胶原／基质增强并伴收缩特征”的功能状态提供了增量支持，最集中于探索性F2相容组。当前支持功能层面的延伸，尚未确认稳定的新F2-C3身份或专属分选亚群。来源内分群也发现胶原增强候选，但小群的边界和患者覆盖仍限制其亚群命名。</div>
<h2>1. 数据与既有基准</h2>
<p>保留原8,174个细胞、F1–F5和33个内部状态。新增NCC 806个细胞，使用作者修订后的7位来源内患者记录；NCI 1,530个作者注释CAF细胞，按已接受的分组对应17位来源内患者记录。合计10,510个细胞，不将来源内患者记录直接相加为已验证的跨研究独立样本量。NCI排除软件双细胞后的主要敏感性层为1,320个细胞；NCC 806个均为当前单细胞层。</p>
<p>本次不是骨／肺差异分析：NCC实际取材部位仍未知，NCI仅有肺病灶。原F2-C3历史校正计数结果与本轮共享基因、原始计数结果分别保留，不将新统计覆盖旧报告。</p>
{table(class_summary,translate=False)}
<h2>2. 新增数据支持F2相容候选中的胶原与收缩程序增强</h2>
<p>NCC的F2相容候选为354个细胞，在7位可评估患者中，纤维胶原、胶原交联、myCAF六基因及F2-C3五基因支持面板均7/7同向升高，平均评分差分别为+0.532、+0.238、+0.425和+0.480。小样本定位为探索性；统计表保留名义显著及未通过FDR的原始数值。</p>
<p>NCI的F2相容候选为225个细胞，纤维胶原和交联程序均在11/11位患者升高，平均效应分别+0.848和+0.328，原始P均0.000977，FDR均0.009766；myCAF面板为9/11同向，效应+0.356，P=0.003906、FDR=0.023872；F2-C3五基因支持面板10/11同向，效应+0.523，P=0.004883、FDR=0.028269。这些统计描述探索性分组内的功能差异。</p>
{table(f2[f2.endpoint.isin(endpoints)],cols)}
<p class="note">Δ单位为患者内“目标细胞平均log1p(CP10k)－其他候选细胞平均log1p(CP10k)”的患者等权均值，不是倍数变化；每个程序为其共同可测成员的均值。myCAF面板、胶原面板和F2-C3面板并非互相独立验证。</p>
<h3>提高最低细胞数后方向保留</h3>
<p>每位患者／部位的目标与对照均至少5个单细胞时，NCC保留6位可评估患者，四个主要程序均6/6升高；NCI保留7位，纤维胶原、交联和F2-C3支持面板均7/7升高，myCAF为6/7。NCI该层相关FDR不再小于0.05，保留名义显著和效应方向。</p>
{table(extra[extra.endpoint.isin(endpoints)],cols)}
<h2>3. F2-C3的完整基因组合未被原样复制</h2>
<p>原报告F2-C3的胶原信号仍保留；本轮显示功能相似，并非五个支持基因逐一完整复现。NCC与NCI主要F2相容组的PCOLCE均未达名义显著，LOXL2也未达名义显著。NCI在最低5细胞敏感性层PCOLCE出现名义信号，而LOXL2仍不支持。不能将多基因均值升高表述为每个基因均升高。</p>
{table(f2[f2.endpoint.isin(['gene__COL3A1','gene__COL5A1','gene__COL5A2','gene__PCOLCE','gene__LOXL2'])],cols)}
<p class="negative">主要并行结果：NCI的F2相容组还显示成骨面板11/11同向升高（Δ=+0.217，P=0.000977，FDR=0.009766），NCC没有同样结果。新增NCI候选应保留“基质与成骨程序重叠”的解释，不据此声称纯myCAF或已确认基质沉积功能。</p>
<h2>4. 不依赖旧标签的来源内分群检出功能候选，但未形成稳定对应亚群</h2>
<p>排除待检验功能及谱系面板基因后，NCC单细胞分为7个计算群，NCI单细胞分为12个计算群；未进行跨来源共同空间整合或重新定义旧标签。主要分群分辨率固定为0.5，另记录不同种子及0.8分辨率敏感性。</p>
<p>NCI_Q10包含30个单细胞，来自6位患者。纤维胶原、交联及F2-C3支持面板均6/6同向升高，平均效应分别+1.582、+0.471、+1.183，为名义显著的探索性基质增强候选。该群不等同于旧F2-C3。</p>
{table(focus,cols)}
<p>Q10的关键敏感性：换随机种子后最佳匹配Jaccard仅0.129，说明其作为独立群的边界不稳；最低5细胞条件只保留2位患者，方向仍正但不能维持统计支持。NCI_Q0的myCAF面板6/6同向升高，但纤维胶原并不相应增强，提示收缩型与基质型应分别评价。NCC_Q4的纤维胶原和F2-C3支持面板4/4同向，但当前检验未达名义显著；其余NCC群存在明显患者主导现象。</p>
{table(cluster_sensitive,cols)}
<details><summary>查看全部来源内计算群稳定性与患者覆盖</summary>{table(stability,translate=False)}</details>
<h2>5. 旧F类标签迁移仍需保留探索性定位</h2>
<p>参考数据按来源分层、按患者划分训练／校准／测试。功能与目标基因未用于映射特征选择。仅凭这一隔离并不能保证生物学身份可靠：独立留出患者中，F2接收标签与原标签一致29/77；在大类一致性门槛也满足后，F2-C3为{int(state.correct_accepted)}/{int(state.accepted_cells)}，样本很少。F3无通过校准的接收标签，新增F3数量不能解释为生物学零。</p>
{table(diag[diag.level.eq('class')|diag.label.eq('F2_C3')],translate=False)}
<p>因此，本轮采用“F2相容候选中的功能增强”“来源内基质增强计算群”等表述。不能把本轮接收数量登记成新增已确认F2-C3患者、或用F4/F5局部较高一致率宣称五类身份整体已验证。模拟遗漏某一参考大类时仍存在较高误接收，未定类识别能力有限。</p>
<details><summary>查看新增细胞的探索性大类映射数量</summary>{table(counts,translate=False)}</details>
<h2>6. Methods</h2>
<p>计数按完整文库归一化至每万计数后log1p转换。基因取所有分析对象中唯一且无歧义符号的共同集合，原队列同时要求所有样本具有测量覆盖。沿用原11个程序、3个胶原结构面板、8个客户程序，并加入由既有F2-C3证据定义的5基因支持面板，共23个程序。共同覆盖不足80%的程序不作推断检验；本轮apCAF仅CD74进入共有空间，保留覆盖信息并不进行apCAF程序检验。</p>
<p>映射模板先按来源内患者／标签求均值，再按来源等权；基因标准化后以模板相关和3个不同患者近邻投票共同判定。训练、校准、测试按固定种子对应的患者哈希次序分离。接收阈值来自校准集中双方法正确记录的相似度和差值第5百分位，测试集不参与阈值拟合。细状态须与独立大类模型父类一致。全参考最近模板仅作附加描述性敏感性，不代替主判定。</p>
<p>来源内分群仅以非双细胞拟合；排除功能、谱系及目标基因和技术基因后选2000个HVG，标准化、PCA30、15近邻及Leiden。双细胞记录保留，并按最近单细胞群质心分配，供全细胞敏感性展示。未生成UMAP或科研图件。</p>
<p>功能检验按来源分别实施：同一患者、同一部位内目标组对其余候选细胞，先计算程序平均分差，再对该患者各部位等权，最后以患者为单位作双侧符号翻转检验。N≤16为精确枚举，其余为65536次蒙特卡洛加观测值。N&lt;2、全零差值或覆盖不足时P缺失并标记原因。主比较与全部保留、排除双细胞、最低5细胞三层并行；程序与目标单基因分别以来源×分析层×分组层级定义BH校正家族。原始P&lt;0.05为名义显著，FDR&lt;0.05才称校正后显著。</p>
<h2>7. Discussion与后续定位</h2>
<p><b>本轮可推进的对象是跨患者重复出现的基质增强功能状态。</b>旧F2功能主线获得新增队列支持；精确F2-C3身份与稳定分群边界仍未闭合。优先将F2相容的胶原增强细胞作为后续验证候选，同时保留收缩、成骨重叠及来源信息；不将NCI_Q10直接命名为新的F2-C3或新CAF亚型。</p>
<p>功能评分来自RNA，不能证明实际胶原蛋白沉积、纤维形成、促转移或可用流式分选。原队列包含混合谱系候选，NCI沿用作者CAF注释且平台为单核RNA；映射和来源内聚类均为同一转录组上的计算分组，排除评分基因减少直接循环但不等于独立实验验证。NCC取材信息、跨研究患者去重及NCI原始逐文库特征掩码仍沿用既有边界。原始P与FDR、反向结果和不可评估记录完整保留。</p>
<h2>8. 交付与验收</h2>
<p>独立复核：{qa['independent_exact_or_exact_MC_comparisons']}项新增主要／最低细胞数层检验与独立枚举结果核对；全部BH家族重算通过；原8,174细胞标签及权威哈希保持，7份NCC计数对象与归档哈希一致。所有原生任务退出码为0。通过表示本轮计算与交付范围完成，不代表CAF身份或纤维沉积机制已验证。</p>
<p><a href="CAF_incremental_statistics.html">打开完整新增统计附表（含全部方向、三种敏感性层、P/FDR与不可评估记录）</a></p>
<p class="note">服务器证据目录：<code>{O}</code>。主要表：functional_statistics.tsv、patient_function_effects.tsv.gz、heldout_label_diagnostics.tsv、query_cluster_stability.tsv、program_coverage.tsv。基线、代码、细胞级数据、日志与验收证据留在服务器。此为内部分析交付，未生成客户改写版。</p>'''
report='<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CAF增量功能评估 Phase115</title><style>'+css+'</style><main>'+body+'</main></html>'
(D/'CAF_incremental_function_report.html').write_text(report,encoding='utf-8')
# Full new-cohort results, searchable with bounded rendering, preserving numeric source precision in JSON.
records=json.loads(newstats.to_json(orient='records',double_precision=15));headers=newstats.columns.tolist()
script='''const rows=DATA, cols=COLS;let offset=0,filtered=rows;const input=document.getElementById('q');function draw(){document.getElementById('count').textContent=`匹配 ${filtered.length} 条；当前 ${offset+1}–${Math.min(offset+100,filtered.length)}`;let t=document.getElementById('data');t.replaceChildren();let h=document.createElement('tr');for(const c of cols){let th=document.createElement('th');th.textContent=c;h.appendChild(th)}t.appendChild(h);for(const r of filtered.slice(offset,offset+100)){let tr=document.createElement('tr');for(const c of cols){let td=document.createElement('td');td.textContent=r[c]===null?'不可评估':r[c];tr.appendChild(td)}t.appendChild(tr)}}input.addEventListener('input',()=>{const terms=input.value.toLowerCase().trim().split(/\\s+/).filter(Boolean);filtered=rows.filter(r=>terms.every(v=>JSON.stringify(r).toLowerCase().includes(v)));offset=0;draw()});document.getElementById('next').onclick=()=>{if(offset+100<filtered.length)offset+=100;draw()};document.getElementById('prev').onclick=()=>{offset=Math.max(0,offset-100);draw()};draw();'''.replace('DATA',json.dumps(records,ensure_ascii=False,allow_nan=False)).replace('COLS',json.dumps(headers))
appendix='<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>CAF完整新增统计表</title><style>'+css+'input{width:75%;padding:10px}button{padding:8px;margin:8px}</style><main><h1>CAF新增队列完整统计附表</h1><p><a href="CAF_incremental_function_report.html">返回主报告</a>｜all=全部候选；singlet=排除双细胞；singlet_min5=目标和对照每患者/部位至少5细胞。class/state为探索性标签，cluster为来源内计算群；N为可评估患者数，cells不是推断N。Δ为平均log1pCP10k差。P、FDR与检验家族按原结果完整保留。</p><p>source×arm×level×family_kind构成FDR家族；低于10例按探索性解读。搜索支持空格分隔多个条件，例如 NCI singlet F2 collagen。</p><input id="q" placeholder="搜索来源、群、程序或分析层"><button id="prev">上一页</button><button id="next">下一页</button><p id="count"></p><div class="scroll"><table id="data"></table></div><script>'+script+'</script></main></html>'
(D/'CAF_incremental_statistics.html').write_text(appendix,encoding='utf-8')
environment={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','scanpy','anndata','scikit-learn','leidenalg']};(O/'resources/SOFTWARE_VERSIONS.json').write_text(json.dumps(environment,indent=2),encoding='utf-8')
for f in D.glob('*.html'):
 t=f.read_text(encoding='utf-8');assert '\ufffd' not in t and 'charset="UTF-8"' in t
archive=D/'CAF_incremental_function_Phase115.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(D.glob('*.html')):z.write(f,f.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
status=dict(status='computed_verified_report_pending_local_render',phase=O.name,completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),cells=10510,added_cells=2336,added_sources=['NCC','NCI'],NCC_source_patients=7,NCI_source_patients=17,figure_queue_status='paused_by_user',original_baseline_unchanged=True,analysis_acceptance='PASS_EXPLORATORY_FUNCTIONAL_ASSESSMENT',server_report=str(D/'CAF_incremental_function_report.html'),server_archive=str(archive),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_bytes=archive.stat().st_size,appendix_rows=len(newstats),all_figure_adaptation_complete=False,remaining_required_work=['local delivery, browser rendering and readback'])
(O/'FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(status))
