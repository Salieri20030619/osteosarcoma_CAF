from pathlib import Path
import pandas as pd,json,html,zipfile,hashlib,datetime
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase116_original33_analogs_20261023';D=O/'delivery'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
qa=json.loads((O/'ACCEPTANCE.json').read_text());assert qa['status'].startswith('PASS')
d=read('all_33_by_19_comparisons.tsv');old=read('original_33_subgroup_dictionary.tsv').fillna('');cov=read('original_marker_coverage.tsv');pathcov=read('inherited_pathway_coverage.tsv');proto=json.loads((O/'resources/PROTOCOL.json').read_text());genes=read('original_marker_gene_evidence.tsv');ps=read('inherited_pathway_newgroup_statistics.tsv');prior=read('inherited_pathway_display.tsv')
c3=d[d.old_state.eq('F2_C3')&d.query_group.eq('NCI_Q10')].iloc[0];c4=d[d.old_state.eq('F1_originalC4')&d.query_group.eq('NCI_Q2')].iloc[0]
css='body{font:16px/1.7 "Microsoft YaHei",sans-serif;margin:0;background:#f3f6f8;color:#203440}main{max-width:1200px;margin:28px auto;padding:36px;background:white}h1{font-size:29px}h2{margin-top:34px;border-bottom:2px solid #dbe5e9;color:#245364}h3{font-size:18px}.lead{background:#eaf4f4;padding:18px;border-left:5px solid #267c86}.note{font-size:14px;color:#536673}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:16px 0}th{background:#254e60;color:white;padding:9px;text-align:left}td{padding:8px;border-bottom:1px solid #dce4e8;vertical-align:top}tr:nth-child(even){background:#f3f7f9}a{color:#136b88}code{overflow-wrap:anywhere}details{margin:14px 0}input{padding:9px;width:70%}button{padding:9px;margin-left:6px}'
def tab(x):return '<div class="scroll">'+x.to_html(index=False,border=0,escape=True,na_rep='不可评估',float_format=lambda v:f'{v:.5g}')+'</div>'
def text(s):return html.escape(str(s))
summary=[]
for _,r in old.iterrows():
 row={'原亚群':r.subgroup_id,'原功能名称':r.working_name,'原代表标记':r.representative_markers,'原功能解释置信度':r.functional_interpretation_confidence}
 for src in ['NCC','NCI']:
  z=d[(d.old_state==r.subgroup_id)&(d.source==src)].sort_values('marker_rank_for_old_state').iloc[0]
  row[src+'标记最近候选']=f'{z.query_group}；标记r={z.marker_similarity:.3f}，功能r={z.function_similarity:.3f}；N={int(z.N)}；一致标记{int(z.consistent_positive_marker_genes)}/{int(z.signature_genes)}'
 summary.append(row)
summary=pd.DataFrame(summary);summary.to_csv(O/'tables/original33_correspondence_summary_ZH.tsv',sep='\t',index=False,encoding='utf-8')
focusstates=['F1_originalC10','F2_C0','F2_C1','F2_C3','F3_C0','F3_C1','F3_C2','F4_C1','F4_C5','F4_C8','F4_C9']
focus=d[d.old_state.isin(focusstates)&d.marker_rank_for_old_state.eq(1)]
displaycols=['source','old_state','query_group','marker_similarity','function_similarity','N','signature_effect','P','FDR','min5_N','seed_jaccard']
g=genes[(genes.old_state=='F2_C3')&(genes.query_group=='NCI_Q10')][['gene','mean_effect','N','positive_patients']]
historical=prior[(prior.contrast=='F2_C3')&(prior.scope=='all')&(prior.NES>0)&(prior.GSEA_model_P<.05)][['library','Term','NES']].drop_duplicates(['library','Term']);historical['pathway']=historical.library+'::'+historical.Term
repro=historical.merge(ps[(ps.query_group=='NCI_Q10')&(ps.minimum==1)],on='pathway',how='left');repro=repro[['Term','NES','N','mean_effect','positive_patients','P','FDR']]
specific=d[((d.old_state=='F2_C3')&(d.query_group=='NCI_Q10'))|((d.old_state=='F1_originalC4')&(d.query_group=='NCI_Q2'))|((d.old_state=='F1_originalC10')&(d.query_group=='NCC_Q1'))|((d.old_state=='F4_C9')&(d.query_group=='NCC_Q4'))]
body=f'''<h1>原33个CAF内部亚群：新增队列相似群检索</h1><p class="note">内部分析｜Phase116｜2026-10-23｜按用户纠正的范围开展：以原亚群画像寻找相似群体。图件队列保持暂停。</p>
<div class="lead"><b>本轮参照单位为原33个内部亚群，而非F1–F5五个父类。</b>直接使用Phase62已建立的标记组合、名称与通路基因集，逐一比较NCC、NCI候选群，并保留细胞级33向相似度检索。没有先通过五类标签的门槛，没有将所有基质增强细胞统称F2。F2-C3在NCI中出现部分相似候选NCI_Q10；其他原亚群亦有可追踪对象，但不强行建立一对一身份对应。</div>
<h2>1. 原分类的33个内部亚群完整保留</h2>
<p>原F1下3个、F2下6个、F3下3个、F4下17个、F5下4个亚群，合计33个。原8,174个细胞的标签与命名依据保留。本轮新增NCC 806个细胞、7位来源内患者记录，NCI 1,530个细胞、17位来源内患者记录；签名检验以排除软件双细胞后的NCC 806个和NCI 1,320个细胞为主。</p>
<p>新队列7个NCC计算群与12个NCI计算群仅作为检索容器：对每个新群比较所有33个原亚群，共627个组合；另对全部2,336个新增细胞直接检索33个原亚群，不受上一轮大类接收结果限制。新群ID沿用Phase115以保持可追溯，不代表已确认的新生物学亚型。</p>
<h2>2. F2-C3在NCI中存在部分标记与功能相似候选</h2>
<p>原F2-C3为“CD248/POSTN胶原基质相关态”，原稳定标记为CD248、MMP9、OLFML3、VCAN、CFH、POSTN、P4HA2、COL16A1、PRSS23、COL5A3、FKBP11、AP001107.3。共有测量空间覆盖其中11个，缺少AP001107.3。这里使用这组原始标记，而非将五基因胶原支持面板当作原亚群定义。</p>
<p>NCI_Q10是NCI各计算群中对原F2-C3标记画像最接近的候选：标记相似度r={c3.marker_similarity:.3f}，在该新群对33个旧亚群的排名为第{int(c3.marker_rank_within_query)}；完整功能画像r={c3.function_similarity:.3f}，排名第{int(c3.function_rank_within_query)}。该群30个单细胞来自6位患者，原标记组合的患者内平均评分差为+{c3.signature_effect:.3f}，构成探索性名义信号，完整P/FDR见下表。</p>
<p>具体原标记中，POSTN、VCAN在5/6位患者同向升高，PRSS23为4/6；CD248为3/6，MMP9仅1/6且平均方向不升。故本轮定位为“F2-C3部分相似候选”，不能将其改名为已确认F2-C3。该新群对其他F2内部状态也有相似性，保留竞争对应。</p>
{tab(specific[displaycols])}
{tab(g)}
<p class="note">r是标准化表达／程序画像的相关相似度，不是身份概率。mean_effect/签名Δ为患者内目标群与其余候选细胞的平均log1p(CP10k)分差，不是fold change。P是当前计算分组条件下的患者级签名检验，原始P&lt;0.05为名义显著，FDR另列。</p>
<h3>原F2-C3通路画像在该候选群中的方向</h3>
<p>以下从原F2-C3归档中提取已有正向通路，使用原完整基因成员在新数据重算分数。原NES为历史基因集置换模型结果；新增P/FDR为患者级平均评分差检验，二者统计量与检验单位不同，不作数值等同。</p>
{tab(repro)}
<p>主要敏感性：NCI_Q10更换随机种子后最佳群匹配Jaccard=0.129；每位患者目标与其余组均至少5细胞时，仅2位可评估。排除参与新群聚类的标记后，F2-C3仍有{int(c3.unused_genes)}个可用标记，效应方向为{c3.unused_effect:+.3f}，P={c3.unused_P:.4g}。这些结果支持保留候选，不支持稳定边界或全部原标记复现。</p>
<h2>3. 同时检索其他原基质／收缩亚群</h2>
<p>F2-C3并非唯一检索对象。F1_originalC10（SFRP2/PLA2G2A基质与分泌相关态）与NCC_Q1在标记和功能上均靠前，但该群由单一患者高度主导，只有2位患者参与比较。F4_C9（TNC/COL3A1基质相关态）在NCC中最接近NCC_Q4，但原标记共同覆盖不足80%，不把部分面板检验作完整签名显著性。F3_C1、F3_C2等也保留候选排序及实际患者覆盖，不因未通过大类映射而从检索中删除。</p>
{tab(focus[displaycols])}
<h2>4. 非基质亚群亦保留逐项对应</h2>
<p>原F1_originalC4（LARGE1/CD44应答相关态）与NCI_Q2具可追踪的对应线索：12位可评估患者，{int(c4.consistent_positive_marker_genes)}个原标记在多数患者中同向，原签名Δ={c4.signature_effect:+.3f}，P={c4.P:.6g}，FDR={c4.FDR:.6g}。最低5细胞层仍有{int(c4.min5_N)}位患者，签名Δ={c4.min5_effect:+.3f}、P={c4.min5_P:.6g}、FDR={c4.min5_FDR:.6g}。这与“每个父类重新找差异”的问题不同，是从特定旧亚群画像出发的候选对应。</p>
<h2>5. 全部33个原亚群的逐项索引</h2>
<p>下表逐项列出每个原亚群在两个新来源中标记画像最近的计算群。<b>最近者不等于确认对应</b>：负相似、功能不一致、少患者和低覆盖均保留。完整附表给出每个旧亚群的前三候选，以及627个组合的全部证据，可检查一对多、多个旧亚群共享同一新群等情况。</p>
{tab(summary)}
<h2>6. Methods</h2>
<p>原参照为Phase62的F1F5_subgroup_summary.tsv与marker_raw_corrected_evidence.tsv。优先采用原stable_markers，无稳定列表时才使用原display_markers并记录。参考表达直接来自原8,174细胞的共享测量基因空间，模板排除26个已有软件双细胞阳性，以部位、患者和来源分级等权求均值；新群以患者等权。联合原标记空间标准化后计算33向相关相似度，另给原始表达Spearman排序作敏感性。</p>
<p>功能画像包括22个覆盖充分的既有程序，以及原Phase62报告归档中{int(((pathcov.coverage>=.8)&(pathcov.measured>=3)).sum())}个覆盖充分的显示通路，均沿用原基因成员重算细胞分数，不仅复用旧通路标题。历史来源、原功能名称、功能解释置信度及原通路证据另外归档。亚群间共享或重叠的基因集不视为独立证据。</p>
<p>原标记签名在来源内、患者内比较目标新群与其他候选细胞，以均值分差为效应，精确双侧符号翻转计算P；不足2位患者、全零效应或基因覆盖不足时不补造P。原完整签名与“未参与新群聚类的标记”分别评价，最低细胞数1与5分别记录。BH家族为来源×签名类型×最低细胞数内全部33签名和新群组合；原归档通路的新增检验以来源×最低细胞数为独立通路家族。以原始P&lt;0.05保留名义候选，FDR不作为唯一候选门槛。</p>
<h2>7. 解释与使用范围</h2>
<p>本轮完成的是以原内部亚群为参照的增量相似群检索。原标记最初多来自同一父类内亚群对照，而新签名检验对照为来源内其他候选细胞，对照范围不同；新群形成也可能使用部分原标记。因此主签名检验属于计算分组条件下的探索性支持，未参与聚类标记的敏感性单列。保留历史标记后再作原队列留一患者检索只是内部再识别诊断，不是独立外部验证。</p>
<p>相似群不要求强行一对一；没有相似候选也不能解释为对应生物学状态不存在。原低跨患者支持、单患者主导、技术敏感或混合谱系状态继续保留原证据等级。RNA基质程序与相似画像不等于实际纤维沉积或已验证分选门。当前应优先追踪F2-C3—NCI_Q10等明确的旧亚群—新候选对应，而不是再次把全部基质信号合并为F2。</p>
<h2>8. 文件与验收</h2>
<p>全部33原亚群、19新群、627组合与2,336新增细胞完整保留；{qa['signature_tests_independently_recomputed']}项签名检验独立重算，最大P误差{qa['max_exact_P_error']}；签名和通路BH家族重算通过；15份输入哈希未变。原标签未修改，图件保持暂停。</p>
<p><a href="CAF_original33_full_comparisons.html">查看全部对应组合、原标记检验与历史通路的新数据结果</a></p><p class="note">服务器证据目录：<code>{O}</code>。Phase115保留为较粗层级分析，不能替代本轮逐亚群检索；本报告为内部分析，不是客户改写版。</p>'''
report='<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>原33亚群相似群检索</title><style>'+css+'</style><main>'+body+'</main></html>'
(D/'CAF_original33_analog_report.html').write_text(report,encoding='utf-8')
frames={'全部627个对应组合':d,'每原亚群前三候选':read('per_original_subgroup_top3.tsv'),'原签名全部检验':read('all33_signature_tests.tsv'),'原通路新增结果':ps,'原标记基因覆盖':cov,'原通路覆盖':pathcov,'参考亚群内部留患者检索':read('reference_retrieval_diagnostics.tsv')}
payload={k:json.loads(v.to_json(orient='records',double_precision=15)) for k,v in frames.items()}
js='''const datasets=PAYLOAD;const sel=document.getElementById('set'),input=document.getElementById('q');let off=0;for(const name of Object.keys(datasets)){let o=document.createElement('option');o.textContent=name;sel.appendChild(o)}function draw(){const data=datasets[sel.value],terms=input.value.toLowerCase().trim().split(/\\s+/).filter(Boolean),rs=data.filter(r=>terms.every(t=>JSON.stringify(r).toLowerCase().includes(t))),cols=Object.keys(data[0]);document.getElementById('count').textContent=`匹配 ${rs.length} 行，显示 ${off+1}–${Math.min(off+60,rs.length)}`;const table=document.getElementById('t');table.replaceChildren();let tr=document.createElement('tr');for(const k of cols){let th=document.createElement('th');th.textContent=k;tr.appendChild(th)}table.appendChild(tr);for(const r of rs.slice(off,off+60)){tr=document.createElement('tr');for(const k of cols){let td=document.createElement('td');td.textContent=r[k]===null?'不可评估':r[k];tr.appendChild(td)}table.appendChild(tr)}return rs.length}input.oninput=()=>{off=0;draw()};sel.onchange=()=>{off=0;draw()};document.getElementById('next').onclick=()=>{if(off+60<draw())off+=60;draw()};document.getElementById('prev').onclick=()=>{off=Math.max(0,off-60);draw()};draw();'''.replace('PAYLOAD',json.dumps(payload,ensure_ascii=False,allow_nan=False))
appendix='<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>原33亚群完整对应表</title><style>'+css+'</style><main><h1>原33亚群：完整检索与统计附表</h1><p><a href="CAF_original33_analog_report.html">返回主报告</a>。old_state为原内部亚群，query_group为新队列检索容器；不使用五类父群门槛。相似度不是身份概率，最近候选不等于确认对应。</p><p>all_original_markers为原标记完整可测部分，unused_in_query_clustering为未进入新群聚类的标记；minimum为患者内每侧最低细胞数。N为可评估患者数，P和FDR不可评估时留空值。retrieval_tier仅供导航，排名阈值不是生物学排除门槛。</p><select id="set"></select><p><input id="q" placeholder="搜索，如 F2_C3 NCI"><button id="prev">上一页</button><button id="next">下一页</button></p><p id="count"></p><div class="scroll"><table id="t"></table></div><script>'+js+'</script></main></html>'
(D/'CAF_original33_full_comparisons.html').write_text(appendix,encoding='utf-8')
for f in D.glob('*.html'):assert '\ufffd' not in f.read_text(encoding='utf-8')
archive=D/'CAF_original33_analogs_Phase116.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(D.glob('*.html')):z.write(f,f.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
status=dict(status='analysis_verified_local_render_pending',phase=O.name,updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),reference_subgroups=33,query_groups=19,pairs=627,query_cells=2336,scope='Original subgroup analog retrieval, no parent gate',server_archive=str(archive),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_bytes=archive.stat().st_size,figure_queue_status='paused_by_user',original_results_modified=False,remaining_required_work=['local delivery and render'])
(O/'FINAL_STATUS.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status))
