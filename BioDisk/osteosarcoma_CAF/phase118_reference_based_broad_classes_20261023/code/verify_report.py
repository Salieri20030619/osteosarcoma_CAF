from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import numpy as np,pandas as pd,json,hashlib,itertools,datetime,html,zipfile
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase118_reference_based_broad_classes_20261023';B=R/'phase117_joint_restratification_20261023'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8')
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
ann=read('J_to_broad_class_FROZEN.tsv');m=read('cell_broad_class_registry.tsv.gz');summ=read('broad_class_summary.tsv');stats=read('broad_class_statistics.tsv');effects=read('broad_patient_effects.tsv.gz');panels=json.loads((O/'resources/PANELS_FROZEN.json').read_text());frozen=json.loads((O/'ANNOTATION_FROZEN.json').read_text())
assert sha(O/'tables/J_to_broad_class_FROZEN.tsv')==frozen['sha256'];assert len(m)==10510 and m.cell_id.is_unique and (~m.doublet).sum()==10274;assert len(ann)==22 and ann.joint_cluster.is_unique;assert summ.all_cells.sum()==10510 and summ.main_cells.sum()==10274
original=pd.read_csv(B/'tables/frozen_registry.tsv.gz',sep='\t',low_memory=False);assert m.cell_id.equals(original.cell_id)
for c in ['original_class','original_state','doublet']:assert m[c].fillna('NA').equals(original[c].fillna('NA'))
raw=np.load(O/'objects/target_gene_values.npz');gs=raw['genes'].tolist();gi={g:i for i,g in enumerate(gs)};max_score_error=0.
for mr in [1000,1500,2000]:
 s=read('cell_rank_scores_'+str(mr)+'.tsv.gz')
 for name,genes in panels.items():
  rr=raw['ranks'][:,[gi[g] for g in genes]];n=np.isfinite(rr).sum(1);maximum=n*mr-n*(n+1)/2;total=np.nansum(np.minimum(rr,mr),axis=1);expected=np.divide(n*mr-total,maximum,out=np.full(len(m),np.nan),where=maximum>0).clip(0,1);expected[n<max(3,int(np.ceil(.8*len(genes))))]=np.nan
  assert np.array_equal(np.isnan(expected),s[name].isna().to_numpy());er=float(np.nanmax(abs(expected-s[name].to_numpy()))) if np.isfinite(expected).any() else 0.;max_score_error=max(max_score_error,er);assert er<1e-12
maxbh=0.
for _,g in stats.groupby(['source','arm','family']):
 g=g[g.P.notna()].sort_values('P');n=len(g)
 if not n:continue
 expected=np.minimum.accumulate((g.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1);er=float(np.max(abs(expected-g.FDR.to_numpy())));assert er<1e-10;maxbh=max(maxbh,er)
assert stats.loc[stats.N.lt(2),'P'].isna().all();assert stats.loc[stats.test.eq('all_zero_not_tested'),'P'].isna().all()
subset=stats[stats.source.isin(['NCC','NCI'])&stats.arm.isin(['singlet','singlet_min5'])&stats.family.eq('signature')];cache={};checked=0;maxerr=0.;maxmc=0.
for key,g in subset.groupby(['source','arm','broad_class']):
 src,arm,cls=key;ee=effects[(effects.source==src)&(effects.arm==arm)&(effects.broad_class==cls)]
 for _,row in g.iterrows():
  v=ee.loc[ee.endpoint.eq(row.endpoint),'effect'].to_numpy();n=len(v);assert n==row.N
  if n:assert abs(v.mean()-row.mean_effect)<1e-10;assert int((v>1e-12).sum())==row.positive_patients
  if pd.isna(row.P):continue
  if n not in cache:cache[n]=np.array(list(itertools.product((-1.,1.),repeat=n)))
  p=float(np.mean(abs(cache[n]@v)>=abs(v.sum())-1e-10));er=abs(p-row.P)
  if row.test=='exact_signflip':assert er<1e-10;maxerr=max(maxerr,er)
  else:assert er<max(.001,6*np.sqrt(p*(1-p)/65536));maxmc=max(maxmc,er)
  checked+=1
hashes=json.loads((O/'resources/INPUT_HASHES.json').read_text())
for h in hashes:assert sha(h['path'])==h['sha256']
manifest=json.loads((O/'references/SOURCE_MANIFEST.json').read_text());assert sha(O/'references/user_four_class_report.docx')==manifest[0]['sha256'];assert sha(O/'references/user_CD36_paper.pdf')==manifest[1]['sha256']
assert (O/'logs/profile.exit_code').read_text().strip()=='0' and (O/'logs/annotate_analyze.exit_code').read_text().strip()=='0'
qa=dict(status='PASS_WORKING_ANNOTATION_AND_STATISTICS',all_cells=10510,main_cells=10274,J_groups=22,functional_classes=4,reserve_J_groups=int(ann.broad_class.eq('U').sum()),source_docs_hash_verified=True,original_F_J_unchanged=True,all_3_maxRank_scores_recomputed=True,max_score_error=max_score_error,independent_new_signature_tests=checked,max_exact_P_error=maxerr,max_MC_P_difference=maxmc,max_BH_error=maxbh,input_hashes_unchanged=len(hashes),annotation_frozen_before_inferential_tests=True,not_validated_CAF_identity=True)
dump('ACCEPTANCE.json',qa)
classes=json.loads((O/'resources/CLASS_DEFINITIONS.json').read_text());cols=['source','broad_class','endpoint','cells','N','mean_effect','positive_patients','P','FDR'];main=stats[stats.arm.eq('singlet')];key=main[main.family.eq('signature')&main.endpoint.isin(['doc_iCAF','doc_myCAF','doc_mCAF','doc_osteogenic','paper_CD36','paper_MYH11','paper_C7','paper_FNDC1'])]
profile=read('balanced_cluster_scores_min5.tsv');summary_table=ann[['joint_cluster','broad_class','name','secondary_state','annotation_support','reason_ZH']].merge(pd.read_csv(B/'tables/joint_cluster_origins_and_stability.tsv',sep='\t')[['joint_cluster','cells','NCC_cells','NCI_cells','dominant_source_fraction']],on='joint_cluster');save(summary_table,'final_J_broad_annotation_evidence.tsv')
byclass=[]
for _,z in summ.iterrows():byclass.append({'新标签':z.broad_class,'工作大类':z.broad_name,'J群':'; '.join(ann.loc[ann.broad_class.eq(z.broad_class),'joint_cluster']),'J群数':int(z.J_groups),'主分析细胞':int(z.main_cells),'全部保留细胞':int(z.all_cells)})
groups=pd.DataFrame(byclass);expr=read('source_marker_expression.tsv');lipid=expr[expr.source.eq('NCI')&expr.joint_cluster.isin(['J07','J10','J12','J20'])&expr.gene.isin(['CD36','FABP4','PTPRC','PECAM1','MYH11','RGS5'])]
css='body{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#213440;background:#f4f7f9;margin:0}main{max-width:1200px;margin:28px auto;padding:36px;background:white}h1{font-size:29px}h2{margin-top:32px;border-bottom:2px solid #dbe4e9;color:#245164}.note{color:#536975;font-size:14px}.lead{padding:16px;background:#eaf4f5;border-left:5px solid #287b86}.scroll{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%;margin:16px 0}th{background:#254e60;color:white;padding:9px;text-align:left}td{padding:8px;border-bottom:1px solid #dbe5e9;vertical-align:top}tr:nth-child(even){background:#f3f7f9}details{margin:14px 0}a{color:#126c89}code{overflow-wrap:anywhere}input{padding:9px;width:65%}button{padding:9px;margin:5px}'
def tab(d):return '<div class="scroll">'+d.to_html(index=False,border=0,escape=True,na_rep='不可评估',float_format=lambda v:f'{v:.6g}')+'</div>'
reference_table=pd.DataFrame([
 ['DOCX四分类报告','炎症/外膜样','APOD IGF1 ITM2A CXCL12 IL6','主轴L1；组织起源未证实'],
 ['DOCX四分类报告','收缩/活化基质','ACTA2 MYL9 TAGLN LRRC15','主轴L2；保留壁细胞重叠'],
 ['DOCX四分类报告','基质合成重塑','MMP11 POSTN COL14A1 COL11A1 LOX','主轴L3；区分强基质和应答过渡'],
 ['DOCX四分类报告','成骨/基质重叠','OMD OGN SOSTDC1 PHOSPHO1 ALPL','主轴L4；CAF身份待定'],
 ['CD36论文Fig1及第6页','C7型','C7 PTGDS TNXB FGL2 IGF2 CXCL12','炎症对照，不与DOCX完全同义'],
 ['CD36论文Fig1及第6页','MYH11型','MYH11 RGS5 MCAM MYH9 TGFBI','收缩/血管对照，单独检查壁细胞'],
 ['CD36论文Fig1及第6页','FNDC1型','FNDC1 COMP MMP2 COL11A1 COL10A1 P4HA3','基质对照'],
 ['CD36论文Fig1及第6页','CD36型','CD36 FABP4 APOE STEAP4 APOLD1','独立脂质轴，本轮不新增同名CAF大类']],columns=['材料','原定义','原核心标记','本轮使用方式'])
body=f'''<h1>基于原四分类报告与CD36论文的CAF大类重新归纳</h1><p class="note">内部工作注释｜Phase118｜2026-10-23｜沿用Phase117的22个J计算群，不重新聚类；原F系列、原33亚群及J标签保留。</p>
<div class="lead"><b>重新归纳为四个文献依据的功能候选方向L1–L4，另保留U混合/竞争/未定层。</b>L1对应J08；L2对应J07；L3对应J03、J19、J21；L4为10个成骨/基质重叠表达群。U有7个J群，不是第五种CAF类型。此为研究工作分类，并非四种已验证CAF生物学身份；L4尤其不能将全部成骨表达细胞直接视为CAF。</div>
<h2>1. 两份原始材料的四分类并不相同</h2><p>用户DOCX明确将四分类作为骨肉瘤观察与泛癌文献整合的工作框架，并非单篇骨肉瘤论文确立的统一四分法。CD36论文则区分C7、MYH11、FNDC1和CD36四种CAF，未另立DOCX中的成骨CAF类。因此本轮以DOCX四个方向为主轴，论文标记为交叉参照，并单独检验脂质方向。</p>{tab(reference_table)}
<p>DOCX中“无CNV即CAF、恶性细胞均有CNV”不能作为绝对判定规则，本轮没有据此修改细胞准入或既有身份。论文中的PD-1应答分组、OxLDL摄取、PPARG活性及ANGPTL4免疫抑制机制不直接迁移到本项目；本轮只使用对应RNA标记及功能方向。</p>
<h2>2. 新的大类结构与细胞范围</h2><p>完整资源仍为10,510个细胞，主分析10,274个非软件双细胞阳性细胞。新增列为broad_class和secondary_state；原original_class、original_state及joint_cluster不变。L编码用于避免把新工作分类与历史F类或J计算群混淆。</p>{tab(groups)}
<p><b>L3内部保留强弱差异：</b>J19和J21更符合基质合成/重塑画像；J03为基质应答过渡候选，不将其与强基质执行群等同。<b>L4的解释层级较低：</b>J01/J05/J06/J15/J17有较明确的成骨核心与背景支持；J09/J11/J13/J14/J16更多依赖成骨谱系/基质背景或混合证据，保留CAF身份未定，并单列应激等状态。细胞数不是已确认CAF数量。</p>
<h2>3. 每个J群的大类与附加状态</h2><p>命名综合原始标记、来源/患者等权评分、至少5目标细胞的患者画像敏感性及既有谱系证据，在计算新的推断P值之前冻结。没有采用“最高分必属某类”的硬性归类，也没有根据显著性反复修改大类。</p>{tab(summary_table)}
<h2>4. J07归入收缩及活化基质 并保留基质重塑叠加</h2><p>J07同时表现myCAF收缩和mCAF基质核心。至少5目标细胞患者画像中，两者平衡分数均约0.31，并伴MYH11/壁细胞成分；因此归入L2，同时保留matrix和mural附加状态。该群仍是较宽的共同群，包含原F2多个内部状态，不能将L2/J07直接改名为原F2-C3。后续内部亚群研究应在这套大类下继续开展，不能把本轮大类归纳当作已完成内部细分。</p>
{tab(profile[profile.joint_cluster.isin(['J07','J08','J03','J19','J21'])][['joint_cluster','doc_iCAF','doc_myCAF','doc_mCAF','doc_osteogenic','paper_MYH11','paper_C7','paper_FNDC1','mural_context']])}
<h2>5. 当前不单列论文中的CD36脂质CAF大类</h2><p>NCI的J10/J12/J20虽有CD36 RNA检出，但FABP4检出很少：患者等权检出比例分别约0.18%、0.88%、0.45%。同时J10/J12存在明显PTPRC信号，J20存在明显PECAM1信号。因此这三个群归入U保留层的免疫/血管竞争状态，不以CD36单基因命名为论文中的脂质CAF。其余部分成骨群的脂质面板信号保留在附表，但未据此新增独立CAF类型。</p>{tab(lipid)}
<p class="note">detection_fraction为患者等权RNA检出比例，不是细胞总体占比、蛋白阳性率或流式结果。低RNA检出不证明蛋白不存在；本轮结论仅是现有表达和身份组合不足以命名对应大类。</p>
<h2>6. 患者级统计支持与反向结果完整保留</h2><p>以下为两个新增来源中四个文献核心签名及论文四面板的主层结果。效应单位为细胞内秩评分的患者内差值，不是log fold change。比较同来源、同患者/部位该大类与其他扩展候选；其中包含身份未定和竞争层，因此这些差异不是纯CAF内部功能特异性的独立证明。小于10例均按探索性呈现，完整P/FDR与相反方向保留。</p>{tab(key[key.source.isin(['NCC','NCI'])][cols])}
<p>统计表不反向决定既已冻结的名字。没有达到名义显著、最低细胞数敏感性变弱，或某文献签名与附加状态冲突，均在完整附表中保留；若后续需改变归类，应形成新版本并说明原因。</p>
<h2>7. 评分 方法与新旧对应</h2><p>按原来源的完整计数及既有canonical feature映射提取130个标记，解决RGS5、IGF2等在合并矩阵中多标识导致的共有基因缺失问题。不可测不填零。采用UCell公开的更新后Mann–Whitney秩归一化公式，Python实现、平均秩处理并列值；主maxRank=1500，敏感性1000/2000。这里没有声称调用原生UCell软件包。已与密集全基因秩抽样核对6,678项，并对全部三种maxRank分数独立反算。公式参见<a href="https://code.bioconductor.org/browse/UCell/">UCell官方源码说明</a>。</p>
<p>先患者内求均值，再来源内患者等权，最后来源等权；同时保留单细胞数加权与每患者至少5目标细胞的画像，避免稀少来源被过度放大。基因表达和RNA检出作为另外的描述性读出。双侧符号翻转N≤16精确枚举，其余固定种子65536次蒙特卡洛加观测值；N&lt;2或全零不补造P。BH在来源×分析层×签名/单基因家族内计算。P&lt;0.05名义显著，FDR&lt;0.05才称校正后显著；不以FDR作为唯一候选门槛。</p>
<p>新旧对应保留“旧F大类→L大类”“原33亚群→L大类”“J群→L大类”三层。新大类的变化反映参考框架与新表达证据，不表示历史F分类被删除。共同J计算群仍作为下一级研究对象。</p>
<h2>8. 身份边界与验收</h2><p>L1–L3为文献支持的功能候选；L4为成骨/基质重叠表达研究层，包含原成骨及间充质未定细胞，尚不能区分全部非恶性成骨CAF与恶性成骨细胞。U保留全部混合、应激、增殖或免疫/血管竞争群，未删除细胞。正常CNV不等于已验证CAF，异常或未评估记录也未改写；各来源原CNV调用及证据等级留表。文献中的纤维沉积、免疫抑制与耐药机制不能由本轮RNA工作分类直接确认。</p>
<p>验收通过：10,510个细胞和22个J群守恒，原F/33/J标签未改；两份用户原文件副本哈希一致；所有秩分数及BH家族反算通过，新增签名检验独立复核{checked}项，最大精确P误差{maxerr:.3g}。命名冻结文件在统计前后哈希一致。通过仅表示工作注释、计算和交付范围完成，不表示CAF身份或机制已验证。</p>
<p><a href="CAF_broad_class_full_tables.html">完整大类对应及统计附表</a></p><p class="note">依据：用户《骨肉瘤CAF四分类分型研究报告.docx》2.1–2.4及5.1；用户CD36论文Fig1及第6页。原文件和提取文本、SHA256清单保存在服务器references/。证据目录：<code>{O}</code>。科学图件继续暂停，未制作客户改写版或对外发送。</p>'''
rp=O/'delivery/CAF_reference_based_broad_classes.html';rp.write_text('<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>CAF文献依据大类重新归纳</title><style>'+css+'</style><main>'+body+'</main></html>',encoding='utf-8')
frames={'大类汇总':summ,'22个J群工作注释':summary_table,'完整患者级统计':stats,'旧F到新大类':read('original_class_to_broad_counts.tsv'),'原33亚群到新大类':read('original_state_to_broad_counts.tsv'),'来源原始CNV证据':read('inherited_CNV_evidence_by_broad_class.tsv'),'面板测量覆盖':read('panel_coverage.tsv'),'maxRank敏感性':read('maxRank_sensitivity.tsv'),'CD36联合RNA检出':read('source_RNA_codetection.tsv')};payload={k:json.loads(v.to_json(orient='records',double_precision=15)) for k,v in frames.items()}
js='''const datasets=PAYLOAD;let off=0;const sel=document.getElementById('set'),q=document.getElementById('q');for(const k of Object.keys(datasets)){let o=document.createElement('option');o.textContent=k;sel.append(o)}function draw(){let data=datasets[sel.value],terms=q.value.toLowerCase().trim().split(/\\s+/).filter(Boolean),rs=data.filter(r=>terms.every(t=>JSON.stringify(r).toLowerCase().includes(t))),cols=Object.keys(data[0]||{});document.getElementById('count').textContent=`匹配${rs.length}行，显示${off+1}至${Math.min(off+60,rs.length)}`;let t=document.getElementById('table');t.replaceChildren();let tr=document.createElement('tr');for(let c of cols){let th=document.createElement('th');th.textContent=c;tr.append(th)}t.append(tr);for(let r of rs.slice(off,off+60)){tr=document.createElement('tr');for(let c of cols){let td=document.createElement('td');td.textContent=r[c]===null?'不可评估':r[c];tr.append(td)}t.append(tr)}return rs.length}sel.onchange=q.oninput=()=>{off=0;draw()};document.getElementById('next').onclick=()=>{if(off+60<draw())off+=60;draw()};document.getElementById('prev').onclick=()=>{off=Math.max(0,off-60);draw()};draw();'''.replace('PAYLOAD',json.dumps(payload,ensure_ascii=False,allow_nan=False))
ap=O/'delivery/CAF_broad_class_full_tables.html';ap.write_text('<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>CAF大类完整统计附表</title><style>'+css+'</style><main><h1>CAF大类归纳完整证据附表</h1><p><a href="CAF_reference_based_broad_classes.html">返回主报告</a>。L1–L4为四个文献功能候选方向，U是保留层而不是第五类CAF。all=所有候选；singlet=排除软件双细胞阳性；singlet_min5=同患者/部位两侧至少5细胞。N为可评估患者数，主效应为秩评分差；gene__行是log1pCP10k差。P/FDR不是命名门槛。</p><select id="set"></select><p><input id="q" placeholder="搜索，如 NCI L2 doc_myCAF"><button id="prev">上一页</button><button id="next">下一页</button></p><p id="count"></p><div class="scroll"><table id="table"></table></div><script>'+js+'</script></main></html>',encoding='utf-8')
for f in [rp,ap]:assert '\ufffd' not in f.read_text(encoding='utf-8')
archive=O/'delivery/CAF_reference_based_classes_Phase118.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in [rp,ap]:z.write(f,f.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
status=dict(status='analysis_verified_local_delivery_pending',phase=O.name,completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),functional_classes=4,reserve_layer='U',J_clusters=22,all_cells=10510,main_cells=10274,summary=summ.to_dict('records'),old_F_J_unchanged=True,annotation_type='literature-informed working functional labels; L4 CAF identity unresolved',server_archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,server_report=str(rp),figure_queue='paused_by_user',remaining_required_work=['local delivery and render'])
dump('FINAL_STATUS.json',status);print(json.dumps(qa));print(json.dumps(status,ensure_ascii=False))
