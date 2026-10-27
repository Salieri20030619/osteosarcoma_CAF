import pathlib,json,base64,hashlib,datetime
import pandas as pd
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase100_NCI_completion_20261020'
b=json.loads((p/'BUILD_COMPLETE.json').read_text(encoding='utf-8'))
qa=json.loads((p/'resources/PROCESSING_QA.json').read_text(encoding='utf-8'))
ncc=json.loads((p/'resources/NCC_SITE_SEARCH_CURRENT.json').read_text(encoding='utf-8'))
assert ncc.get('supplementary_pdf_retrieval_finalized',False),'NCC supplementary retrieval still pending'
s=pd.read_csv(p/'tables/final_library_disposition.tsv',sep='\t',encoding='utf-8')
r=pd.read_csv(p/'tables/cluster_review_decisions.tsv',sep='\t',encoding='utf-8')
gs=pd.read_csv(p/'tables/NCI_14gene_descriptive_summary.tsv',sep='\t',encoding='utf-8')
patient=pd.read_csv(p/'tables/NCI_sample_patient_crosswalk_final.tsv',sep='\t',encoding='utf-8',keep_default_na=False)
expr=pd.read_csv(p/'tables/NCI_14gene_library_descriptive.tsv',sep='\t',encoding='utf-8')
def tab(x):return x.to_html(index=False,escape=True,float_format=lambda x:f'{x:.5g}',na_rep='未评估',border=0)
def pic(name,alt):return '<img alt="'+alt+'" src="data:image/png;base64,'+base64.b64encode((p/'figures'/name).read_bytes()).decode()+'">'
col={'sample':'文库','cells':'输入细胞核','reference_pass':'参考诊断通过','expanded_candidates':'拓展候选','ambient_supported':'去污染支持子集','author_CAF_agreement':'作者CAF一致子集','author_tumor_conflict':'作者肿瘤标签冲突','author_tumor_excluded':'不含作者肿瘤标签','putative_group':'编号归并假设'}
other=b['expanded_candidates']-b['author_CAF_agreement']-b['author_tumor_conflict']
h=['''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>NCI肺转移候选资源整理与NCC取材信息核查</title><style>body{font:16px/1.75 "Microsoft YaHei","Noto Sans CJK SC",sans-serif;color:#202a35;max-width:1200px;margin:34px auto;padding:0 26px}h1{font-size:28px}h2{font-size:22px;margin-top:32px}h3{font-size:18px}img{width:100%;height:auto;margin:16px 0}table{width:100%;border-collapse:collapse;font-size:13px}th,td{padding:7px;border-bottom:1px solid #dbe2ea;text-align:left;vertical-align:top}th{background:#edf3f8}.lead{background:#f0f5f8;padding:16px}.wide,details{overflow-x:auto}small{color:#566574}</style></head><body>
<h1>NCI肺转移候选资源整理与NCC取材信息核查</h1><p>2026年10月20日｜内部分析记录｜phase100版本</p>
<div class="lead"><b>NCI现有计数处理与候选复核已完成。</b>19份肺转移单核文库共123,255个细胞核，14份通过冻结的参考诊断，形成10,627个身份分层的拓展候选及配套计数对象。独立患者对应关系尚未确认，本版不把文库数或编号归并组数作为生物学样本量，也未生成加入NCI后的骨肺推断P值。</div>
<h2>结果</h2><h3>14份文库形成可追溯的拓展候选资源</h3>
<p>18份文库完成原生CopyKAT计算，1份因参考不足未进入该计算；在已计算文库中，14份通过参考诊断，4份未通过。全部19份已有明确处理状态。295个群均已记录保留、排除或因文库参考不合格而不释放的处理决定；14个候选计数对象均通过保存后回读核查。</p>''',
tab(s[['sample','cells','reference_pass','expanded_candidates','ambient_supported','author_CAF_agreement','author_tumor_conflict']].rename(columns=col)),
'<p>去污染支持9,710个，是拓展候选内部的证据子集，不能与总数相加。5份未释放候选的文库按既定参考诊断保留排除或覆盖不足结论，没有为增加样本量放宽门槛。</p>',
pic('NCI_candidate_disposition.png','NCI各文库候选数及与作者注释的一致性'),
f'''<h3>拓展候选与作者CAF注释存在大量分歧</h3><p>10,627个拓展候选中，478个与作者CAF标签一致；8,556个（{b['author_tumor_conflict']/b['expanded_candidates']*100:.2f}%）与作者肿瘤标签冲突；其余{other}个属于作者其他宿主标签。CCDI0097_8978贡献5,702个候选，其中5,367个带作者肿瘤标签冲突。该差异已逐细胞保留。</p>
<p>本版同时提供拓展层、去污染支持层、作者CAF一致层、排除作者肿瘤标签层的描述性表达表。原作者标签与现有门控均作为证据，不把任一方单独视为已确认细胞真值。</p>
<h3>14个候选基因已形成文库级与编号归并假设下的描述性表达表</h3><p>作者计数对象包含全部14个目标基因条目。拓展候选中，各基因获得非零读数的文库数如下；该数字表示表达覆盖，不是独立患者数或检验显著性。</p>''',
tab(gs[gs.stratum.eq('expanded')][['gene','libraries','libraries_with_positive_candidate_counts','min_expression','median_expression','max_expression']].rename(columns={'gene':'基因','libraries':'有候选文库数','libraries_with_positive_candidate_counts':'候选中非零读数文库数','min_expression':'表达最小值','median_expression':'文库表达中位数','max_expression':'表达最大值'})),
pic('NCI_14gene_descriptive.png','14个候选的文库级表达以及作者CAF一致子集'),
'''<p>热图单位为log2(CP10k＋0.1)，使用原作者现有计数矩阵的文库总量归一化。作者CAF一致子集在部分文库中仅有少量细胞，其完整细胞数见上表。热图不包含患者级显著性检验。</p>
<h3>肺取材已确认，独立患者数仍需逐样本对应证据</h3><p>导出的19份文库均与作者样本表匹配，location均为lung，method均为single_nucleus。19个文库编号按后缀形成17个归并假设；有拓展候选的14份文库对应12个归并假设。16178/16178b、16183/16183b的重复关系不能仅凭名字确认。</p>
<p>先前对两对文库共同条形码的21,332特征计数比较未发现完全相同向量，说明未见这些条形码的原样复制，不能据此证明来自不同患者。本轮查到phs004490.v1.p1及相应NCI单核多组学预印本的公开研究线索；研究级22位受试者、24份标本不能直接套用到当前19份文库。所列代码仓库公开API返回404，尚未取得逐样本患者对应表。</p>
<h3>NCC仍缺少逐BS编号的实际解剖取材对照</h3><p>已核对当前作者仓库、Patient metadata、作者Issue2回复、2024年正文及补充方法、2025年预印本v2及本轮补充图获取记录。现有Location条目为Primary/Recurrent，表示肿瘤状态；预印本中的骨化/未骨化区域描述属于组织学，不构成骨／肺解剖取材对应。</p>
<p>作者已澄清BS14属于Patient14且为治疗前；BS4与BS15才是Patient4的治疗前后样本。当前纳入的7份NCC文库不含BS15，因此旧“6位”登记应按作者更正改为7位来源内患者。本轮保留phase99更正对照，不覆盖历史入库基准。仍没有新增可确认进入骨肺分组的NCC样本。</p>''',
'<p>补充图检查状态：'+str(ncc.get('supplementary_pdf_result','未记录'))+'</p>',
'''<h2>方法</h2><p>输入为phase95已导出的19份NCI肺转移文库，共123,255个细胞核。复用该阶段的谱系评分、scDblFinder、CopyKAT和decontX结果。机械门控为：唯一正向最高Phase27成纤维评分、scDblFinder单细胞、CopyKAT二倍体且文库参考诊断通过。正常参考诊断与后续≥20个正常参考细胞的背景敏感性标准分别记录。</p>
<p>群级复核使用描述性标志基因、竞争谱系评分、作者注释背景及候选子集内38个身份标志的检出比例和平均表达。明确的免疫、内皮、上皮或血管壁竞争谱系群不纳入；符合历史拓展层口径的基质、成骨或混合状态子集保留，并附身份冲突标记。判定清单、证据表和逐细胞处理结果分别保存。</p>
<p>按文库和身份分层汇总原始计数，计算log2(CP10k＋0.1)。另按编号后缀假设归并计数形成描述性表，归并状态与患者核实状态分别编码。14个输出计数对象回读后与对应原始子矩阵逐值比较，输出和来源证据记录SHA-256。</p>
<h2>讨论与后续分析边界</h2><p>本版完成的是现有输入下的候选资源处理。群级判定属于AI辅助表达证据复核，没有独立病理或实验确认。CopyKAT二倍体、成纤维表达或去污染支持都不能单独证明纯CAF身份；本批与作者肿瘤标签的分歧尤其需要保留。原作者已处理计数对象不等于完整provider-called原始矩阵，原始逐文库测量特征掩码仍不可用。</p>
<p>NCI本批只有肺端。即使后续独立患者关系得到确认，也不能简单与另一个来源的骨端拼接后，把来源或技术差异当成部位效应。正常参考阈值仍保留每位供者至少20个所选参考细胞，但本轮未对NCI生成该背景扣除后的患者级P值。原骨肺主分析与此前正常参考敏感性分析保持原版本；数值未变不能被解释为新增患者验证通过。</p>
<h2>剩余条件与证据</h2><table><tr><th>编号</th><th>级别／影响</th><th>证据及需要的行动</th></tr>
<tr><td>NCI-PAT-01</td><td>阻断新增患者级推断</td><td>作者样本表无独立patient ID列，两个重复编号配对未证实。需19份文库到匿名患者ID、标本ID、取材时间点及技术重复的对应表，并核查与原队列重叠。验收为每份文库有可追溯患者映射，非由表达或编号猜测。</td></tr>
<tr><td>NCI-ID-01</td><td>高／限制纯CAF与临床解释</td><td>8,556个拓展候选与作者肿瘤标签冲突。逐细胞冲突标记及作者CAF一致层已保留；不得将全部10,627个写成确认CAF。强化身份结论需独立证据，而非删除冲突记录。</td></tr>
<tr><td>CAF-ADD-01</td><td>待补信息／阻断NCC骨肺分组</td><td>BS1、BS3、BS4、BS5、BS7、BS12、BS14缺实际测序组织部位对照。需对应Bone/Lung/Other的标本级证据；Primary/Recurrent或骨化组织学不能替代。</td></tr>
<tr><td>CAF-ADD-02</td><td>作者证据已明确，历史基准未覆盖</td><td>BS14映射Patient14；BS4/BS15才为Patient4配对。phase99已有更正工作对照；未来NCC重分析应用该有依据映射，不能沿用BS4/BS14错误合并。</td></tr></table>
<h2>验收与文件位置</h2><p>已检查19份输入文库处理状态、295个群的处理决定、14个计数对象回读与哈希、原始证据文件未变化、候选数量对账以及新旧准入状态字段一致性。两张图已检查标签与版式。该验收仅针对本版处理范围，不代表独立患者身份或纯CAF身份已验证。</p>
<p>服务器目录：/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase100_NCI_completion_20261020。原始计数、代码、日志、身份分层、完整表达表、患者对应状态和NCC搜索证据留在该目录；本HTML为自包含内部报告。</p>
<h2>公开来源</h2><p><a href="https://github.com/kidcancerlab/Osteo_atlas">OstEoCAR作者代码及样本表</a>；<a href="https://doi.org/10.64898/2026.01.13.698472">OstEoCAR预印本</a>；<a href="https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs004490.v1.p1">NCI snMultiome研究登记</a>；<a href="https://doi.org/10.64898/2026.08.11.744222">NCI直接研究预印本线索</a>；<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11337771/">NCC原文</a>；<a href="https://doi.org/10.1101/2025.06.27.661729">NCC后续预印本</a>；<a href="https://github.com/zhengxj1/A-Single-Cell-and-Spatially-Resolved-Atlas-of-Human-Osteosarcomas/issues/2#issuecomment-6194375996">NCC作者更正</a>。</p>
<details><summary>逐文库患者对应状态</summary>''',tab(patient),'</details><details><summary>完整逐群处理决定</summary>',tab(r),'</details><details><summary>14基因文库级描述性表达</summary>',tab(expr),'</details></body></html>']
out=p/'delivery/NCI_completion_and_NCC_site_review_20261020.html';out.write_text('\n'.join(h),encoding='utf-8')
status=dict(status='existing_NCI_input_processing_complete_verified;patient_inference_not_released',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase=p.name,report=str(out),report_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),processing= b,QA=qa,patient_crosswalk='not_verified',NCC_anatomical_crosswalk='not_found',new_inferential_tests=False,normal_reference_threshold_retained=20)
(p/'FINAL_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in status.items() if k not in ['processing','QA']},ensure_ascii=False,indent=2))
