import pathlib,json,datetime,html,hashlib
import numpy as np,pandas as pd
p=pathlib.Path(__file__).resolve().parents[1]
d=pd.read_csv(p/'tables/normal_background_all_results.tsv',sep='\t',encoding='utf-8')
main=d[d.specification.eq('primary')].copy()
b=pd.read_csv(p/'tables/bone_reference_pseudobulk.tsv',sep='\t',encoding='utf-8')
l=pd.read_csv(p/'tables/lung_reference_pseudobulk.tsv',sep='\t',encoding='utf-8')
def tab(z):return z.to_html(index=False,escape=True,float_format=lambda x:f'{x:.6g}',na_rep='不可评估',border=0)
def join(v):return '、'.join(v) or '无'
records=[]
for s,g in main.groupby('source'):
 records.append(dict(source=s,evaluable=int(g.P.notna().sum()),same_direction_nominal=g.loc[g.nominal & g.same_direction.eq(1),'gene'].tolist(),nominal_reverse=g.loc[g.nominal & g.same_direction.eq(0),'gene'].tolist(),FDR_significant=g.loc[g.FDR_significant,'gene'].tolist(),unmeasured=g.loc[g.P.isna(),'gene'].tolist()))
both=set(records[0]['same_direction_nominal'])&set(records[1]['same_direction_nominal'])
bq=b[b.specification.eq('Fibro-MSC')].drop_duplicates('donor');lq=l[l.specification.eq('all_fibroblasts')].drop_duplicates('donor')
qa=json.loads((p/'resources/CONTRAST_QA.json').read_text())
original=json.loads((p/'resources/ORIGINAL_INPUT_MANIFEST.json').read_text())
assert all(hashlib.sha256(pathlib.Path(r['path']).read_bytes()).hexdigest()==r['sha256'] for r in original)
rchecks=pd.read_csv(p/'tables/independent_R_HC2_checks.tsv',sep='\t',encoding='utf-8')
assert rchecks[['effect_error','SE_error','P_error']].max().max()<1e-10
h=['''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CAF十四候选基因 非肿瘤背景敏感性分析</title><style>body{font:16px/1.7 "Microsoft YaHei",sans-serif;max-width:1250px;margin:36px auto;padding:0 24px;color:#202a35}h1{font-size:27px}h2{font-size:21px;margin-top:30px}table{border-collapse:collapse;width:100%;font-size:13px;margin:16px 0}th,td{border-bottom:1px solid #dbe2ea;padding:7px;text-align:left}th{background:#eef3f8}.notice{padding:16px;background:#f3f6fa}details,.wide{overflow-x:auto}</style></head><body>
<h1>CAF十四个候选基因的非肿瘤组织背景敏感性分析</h1><p>2026年10月20日｜补充分析｜原报告及原始结果完整保留</p>
<div class="notice"><b>本次完成的是外部非肿瘤参考背景敏感性分析，不能等同于匹配正常骨、正常肺后的无偏组织校正。</b>正常组织来自不同研究，组织与批次不可完全分离；P值以参考差异可以用于背景扣除为条件，属于探索性模型证据。</div><h2>主要结果</h2><ul>''']
for r in records:
 h.append(f"<li><b>{r['source']}</b>：可评估{r['evaluable']}/14；扣背景后同向且P&lt;0.05：{join(r['same_direction_nominal'])}；反向且P&lt;0.05：{join(r['nominal_reverse'])}；14基因家族内FDR&lt;0.05：{join(r['FDR_significant'])}。</li>")
h.append('</ul><p>两个肿瘤来源均保留原方向且扣背景后P&lt;0.05：<b>'+join(sorted(both))+'</b>。</p>')
h.append('<p><b>关键敏感性结果：</b>排除3位仅有1–2个骨Fibro-MSC的供者后，RGS2-AS1在SCPCP000017中P＝0.049604、FDR＝0.119049，在SCPCP000023中P＝0.000396、FDR＝0.004754。它在该方案下为双来源同向名义显著，但未在两个来源均通过FDR；主方案中第一个来源P＝0.078650，说明其跨来源显著性依赖骨参考覆盖范围。</p>')
h.append('<p><b>方法比较：</b>同一Welch框架在未扣背景前，也没有双来源同时P&lt;0.05的基因，不能把与原精确置换检验相比的显著性变化全部归因于正常背景。</p>')
h.append('<p>当前参考无法评估的基因：'+join(sorted(set(main.loc[main.P.isna(),'gene'])))+'。缺失特征未当作零表达补齐。</p><h2>数据与范围</h2>')
h.append(tab(pd.DataFrame([{'参考':'骨髓 Fibro-MSC','来源':'GSE253355','供者数':len(bq),'细胞数':int(bq.cells.sum())},{'参考':'肺成纤维细胞合并','来源':'Travaglini2020','供者数':len(lq),'细胞数':int(lq.cells.sum())}])))
h.append('''<p>肿瘤端沿用冻结的历史拓展层CAF候选池及原独立组患者名单：SCPCP000017骨12例、肺3例；SCPCP000023骨4例、肺3例。未加入后续NCC/NCI扩展数据。拓展层仍包含身份未解决或成骨样候选，不能等同于纯CAF。</p>
<p>骨主参考供者的细胞数为1057、126、31、55、1、1、2、170，覆盖明显不均衡；另报告排除1–2细胞供者的结果。肺主参考三位供者分别为247、525、1364个细胞。</p><h2>方法</h2><p>保留参考原作者细胞注释，按供者汇总所选细胞原始计数，以完整文库计数计算log2(CP10k＋0.1)。同一肺供者的不同区域合并，不作为独立重复。肿瘤端使用冻结患者表达矩阵。效应＝（肿瘤肺－肿瘤骨）－（参考肺－参考骨）。</p>
<p>每组至少2名有效独立供者且方差可计算时，采用四组均值的Welch–Satterthwaite对比，计入肿瘤与参考供者方差，报告双侧P及95%CI。原报告采用精确置换检验；另列未扣背景的肿瘤两组Welch P，以区分方法变化和背景扣除的影响。</p>
<p>每个肿瘤来源、每种参考方案分别在预设14基因中对有效P作BH校正。原始P&lt;0.05为名义显著；FDR&lt;0.05仅表示该14基因家族内通过校正，不能与原全转录组FDR直接比较。这是原候选筛选后的补充敏感性检验，不是独立发现验证。</p>''')
for s,g in main.groupby('source'):
 z=g[['gene','original_effect','normal_background_effect','effect','CI_low','CI_high','original_P','uncorrected_Welch_P','P','FDR_14']].copy()
 z.columns=['基因','原病灶肺−骨','参考肺−骨','扣背景后效应','95%CI下限','95%CI上限','原置换P','未扣背景Welch P','扣背景后P','14基因FDR']
 h.append('<h2>'+s+'</h2><div class="wide">'+tab(z)+'</div>')
sens=[]
for (s,spec),g in d.groupby(['source','specification']):
 sens.append({'肿瘤来源':s,'参考方案':spec,'可评估基因':int(g.P.notna().sum()),'同向且P<0.05':join(g.loc[g.nominal & g.same_direction.eq(1),'gene']),'FDR<0.05':join(g.loc[g.FDR_significant,'gene'])})
h.append('''<h2>参考选择敏感性</h2><p>primary：骨Fibro-MSC与肺全部成纤维细胞；bone_broad_MSC：较宽骨间质细胞范围（排除RNAlo）；lung_alveolar与lung_adventitial：肺分别限定肺泡型或外膜型成纤维细胞。前四个方案在查看结果前登记。bone_fibro_exclude_1to2为质量核查后追加：排除仅有1–2个Fibro-MSC的H33、H35、H39，剩余5位供者各有31个及以上细胞；主方案不变。</p>'''+tab(pd.DataFrame(sens)))
h.append('''<h2>检验方法与背景扣除的区别</h2><p>在同一个Welch框架内，尚未扣除背景时，也没有基因在两个肿瘤来源中同时P&lt;0.05。因此，本报告主方案的双来源未同时显著，不能全部归因于正常组织背景；必须同时看效应变化、检验方法及参考方差。原精确置换结果仍保留。</p><h2>解释边界</h2><p>本分析回答“减去这组外部非肿瘤参考的肺骨差异后，候选信号怎样变化”。骨参考为关节置换患者的形态正常骨髓，来源队列52–80岁，并非年龄匹配健康骨皮质。肺参考为3名局灶肺肿瘤手术患者的病理正常未受累组织，1名有轻度肺气肿，并非健康志愿者。年龄、组织区室、取材及技术批次未匹配。</p><p>正常组织部位与研究来源混杂，不能靠一个批次协变量或整合表达值得到可识别的无偏肿瘤×组织交互效应。结果只能用于探索性背景敏感性评价，不能证明肿瘤因果作用、CAF特异性或已消除所有正常组织差异。肺参考仅3位供者，应结合方向、效应及参考方案稳定性解释。两个肿瘤来源复用了同一正常参考，不能作为两次独立正常参考验证。Welch P依赖小样本下的近似分布假设。</p><h2>复核</h2>''')
h.append(f"<p>原患者效应复算最大误差{qa['original_effect_max_error']:.3g}。R线性模型与HC2异方差协方差独立复核{len(rchecks)}项有效检验，效应、标准误及P误差均小于1×10⁻¹⁰。原输入SHA-256不变。原始参考通过压缩完整性及计数、文库分母核查。完整统计家族与缺失记录保存在服务器项目目录。</p>")
h.append('''<h2>来源</h2><p><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11162340/">骨髓图谱：Bandyopadhyay等，Cell 2024，GSE253355</a><br><a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7704697/">肺图谱：Travaglini等，Nature 2020</a><br><a href="https://cells.ucsc.edu/stanford-czb-hlca/">UCSC原始UMI与注释</a></p><details><summary>完整补充统计表</summary>'''+tab(d)+'</details></body></html>')
out=p/'delivery/CAF_normal_background_20261020.html';out.write_text('\n'.join(h),encoding='utf-8')
status=dict(status='completed_verified_external_reference_sensitivity_only',phase=p.name,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),report=str(out),primary_results=records,both_sources_same_direction_nominal=sorted(both),strict_matched_normal_adjustment=False,original_inputs_unchanged=True,independent_R_checked_tests=len(rchecks),normal_bone_donors=len(bq),normal_lung_donors=len(lq))
(p/'COMPLETE.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(status,ensure_ascii=False,indent=2))
