import pathlib,json,datetime,hashlib
import numpy as np,pandas as pd
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase98_reference_thresholds_20261020';old=root/'phase97_normal_background_20261020'
d=pd.read_csv(p/'tables/normal_background_all_results.tsv',sep='\t',encoding='utf-8')
cov=pd.read_csv(p/'tables/reference_threshold_coverage.tsv',sep='\t',encoding='utf-8')
checks=pd.read_csv(p/'tables/independent_R_HC2_checks.tsv',sep='\t',encoding='utf-8')
assert len(checks)==96 and checks[['effect_error','SE_error','P_error']].max().max()<1e-10
prior=pd.read_csv(old/'tables/normal_background_all_results.tsv',sep='\t',encoding='utf-8')
for a,b in [('min1','primary'),('min20','bone_fibro_exclude_1to2')]:
 x=d[d.specification.eq(a)].set_index(['source','gene']);y=prior[prior.specification.eq(b)].set_index(['source','gene']).loc[x.index]
 assert np.allclose(x[['effect','P','FDR_14']],y[['effect','P','FDR_14']],rtol=0,atol=1e-12,equal_nan=True)
for r in json.loads((p/'resources/INPUT_MANIFEST.json').read_text(encoding='utf-8')):
 assert hashlib.sha256(pathlib.Path(r['path']).read_bytes()).hexdigest()==r['sha256']
def tab(z):return z.to_html(index=False,escape=True,float_format=lambda x:f'{x:.7g}',na_rep='不可评估',border=0)
def join(x):return '、'.join(x) or '无'
summary=[];both={}
for n in [1,20,50,100]:
 g=d[d.specification.eq('min'+str(n))];ss=[]
 for src,z in g.groupby('source'):
  hit=z.loc[z.nominal & z.same_direction.eq(1),'gene'].tolist();ss.append(set(hit))
  summary.append({'阈值':n,'来源':src,'原方向且P<0.05':join(hit),'FDR<0.05':join(z.loc[z.FDR_significant,'gene']),'反向且P<0.05':join(z.loc[z.nominal & z.same_direction.eq(0),'gene'])})
 both[str(n)]=sorted(set.intersection(*ss))
rg=d[d.gene.eq('RGS2-AS1')][['specification','source','reference_bone_N','reference_lung_N','effect','P','FDR_14']].copy()
rg.columns=['方案','肿瘤来源','骨参考供者','肺参考供者','扣背景后效应','P','候选家族FDR']
h=['''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CAF参考细胞覆盖阈值敏感性分析</title><style>body{font:16px/1.7 "Microsoft YaHei",sans-serif;color:#202a35;max-width:1250px;margin:32px auto;padding:0 24px}h1{font-size:27px}h2{font-size:21px;margin-top:28px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px;border-bottom:1px solid #dbe2ea;text-align:left}th{background:#eef3f8}.wide{overflow-x:auto}.note{background:#f3f6fa;padding:14px}</style></head><body><h1>CAF参考细胞覆盖阈值敏感性分析</h1><p>2026年10月20日｜独立补充版本｜肿瘤骨肺比较的原始结果保留</p>
<p class="note"><b>主阈值为每位正常供者至少20个所选参考细胞，同时报告50、100及原始覆盖方案。</b>阈值限制的是供者的参考细胞覆盖量；没有依据基因表达或显著性选择参考。原始骨肺病灶区分仍为主问题，正常参考分析不作为剔除原候选的门槛。</p>
<h2>主要结论</h2><p>≥20个细胞时，RGS2-AS1在两个肿瘤来源中均保留原方向且P&lt;0.05；只有SCPCP000023通过候选集内FDR。提高到≥50或≥100时，SCPCP000017的P分别为0.07409和0.06669，不再达到名义显著。双来源显著性依赖阈值，不能只展示≥20的有利结果。</p><p>SCPCP000023的RGS2-AS1在全部覆盖方案下均通过候选集内FDR，是该来源中相对稳定的探索性信号。提高阈值同时使骨参考供者从5减至4、3，不能将更高阈值理解为必然更可靠。</p>
<h2>阈值及供者覆盖</h2>''',tab(cov.rename(columns={'tissue':'组织','min_cells':'最低细胞数','donors':'供者数','selected_cells':'参考细胞数','donor_ids':'纳入供者','donor_cells':'逐供者细胞数','excluded_donors':'排除供者'})),
'''<p>骨参考固定为GSE253355作者注释的Fibro-MSC；肺参考固定为Travaglini2020作者注释的成纤维细胞合并池。≥20保留5位骨供者、1439个细胞；肺端在所有阈值下均为3位供者、2136个细胞。骨阈值≥20与上轮排除仅含1–2个细胞供者的名单完全相同，重算效应、P和FDR一致。</p>
<h2>RGS2-AS1的跨阈值结果</h2>''',tab(rg),'<h2>各阈值的候选结果</h2>',tab(pd.DataFrame(summary)),
'''<h2>≥20个参考细胞时的完整结果</h2><p>效应为（肿瘤肺－肿瘤骨）－（参考肺－参考骨），单位为log2(CP10k＋0.1)尺度的均值差；正值表示扣背景后肺端较高。EBF1的扣背景后方向与原结果相反，不计入保留原方向的候选。</p>''']
for src,g in d[d.specification.eq('min20')].groupby('source'):
 z=g[['gene','original_effect','normal_background_effect','effect','CI_low','CI_high','original_P','uncorrected_Welch_P','P','FDR_14']].copy()
 z.columns=['基因','原病灶效应','参考背景效应','扣背景后效应','95%CI下限','95%CI上限','原置换P','未扣背景Welch P','本轮P','候选家族FDR']
 h+=['<h3>'+src+'</h3><div class="wide">',tab(z),'</div>']
h+=['''<h2>方法及解释范围</h2><p>沿用冻结的患者级表达：SCPCP000017骨12/肺3；SCPCP000023骨4/肺3。每位供者的选定参考细胞先汇总原始计数，再按完整文库归一化为log2(CP10k＋0.1)。在供者层面应用细胞数阈值，独立生物学样本量仍是供者数，不是细胞数。同一肺供者的区域不拆成独立重复。</p>
<p>沿用上轮四组Welch–Satterthwaite对比，包含肿瘤与参考供者间方差。每个来源、每个阈值分别对预设14候选中12项有效检验作BH校正；PEF1-AS1、AL096854.1因参考覆盖不足为不可评估，不按零表达或不显著处理。阈值间是敏感性方案，不视为独立验证，也不选择最小P作为最终结论。</p>
<p>这里的阈值不是新的逐细胞质量控制指标：原作者细胞QC与注释保留。供者≥20为本次默认覆盖标准，≥50和≥100用于检查依赖性。既往结果已知，本轮不声称盲法；阈值方案在本轮重算前记录。</p>
<p>本分析仍属于外部非肿瘤参考的探索性背景敏感性分析。骨参考来自关节置换患者的正常形态骨髓，肺参考来自局灶肺肿瘤患者的正常未受累组织；年龄、组织区室与研究来源未匹配。参考组织与研究批次不可完全分离。阈值筛选不能解决这些混杂，也不能确认纯CAF身份或肿瘤因果作用。两个肿瘤来源共用同一正常参考。</p>
<p>原报告采用精确置换检验，补充分析采用Welch近似检验。未扣背景的Welch框架本身也没有双来源同时P&lt;0.05的基因，因此不能将显著性变化全部归因于正常背景。区分骨肺病灶时，组织本底差异仍可以是有用信息；不依据本补充分析删除原14个候选。</p>
<h2>复核与保存</h2><p>96项有效对比经独立R线性模型及HC2协方差复核，效应、标准误、P值误差均小于10⁻¹⁰。原始覆盖与≥20方案与上一版衔接一致；使用的输入文件SHA-256不变。原始数据、代码、完整统计表及检查记录保存于服务器phase98_reference_thresholds_20261020。</p>
<p>参考来源：<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC11162340/">GSE253355骨髓图谱</a>；<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7704697/">Travaglini肺图谱</a>。</p><details><summary>完整统计表及检验家族</summary>''',tab(d),'</details></body></html>']
report=p/'delivery/CAF_reference_thresholds_20261020.html';report.write_text('\n'.join(h),encoding='utf-8')
status=dict(status='complete_verified',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),phase=p.name,primary_threshold=20,thresholds=[1,20,50,100],both_sources_same_direction_nominal=both,independent_checks=96,baseline_and_min20_match_previous=True,original_inputs_unchanged=True,report=str(report),report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),interpretation='Supplementary reference-coverage sensitivity; raw tumor-site discrimination remains primary')
(p/'COMPLETE.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
(root/'CAF_REFERENCE_THRESHOLD_CURRENT.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(status,ensure_ascii=False,indent=2))
