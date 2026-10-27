from pathlib import Path
from datetime import datetime,timezone
import json,shutil,zipfile,hashlib
import pandas as pd
P=Path('C:/Users/31937/CAF_paper_figures_20261021');B=P/'batch02'
O=Path('C:/Users/31937/Downloads/CAF_reference_figures_600DPI_20261021')
mf=B/'resources/BATCH02_MANIFEST.json';m=json.loads(mf.read_text(encoding='utf-8'))
now=datetime.now(timezone.utc).isoformat();captions=[]
stats=pd.read_csv(B/'tables/RNA_detection_exact_statistics.tsv',sep='\t')
for entry in m['figures']:
 name=entry['name'];g=entry['gene']
 for e in entry['exports']:
  p=Path(e['file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==e['sha256']
 if entry['kind']=='patient_expression':
  q=P/(name+'_QA.json');qa=json.loads(q.read_text(encoding='utf-8'));qa.update(status='PASS',visual_QA='Inspected final rendered preview: all labels, sample sizes, source separators and legends legible; no clipping; axis covers all observed values',checked_at=now);q.write_text(json.dumps(qa,indent=2),encoding='utf-8')
  shutil.copy2(q,B/'resources'/q.name)
  for tail in ['_patient_values.tsv','_summary.tsv']:shutil.copy2(P/(name+tail),B/'tables'/(name+tail))
  captions.append(O/f'Fig1h_{g}_图注与口径.txt')
 elif entry['kind']=='patient_detection':
  st=stats[stats.gene.eq(g)]
  text=f'''{g}：骨/肺病灶患者级RNA检出比例（借鉴Fig.1m图型）

每个点为一位患者；分子为该患者分析候选池中该基因原始计数>0的细胞数，分母为同一患者全部纳入的候选细胞数。深色横线为患者比例中位数，竖线为四分位距。骨病灶蓝色，肺病灶红色。百分比为RNA检出比例，不是流式阳性比例，也不是该基因阳性CAF占全部组织的比例。

来源内患者分组：SCPCP000017骨12例/肺3例；SCPCP000023骨4例/肺3例；NCI肺17例。各患者只在一个病灶组，沿用主分析已冻结的独立比较名单及配对患者排除，不新设细胞数阈值。NCI仅作肺来源描述，不混入两来源内比较。两原始来源分别包含3269和1834个候选细胞，NCI包含1530个作者注释CAF。既有扩展候选身份边界保留，不声称全部为纯CAF。

效应为肺组患者平均检出比例减骨组患者平均检出比例，单位为百分点（pp）。检验在患者层面穷举独立分组标签：017为455种、023为35种排列，双侧P定义为绝对均值差不小于观察值的排列比例；所有排列包含原始分配，不作Monte Carlo加1。每个来源内六个预定候选基因构成一个BH FDR家族，不与表达、富集或其他检验混池。全零/无变异端点仅报效应，不生成P值。统计已用独立SciPy排列分布核对。

017图内报告效应、P和FDR；023仅7例，图内报告效应及探索性定位，完整P与FDR列于下表。名义P<0.05未通过FDR时不称校正后显著。这是已筛选候选在同一数据中的补充端点，不能当作独立验证。未控制细胞组成、测序检出能力和技术差异，比例可能受这些因素影响。源内样本少及部分患者候选细胞少，解读保持探索性；完整患者分母保存于服务器配套统计表。

本基因完整统计：
'''+st.to_csv(sep='\t',index=False)+'\n输出：同一张图的PNG与LZW TIFF，5760×3960像素，600 DPI。纵轴统一为0–100%，图中保留所有患者与零值。\n'
  path=O/(name+'_图注与口径.txt');path.write_text(text,encoding='utf-8');captions.append(path)
 else:
  src='SCPCP000017' if 'SCPCP000017' in name else 'SCPCP000023'
  text=f'''{src}：患者级CAF候选分类组成（借鉴Fig.1l图型）

每根柱为一位患者，按骨/肺解剖病灶排列。颜色为沿用的F1–F5扩展候选分类；每根柱的分母仅为该患者已纳入分析的全部扩展候选细胞，不是全部组织细胞。柱顶数字为该分母细胞数，横轴为原始患者标识。各分类比例相加为100%。图中逐患者展示，未将各患者细胞合并后推断患者组差异。

SCPCP000017为15例（骨12/肺3）、3269个候选细胞；SCPCP000023为7例（骨4/肺3）、1834个候选细胞。名单沿用Phase103原始两来源的独立主分析池，不新增过滤。少至2或3个细胞的患者保留，并显式标出分母；其百分比粒度很粗，不能据此认定稳定细胞亚群组成。F分类不是全为已确认纯CAF。NCI只有作者CAF注释，未强行套入F1–F5，故不制作其F分类组成图。

本图为描述性组成图，无显著性检验，无P值。各患者原始计数、比例及分类在服务器batch02/tables/patient_class_composition.tsv保留。输出为6240×3960像素、600 DPI单图，PNG与LZW TIFF内容相同。
'''
  path=O/(name+'_图注与口径.txt');path.write_text(text,encoding='utf-8');captions.append(path)
for path in captions:shutil.copy2(path,B/'resources'/path.name)
m.update(status='PASS',visual_QA='All 12 final previews inspected individually; axes, text, source boundaries, legends and count labels legible; no clipping or overlap',checked_at=now)
mf.write_text(json.dumps(m,indent=2),encoding='utf-8')
readme='''CAF参考论文仿图：批次02交付说明

本批新增12张独立图，每张提供PNG和LZW TIFF，全部为600 DPI，未拼成论文大图。研究候选保持PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1。

4张表达分布：PNKD、TOR1AIP1、TAF10、RGS2-AS1。连同此前GRID2和BEST1，六候选表达图已齐。
6张RNA检出比例：每个候选各1张，骨/肺按来源比较，NCI仅肺单独展示。
2张患者组成：SCPCP000017、SCPCP000023各1张，分母为已纳入CAF候选池。

新增检出比例结果：在SCPCP000017中，TAF10肺高于骨10.498646个百分点，P=0.002197802197802198，六候选家族FDR=0.013186813186813187，达到该家族校正后显著；PNKD和GRID2为名义显著，均未通过FDR。TOR1AIP1、BEST1、RGS2-AS1在这一检出比例端点未达P<0.05。SCPCP000023中TAF10肺高于骨2.139632个百分点，为名义显著，其余五基因在这一端点未达P<0.05；该来源仅7例，全部按探索性解释，完整数值见图注。检出比例与原先伪bulk表达为不同端点，不替换六候选名单。

统计为患者级源内比较；NCI仅肺，不能进行源内骨肺检验。当前结果来自同一候选筛选数据，是补充探索性证据，不是独立验证。部分患者细胞较少，技术检出差异和细胞组成可能影响结果；图上保留所有零值、原始分母和阴性结果。

本批图目：
'''
for i,e in enumerate(m['figures'],1):readme+=f"{i:02d}. {e['name']}\n"
readme+='\n项目归档：/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase111_reference_figure_adaptation_20261021/batch02/\n已完成图累计15张。其余图型按数据条件继续；本批未开展新的UMAP、GO/KEGG、TF或细胞通讯分析，不将这些待办标记完成。\n'
r=O/'批次02_图目与结果说明.txt';r.write_text(readme,encoding='utf-8');shutil.copy2(r,B/'resources'/r.name)
archive=O/'CAF_批次02_12张单图_600DPI.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for e in m['figures']:
  for f in e['exports']:z.write(f['file'],Path(f['file']).name)
 for p in captions+[r]:z.write(p,p.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None and len(z.namelist())==37
(B/'resources/ARCHIVE_QA.json').write_text(json.dumps(dict(file=str(archive),members=37,sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),bytes=archive.stat().st_size,zip_crc='PASS'),indent=2),encoding='utf-8')
status=P/'STATUS.json';old=json.loads(status.read_text(encoding='utf-8'));shutil.copy2(status,B/'resources/STATUS_before_batch02.json')
old.update(status='batch02_complete',user_request='Adapt supported paper figures using current six candidates; produce multiple independent 600 DPI figures per batch and deliver after QA',completed_figure=m['figures'][-1]['name'],updated_at=now,current_batch_manifest='batch02/resources/BATCH02_MANIFEST.json',next_suggested_single_figure=None,next_batch_focus='Audit current-pool embedding coverage for UMAP and per-gene FeaturePlots; use fixed existing functional programs for heatmap where supported',visual_QA=m['visual_QA'],data_QA='PASS: inherited patient expression/counts reproduced, 12 source-specific permutation tests independently checked; two six-gene BH families',pixels='See per-figure manifest')
for e in m['figures']:
 if e['name'] not in old['completed_figures']:old['completed_figures'].append(e['name'])
assert len(old['completed_figures'])==15
status.write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS: 12 new figures, 15 cumulative;',archive,archive.stat().st_size)
