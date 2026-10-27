from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile
import numpy as np,pandas as pd
P=Path('C:/Users/31937/CAF_paper_figures_20261021/CK');O=Path('C:/Users/31937/Downloads/CAF_Fig1c_Fig1k_600DPI_20261021')
m=json.loads((P/'QA/MANIFEST.json').read_text(encoding='utf-8'));cfg=json.loads((P/'PANELS.json').read_text(encoding='utf-8'));qa=json.loads((P/'DATA_QA.json').read_text(encoding='utf-8'))
c=pd.read_csv(P/'C_dotplot_values.tsv',sep='\t');assert c.mean_log1pCP10k.ge(0).all() and c.detection_percent.between(0,100).all() and c.cells.eq(c.measured_cells).all()
for src,total in [('SCPCP000017',167180),('SCPCP000023',178039)]:assert c[c.source.eq(src)].groupby('gene').cells.sum().eq(total).all()
for typ in ['gene','program']:
 d=pd.read_csv(P/f'K_{typ}_patient_class_profiles.tsv',sep='\t');assert len(d)==59 and d.cells.sum()==5103
 assert not d[['source','class','donor']].duplicated().any()
for zfile in (P/'QA').glob('*zscores.tsv'):
 z=pd.read_csv(zfile,sep='\t',index_col=0).to_numpy();axis=0 if '_C_' in zfile.name else 1
 assert np.isfinite(z).all() and np.max(np.abs(np.mean(z,axis=axis)))<1e-10
 sd=np.std(z,axis=axis);assert np.all((sd<1e-10)|(abs(sd-1)<1e-10))
assert len(m['figures'])==6
files=[]
for f in m['figures']:
 for e in f['exports']:
  p=Path(e['file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==e['sha256'];files.append(p)
text='''Fig.1c / Fig.1k对应图：中文图注与方法（2026-10-21）

交付范围
SCPCP000017与SCPCP000023各3张：全细胞marker气泡图、CAF候选功能基因热图、固定功能表达程序热图，共6张。均为600 DPI单图，PNG与LZW TIFF为同一图的两种格式。NCI未混入本批；未修改原始注释、候选名单或旧统计结果。

Fig.1c：全细胞marker气泡图
与上一套Fig.1b全细胞UMAP使用同一细胞范围和冻结表达谱系注释：SCPCP000017为167180细胞、20个来源内患者记录、19类注释；SCPCP000023为178039细胞、26个来源内患者记录、15类注释。未定、混合谱系及稀少类型均保留。每个来源展示相同的39个marker，来自项目已冻结的13个谱系面板，基因是在此次查看表达结果前固定，不按图上分离效果反复筛选。
每个细胞按全转录组原始UMI归一化为CP10k后取自然对数log1p。先在同一患者/细胞类型内汇总全部样本和病灶的细胞，计算平均log1p(CP10k)及原始计数>0的细胞比例；再对含有该类型细胞的患者等权平均。同一患者的多个样本不会作为独立患者。不存在该细胞类型的患者不被虚构为“该类型表达为零”的患者。纵轴n表示含该类型细胞的患者数。
气泡面积与患者等权RNA检出比例成正比；颜色为每个基因在该来源不同细胞类型间平均表达的z-score（总体标准差ddof=0），显示范围-2至2，超出范围使用端点颜色；原始值和未截断z-score保留。色标表示同一基因在类型之间的相对分布，不能由颜色比较不同基因的绝对表达量，也不能跨来源直接比较颜色强度。零检出不画气泡。所有39个marker在本次样本均有完整测量覆盖，不以未测量代替零值。
本图为冻结注释的marker展示，不独立证明恶性身份或确诊CAF身份；若细分类型共享marker或表达不典型，仍保留冻结标签，不在本次画图中静默改注释。图上13个marker面板与19/15个细分注释不是一一对应。稀少患者类型仅作描述，无细胞级P值。

C图marker分组：
'''
for k,gs in cfg['C'].items():text+=k+': '+', '.join(gs)+'\n'
text+='''
Fig.1k：CAF候选功能基因热图
沿用已验收六候选绘图使用的扩展候选池，选其中两个原始来源：SCPCP000017为3269细胞/15例，SCPCP000023为1834细胞/7例，共5103细胞。NCI的1530细胞没有既有F1–F5等价分类，本批不混入。K图不是把全细胞都当作CAF，也不是C图的同一患者名单。F1–F5沿用既有扩展候选分类，不声称全部为已确认纯成纤维细胞；双细胞及身份边界保持原分析口径。
每列为一个患者在对应F分类中的平均log1p(CP10k)表达谱，先逐细胞使用完整转录组UMI归一化，再求患者/分类内均值。017为45列，023为14列，共59个患者/分类表达谱，不能把59当作独立患者数。相同患者可跨分类重复出现。列按既有F类别、患者标识排序，未按结果调聚类；上方分别标记分类和骨/肺病灶，蓝为骨、红为肺。各分类n为该分类中出现的患者数，少量细胞及单患者分类全部如实保留，不能据此作总体推断。每个表达谱的实际细胞数保存在配套表中。
功能示例基因从已冻结程序成员中选择，覆盖胶原成分、胶原加工、交联、ECM重塑、收缩、矿化和炎症7组，共26个；另加6个当前候选基因，合计32行。TAGLN与SPP1不属于本次采用的相应冻结程序成员，未作为该组示例加入；这不是根据阴性结果排除。六候选均保留原始标识，不改变研究名单。
颜色为每行在同来源所有患者/分类表达谱间的z-score，ddof=0，显示范围-2.5至2.5。常数行z-score置0，其原始零值或常数值在数据表中保留。颜色仅反映该行相对差异；不能把红色等同统计显著、通路活性或临床作用，也不能跨基因/跨来源直接比较绝对表达。原始均值和未截断z-score均保留。

K图行分组：
'''
for k,gs in cfg['K'].items():text+=k+': '+', '.join(gs)+'\n'
text+='''
Fig.1k补充：固定功能表达程序热图
使用同一患者/分类列顺序，继承Phase66冻结的11个raw表达程序分数，先在患者/分类内求均值，再按行作同来源z-score。程序分数是冻结成员的log1p(CP10k)均值，沿用原程序覆盖与计分规则，不以程序名称代替实测功能。程序之间有基因重叠，程序图与成员基因图不是独立验证；候选基因若为某程序成员，其关联也不能当作独立功能证据。LOX_family是靶向RNA面板，不是酶活性测量。所有图为描述性展示，没有重新生成P或FDR。

冻结程序及成员数：
'''
for k,gs in cfg['programs']['sets'].items():text+=f'{k}: {len(gs)}\n'
text+='''
验收
核对全细胞样本细胞ID与完整计数矩阵一一匹配，逐样本完整UMI与既有n_counts一致。C图各基因覆盖与分母一致，无缺失测量。K图六候选的逐细胞表达重建与此前已验收数据一致；59列表达谱与程序谱列标识相同，合计5103细胞。z-score均值与标准差完成数值检查。6张图均回读像素和DPI，并逐张检查标签、色标和注释条；已修正基因标签与行分组条的遮挡。

可复核文件均在sgs10项目：
/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase113_marker_function_heatmaps_20261021/
tables/C_patient_type_gene.tsv.gz：患者/类型/marker表达与检出率。
tables/C_dotplot_values.tsv：C图全部绘图值。
tables/K_gene_patient_class_profiles.tsv：32基因患者/分类原始均值及细胞数。
tables/K_program_patient_class_profiles.tsv：11程序患者/分类均值及细胞数。
resources/PANELS.json：本批marker、功能示例及全部冻结程序成员。
resources/C_feature_mapping.tsv、K_feature_coverage.tsv：特征映射与测量覆盖。
QA/：未截断z-score、输出哈希、像素和DPI。
'''
cap=O/'Fig1c_Fig1k_中文图注与方法.txt';cap.write_text(text,encoding='utf-8');files.append(cap)
index=O/'图件索引.txt';index.write_text('本批6张独立单图，600 DPI；每张PNG与TIFF内容相同。\n\n'+'\n'.join(f['name'] for f in sorted(m['figures'],key=lambda x:x['name']))+'\n\nC图使用全细胞；K图使用既有扩展CAF候选池。请结合中文图注阅读患者数、色标和数据边界。\n',encoding='utf-8');files.append(index)
m.update(status='PASS',data_QA=qa,visual_QA='All final six previews individually inspected; gene labels, marker labels, group boundaries, colorbars and legends legible without clipping or overlap',numeric_QA='Counts/coverage and z-score normalization checked',completed_at=datetime.now(timezone.utc).isoformat())
(P/'QA/MANIFEST.json').write_text(json.dumps(m,indent=2),encoding='utf-8')
archive=O/'CAF_Fig1c_Fig1k_6张单图_600DPI.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,f.name)
with zipfile.ZipFile(archive) as z:assert len(z.namelist())==14 and z.testzip() is None
receipt=dict(status='completed_verified',phase='phase113_marker_function_heatmaps_20261021',completed_at=m['completed_at'],figures=6,dpi=600,C_cells=345219,K_cells=5103,sources=['SCPCP000017','SCPCP000023'],C_markers=39,K_function_genes=26,K_candidates=6,K_programs=11,visual_QA='PASS',numeric_QA='PASS',archive=str(archive),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_bytes=archive.stat().st_size,archive_members=14,original_results_modified=False)
(P/'QA/FINAL_STATUS.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(receipt,ensure_ascii=False))
