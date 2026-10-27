from pathlib import Path
import json,hashlib,re
import pymupdf
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026')
source=Path('C:/Users/31937/caf122_build_report.py').read_text(encoding='utf-8')
intro="""
p('骨肺 CAF 四项候选<br/>HLCA 补充分析','title')
heading('本补充要回答什么问题？')
p('主分析的目标，是在骨、肺病灶队列中寻找基质生成、胶原表达或纤维沉积相关的CAF亚群。原分析已得到Lung_04、Bone_07、Bone_07_s1和Lung_04_s2四项候选，并依据功能信号及患者支持优先推荐Lung_04。本补充引入外部肺组织参照，进一步解释这些候选的身份和功能。')
heading('HLCA是什么，在这里起什么作用？')
p('HLCA是Human Lung Cell Atlas，即人类肺细胞图谱，是整合多个研究的肺组织单细胞数据形成的外部参考。本次使用其健康肺核心图谱判断细胞身份，并使用包含健康肺及疾病样本的完整图谱比较表达状态。它为本项目提供已知肺细胞类型及其表达背景。')
table([['参照用途','与主分析衔接的问题','本报告采用的比较'],['身份参照','候选更接近成纤维细胞、周细胞，还是平滑肌细胞？','将候选映射到HLCA健康核心图谱，按患者汇总身份比例。'],['正常肺背景','胶原和基质表达是否超过肺组织原有成纤维细胞的水平？','比较健康肺、肿瘤邻近肺成纤维细胞的基质功能分数。'],['纤维化功能参照','候选的基质活化与已知纤维化状态重叠到什么程度？','比较IPF相关表达，并检查高分是否主要由胶原基因贡献。']],[76,207,207])
story.append(Spacer(1,12))
heading('为什么引入IPF？')
p('IPF是Idiopathic Pulmonary Fibrosis，即特发性肺纤维化。该疾病具有肺组织瘢痕化及细胞外基质沉积特征，因此可作为“纤维化相关成纤维细胞状态”的功能参照。HLCA研究曾分析IPF相关成纤维细胞，并报告CCL2、COL1A1、CTHRC1、MMP19、SERPINE1和HIF1A等相关基因；本补充据此设置六基因比较，再去除COL1A1检验胶原表达的贡献。')
p('IPF在此承担外部功能参照的角色，不能据此认定骨肉瘤患者患有IPF，也不要求CAF候选必须达到IPF状态才有研究价值。对肺侧候选，HLCA同时提供组织背景；对骨侧候选，它提供跨器官的功能参照。外部供者不计入本项目的候选患者支持数。')
p('分析关系：原骨肺队列筛选候选 → HLCA补充身份和背景解释 → IPF参照检验纤维化程序重叠 → 综合原队列患者支持确定后续验证优先级。Lung_04的推荐仍以原队列的基质功能及患者支持为主。','small')
story.append(PageBreak())
"""
anchor="p('骨肺 CAF 四项候选<br/>HLCA 补充分析','title');"
assert source.count(anchor)==1
source=source.replace(anchor,intro+"p('主要发现与候选推荐','title');",1)
source=source.replace('骨肺CAF四项候选_HLCA补充分析报告_20261027.pdf','骨肺CAF四项候选_HLCA补充分析报告_说明完善版_20261027.pdf')
source=source.replace("B/f'report_page_{i+1}.png'","B/f'explained_page_{i+1}.png'").replace("B/'REPORT_RENDER_QA.json'","B/'REPORT_EXPLAINED_RENDER_QA.json'")
source=source.replace("heading('三、IPF功能比较：分开呈现参考验证与候选表现');figure('ipf.png')","heading('三、IPF功能参照：候选与纤维化程序的重叠');p('这一步承接前面的基质升高结果：进一步判断升高是否伴随更广泛的纤维化相关表达。IPF样本提供外部疾病状态参照，候选的研究价值仍结合原骨肺队列证据判断。');figure('ipf.png')")
source=source.replace('图3｜HLCA论文Fig.6h所列6基因：','图3｜以HLCA中IPF（特发性肺纤维化）相关成纤维细胞作功能参照，使用原论文Fig.6h所列6基因：')
source=source.replace('图1｜每位患者先计算候选细胞的身份比例，','图1｜HLCA健康核心图谱用于提供已知肺细胞身份参照。每位患者先计算候选细胞的身份比例，')
source=source.replace('图2｜候选患者功能分数中位数减去健康参考的研究等权平均分数。','图2｜HLCA健康肺成纤维细胞用于判断肺组织固有基质表达背景。数值为候选患者功能分数中位数减去健康参考的研究等权平均分数。')
target=Path('C:/Users/31937/caf122_build_report_explained.py');target.write_text(source,encoding='utf-8')
exec(compile(source,str(target),'exec'),{'__name__':'__main__'})
old=Path('C:/Users/31937/Downloads/骨肺CAF四项候选_HLCA补充分析报告_20261027.pdf');new=Path('C:/Users/31937/Downloads/骨肺CAF四项候选_HLCA补充分析报告_说明完善版_20261027.pdf')
a=pymupdf.open(old);b=pymupdf.open(new);assert len(a)==6 and len(b)==7
text=''.join(page.get_text() for page in b);assert all(t in text for t in ['Human Lung Cell Atlas','Idiopathic Pulmonary Fibrosis','外部供者','0.375','0.711','97.4'])
def nums(doc,start):return re.findall(r'(?<![A-Za-z_])\d+\.\d+%?',re.sub(r'\s+','',''.join(page.get_text() for page in list(doc)[start:])))
assert sorted(nums(a,0))==sorted(nums(b,1)),'Original numerical results changed'
rec=dict(status='PASS',scope='Add HLCA/IPF definitions, rationale and connection to primary CAF analysis; clarify captions',original_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),revised_sha256=hashlib.sha256(new.read_bytes()).hexdigest(),pages=7,original_decimal_values_unchanged=True,analysis_rerun=False,original_preserved=True)
(B/'REPORT_EXPLANATION_REVISION_QA.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(rec))
