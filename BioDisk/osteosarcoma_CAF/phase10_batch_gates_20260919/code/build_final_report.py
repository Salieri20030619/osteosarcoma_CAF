from pathlib import Path
import json,html,tarfile,hashlib,io
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4,A3,landscape
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfgen import canvas
from pypdf import PdfReader,PdfWriter
from PIL import Image
import pypdfium2 as pdfium
B=Path(r'C:\Users\31937\caf_batch_gates_20260919');O=B/'results/phase10_batch_gates_20260919';D=O/'delivery';Q=B/'pdf_qa'
if not (O/'MANIFEST_SHA256.json').exists():
 with tarfile.open(B/'analysis_delivery.tar.gz',encoding='utf-8') as t:t.extractall(B/'results',filter='data')
m=json.loads((O/'MANIFEST_SHA256.json').read_text(encoding='utf-8'))
for name,h in m.items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==h,name
D.mkdir(exist_ok=True);Q.mkdir(exist_ok=True);status=json.loads((O/'STATUS.json').read_text());changes=pd.read_csv(O/'tables/cohort_gate_changes.tsv',sep='\t');metrics=pd.read_csv(O/'tables/batch_model_comparison.tsv',sep='\t').set_index('analysis');support=pd.read_csv(O/'tables/primary_cluster_review.tsv',sep='\t');loo=pd.read_csv(O/'tables/leaveout_summary.tsv',sep='\t')
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'));styles={k:ParagraphStyle(k,fontName='STSong-Light',fontSize=sz,leading=ld,spaceAfter=7,wordWrap='CJK',keepWithNext=k in ['title','head']) for k,sz,ld in [('title',19,27),('head',13,19),('body',10,15),('cell',8,11)]};story=[];md=[]
def add(t,k='body'):story.append(Paragraph(html.escape(t),styles[k]));md.append(('# ' if k=='title' else '## ' if k=='head' else '')+t+'\n')
def table(headers,rows,widths):
 md.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)+'\n');t=Table([[Paragraph(html.escape(str(x)),styles['cell']) for x in row] for row in [headers]+rows],colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e4edf3')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,-1),.25,colors.HexColor('#cbd5df')),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,8)])
def page():story.append(PageBreak())
add('骨肉瘤 CAF 批次与排除线统一处理报告','title');add('2026-09-19 | phase10_batch_gates_20260919 | 完成准入统一、补充环境RNA检查及34组模型比较')
add('统一规则保留3,976个严格候选，原始结果完整保留','head')
add('对6项研究、40条样本记录中的255,241个既有QC后细胞重新计算统一证据。严格主集包含5项研究、14条样本记录、3,976个候选细胞：原4,300个联合候选保留3,937个，另有39个通过统一规则进入；363个原候选未进入新严格集，其证据与排除原因仍保留。')
table(['研究','原联合候选','统一严格集','成骨敏感集','竞争程序敏感集'],[[r.dataset,int(r.old_joint_cells),int(r.unified_strict),int(r.osteo_sensitive),int(r.competitor_sensitive)] for r in changes.itertuples()],[112,83,90,95,107])
add('最早7个样本已补做DecontX；合计36条样本记录有可用去环境RNA结果，其余4条在既有FibLike范围内没有候选并保留未评估标记。严格主集每个细胞均有原始及去环境RNA证据、双CNV低信号与参考检查记录，并通过已有独立双细胞排除保护。')
add('Methods：各来源使用同一判定代码','head')
add('保留既有宽谱系注释的FibLike范围，统一重新计算相同Fib、竞争谱系、成骨及CTA面板。严格规则为检出基因≥500、原始和DecontX后均至少2个Fib基因计数≥1、通过双CNV低信号与正常参考检查、群级任一/双方CNV升高比例均<20%，并通过竞争、混合、CTA及独立双细胞保护。成骨冲突定义为≥3个程序基因且占完整文库计数>0.2%；完整面板和全部判定代码随包保存。')
add('各来源原始完整文库作为归一化分母；用15,391个共同测量Ensembl基因构建联合对象。主模型使用log1p(CP10k)、按研究×技术选择2,000个HVG、PC30、Harmony theta=2及动态lambda、15近邻、Leiden分辨率0.5。表达与marker描述来自未改动原始counts，批次修正作用于低维坐标。')
page()
add('排除线已统一，敏感细胞单列保存','head')
add('原GSE152048减少226个候选，最早GSE270231中137个原候选未入主集，同时39个先前受不同marker规则影响的细胞通过统一规则进入。后续GEO的原2,065个主候选保持在严格集中。旧批次的群级CNV风险、混合/竞争程序和去环境RNA支持是主要变化来源；排除原因可重叠，不能简单相加。')
table(['比较集合／规则','候选数','相对严格集','用途'],[
['统一严格集','3976','基线','当前工作主集'],
['仅成骨比例改为0.5%','4193','+217','成骨相关准入敏感性'],
['竞争/混合/成骨比例改为0.5%','4527','+551','宽松规则敏感性；不自动晋升'],
['群级CNV风险上限改为30%','4067','+91','排除线边界敏感性'],
['校正后检出由≥1改为>0','4193','+217','微量残留计数敏感性'],
['检出基因下限改为300','3994','+18','低复杂度诊断'],
['Fib基因门槛改为≥3','3303','-673','阳性证据门槛诊断']],[188,64,76,159])
add('两个4,193集合只是细胞总数相同，规则和成员不应视为相同。另保留1,769个满足所定义证据要求的Mural参考候选，供后续MYH11/血管壁相关程序对照；该集合与CAF敏感性集合可能重叠，不计入严格主集数量。')
add('严格线0.2%暂作为保守工作规则，0.5%作为并列敏感性。所有阈值均记录适用条件和判定原因，不把单条程序线等同于肿瘤身份真值，也不为增加marker阳性数量而修改阈值。')
add('阴性参考检查支持保留≥1的校正计数工作门槛','head')
add('在3,109个已登记的留出免疫／内皮参考细胞中，去掉FibLike范围限制后，严格、成骨0.5%和竞争0.5%的完整分子/质量规则均无细胞通过；校正计数改为>0时有2个通过。这是内部参考中的描述性泄漏检查，不等于独立证明CAF分类准确率。')
add('所有原候选的归属、原群标签、CNV结果、原始和校正面板证据及剔除原因都保留在逐细胞表中。无CNV结果或无环境RNA证据的细胞未被当作阴性结果放行。')
page()
add('批次分离减弱，技术与供者依赖仍可见','head')
p=metrics.loc['unified_primary'];u=metrics.loc['unified_uncorrected'];add(f'在同一个3,976细胞集合内，同批次近邻比例由{u.same_batch_neighbors:.2%}降至{p.same_batch_neighbors:.2%}，批次轮廓系数由{u.batch_silhouette:.3f}降至{p.batch_silhouette:.3f}。主分析得到11个计算分群。主模型是预先保留的工作基线，没有以混合度最大或群数贴合论文作为选参目标。')
rows=[]
for name,title in [('unified_primary','主分析'),('theta1','theta=1'),('sample_batch','按样本校正'),('fixed_lambda1','固定lambda=1'),('balanced_PCA_fit','供者记录均衡PCA拟合'),('seed_sensitivity','另一随机种子'),('no_FLEX','去除FLEX'),('corrected_counts','使用DecontX计数'),('osteo_sensitive','成骨阈值0.5%'),('competitor_sensitive','竞争阈值0.5%'),('CNV_risk_0_3','群CNV风险线30%')]:
 r=metrics.loc[name];rows.append([title,int(r.cells),int(r.clusters),'-' if pd.isna(r.ARI_vs_unified_primary) else f'{r.ARI_vs_unified_primary:.3f}'])
table(['条件','细胞数','群数','共同细胞ARI'],rows,[223,80,65,119])
h=metrics.loc[metrics.index.str.startswith('hvg')];add(f'另对HVG1500/2000/3000与PC20/30/40进行了完整组合检查，工作主模型以外8组的ARI范围为{h.ARI_vs_unified_primary.min():.3f}-{h.ARI_vs_unified_primary.max():.3f}。14次供者记录留出重拟合的ARI范围为{loo.ARI.min():.3f}-{loo.ARI.max():.3f}。两项结果均为参数或记录依赖诊断，不是独立外部验证。')
add('成骨阈值0.5%时与主集较一致（ARI=0.918）；去除FLEX后的变化较大（ARI=0.534）。按样本校正使分群明显改变（ARI=0.279），不能只凭混合更充分便改作主结果。保留共同分群、技术分层和敏感性对照，用于区分可重复状态与依赖来源的结构。')
page()
add('交付的是统一分析基础，亚型命名继续受证据约束','head')
dominated=support.loc[support.interpretation=='donor_dominated','cluster'].astype(str).tolist();add('主分群用U编号保存。按最大已映射单供者占比>70%的复核规则，需重点检查的群为：'+', '.join('U'+c for c in dominated)+'。完整供者细胞数、跨研究支持及逐群留出Jaccard均已交付；70%仍是风险提示，不是有效亚型的通过线。')
add('逐群留出Jaccard仅统计移除供者后仍有至少20个该群细胞的轮次；附图列出有效轮次数。因此，较高的条件稳定性不能代替多供者支持。')
add('Discussion：宽谱系注释范围、来源特定的CNV检测器和参考定义保留追踪记录；本轮统一的是证据计算、准入决策和批次处理，并未获得独立的CAF身份真值。DecontX使用样本内已有宽谱系标签、无空液滴背景，其结果作为支持证据。FLEX与全转录组测量差异、患者/部位混杂及跨队列身份重叠仍影响细分群的外推。')
add('14条主分析样本记录不能自动等同于14名跨队列独立患者。GSE250015跨队列重叠未解决，仍不进入患者推断统计；未知患者身份不用于增加独立样本数。本轮没有新增疗效、肺转移或亚型临床价值的推断性结论。')
add('SR1004按用户原名称保留。本轮不以该marker的检出结果选参数，也不将它改写为S100A4；后续精确marker核查应同时检查名称、原始表达和各平台是否实际测量。')
add('交付内容与验收','head')
add('objects/：统一严格集、共同特征候选联合对象，以及去除FLEX、DecontX、成骨/竞争/CNV风险敏感性的独立对象。tables/：255,241细胞的完整准入证据、集合变化、全来源阈值扫描、34组模型比较、逐群/供者支持、旧新群对应及描述性marker。原始完整计数与额外来源级对象留在服务器，交付包避免重复复制输入。')
add('验收通过：最早7个DecontX任务完成；7个主要执行阶段退出码均为0；40个完整来源文件SHA256与处理时记录一致；3,976个严格细胞原始共同基因counts与来源候选对象逐值一致；主集无已复核排除的双细胞；上一版4,300细胞联合对象保持原样。')
add('处理流程：complete_early_ambient.py → unify_gates.py → apply_quality_guard.py → batch_models.py；all_source_thresholds.py与cnv_risk_model.py补充全来源阈值及CNV边界检查；review_batch_gates.py完成必要验收。PROTOCOL.json及resources/batch_protocol.json记录主规则、默认值与审计后增加的敏感性。')
add('图件附后：统一集合整合前后UMAP、各队列准入数量、批次模型比较、分群来源与稳定性。所有图件只描述已计算结果，未以视觉混合代替生物学验证。')
def footer(c,d):c.setFont('STSong-Light',8);c.setFillColor(colors.HexColor('#657180'));c.drawString(54,26,'CAF 批次与排除线统一处理 | 2026-09-19');c.drawRightString(A4[0]-54,26,str(d.page))
body=Q/'body.pdf';SimpleDocTemplate(str(body),pagesize=A4,rightMargin=54,leftMargin=54,topMargin=42,bottomMargin=43,title='骨肉瘤CAF批次与排除线统一处理').build(story,onFirstPage=footer,onLaterPages=footer)
figs=[('01_unified_umap.png','统一严格集：整合前后UMAP'),('02_gate_changes.png','原候选与统一准入集合'),('03_batch_sensitivity.png','批次处理及排除线敏感性'),('04_cluster_support.png','分群来源与重拟合稳定性')]
f=Q/'figures.pdf';c=canvas.Canvas(str(f),pagesize=landscape(A3));w,hh=landscape(A3)
for i,(name,title) in enumerate(figs,1):
 c.setFont('STSong-Light',16);c.drawString(35,hh-35,'附图 A'+str(i)+' '+title);im=Image.open(O/'figures'/name);iw,ih=im.size;scale=min((w-70)/iw,(hh-100)/ih);c.drawImage(str(O/'figures'/name),(w-iw*scale)/2,(hh-ih*scale)/2-12,iw*scale,ih*scale);c.showPage();md.append(f'![{title}](../figures/{name})\n')
c.save();writer=PdfWriter();writer.append(PdfReader(body));writer.append(PdfReader(f));out=D/'骨肉瘤CAF_批次与排除线统一处理_20260919.pdf';writer.write(out);(D/'batch_gate_remediation_report.md').write_text('\n'.join(md),encoding='utf-8');pdf=pdfium.PdfDocument(str(out))
for i in range(len(pdf)):pdf[i].render(scale=1.35).to_pil().save(Q/f'page-{i+1}.png')
text='\n'.join(p.extract_text() for p in PdfReader(out).pages);assert all(x in text for x in ['3,976','255,241','SR1004','S100A4','0.918'])
print(json.dumps({'pdf_pages':len(pdf),'body_pages':len(PdfReader(body).pages),'files_verified':len(m),'text_checks_passed':True}))