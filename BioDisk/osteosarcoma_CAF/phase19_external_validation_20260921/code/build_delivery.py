from pathlib import Path
import json,shutil,hashlib,zipfile,datetime,html,re
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
B=Path('C:/codex_test');O=B/'CAF_external_delivery_20260921';D=O/'delivery';D.mkdir(exist_ok=True)
F=pd.read_csv(O/'tables/external_primary_associations.tsv',sep='\t');S=pd.read_csv(O/'tables/external_sensitivity_associations.tsv',sep='\t');QA=json.loads((D/'VALIDATION_QA.json').read_text())
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('CNH',r'C:\Windows\Fonts\simhei.ttf'))
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='CNBody',fontName='CN',fontSize=11,leading=17,spaceAfter=7));styles.add(ParagraphStyle(name='CNTitle',fontName='CNH',fontSize=20,leading=26,spaceAfter=12));styles.add(ParagraphStyle(name='CNHead',fontName='CNH',fontSize=13,leading=19,spaceBefore=10,spaceAfter=7));styles.add(ParagraphStyle(name='CNCell',fontName='CN',fontSize=9.5,leading=14));styles.add(ParagraphStyle(name='CNNote',fontName='CN',fontSize=9,leading=14,spaceAfter=6,textColor=colors.HexColor('#444444')))
story=[];md=[]
def formatted(t):return re.sub(r'[A-Za-z][A-Za-z0-9_]*(?:[./][A-Za-z0-9_]+)*',lambda m:'<nobr>'+m.group(0)+'</nobr>',html.escape(t))
def para(t,style='CNBody'):story.append(Paragraph(formatted(t),styles[style]));md.append(t+'\n')
def head(t):story.append(Paragraph(html.escape(t),styles['CNHead']));md.append('## '+t+'\n')
def table(headers,rows,widths):
 data=[[Paragraph(formatted(str(x)),styles['CNCell']) for x in row] for row in [headers]+rows];tb=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
 tb.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eef2')),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#708390')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#d4dadd')),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]));story.append(tb);story.append(Spacer(1,8));md.append('| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(str(x) for x in row)+' |' for row in rows)+'\n')
def newpage():story.append(PageBreak())
para('骨肉瘤CAF阶段结果与外部验证','CNTitle');para('阶段性交付  2026年9月21日','CNNote')
para('本轮已完成全队列分层探索、逐群证据登记及两个独立组织队列的候选特征验证。C7和C9保留为有跨供者表达支持的候选线索；当前外部数据不支持其与转移或复发结局存在明确关联。最终“转移相关新CAF亚群及可用流式分选panel”尚未确立，本次交付用于客户审阅阶段结果和确定下一步。')
head('全队列探索形成两个优先复核候选')
table(['项目','已完成结果'],[['数据范围','52条样本记录，8个来源，299,483个QC后细胞'],['分层分析','严格层4,099个；支持层54,997个；广义拓展层140,135个候选细胞'],['逐群证据','三层共52个表达群、260项程序证据；19项达到跨层探索性支持，13项达到跨供者重复支持'],['支持层C7','2,509个细胞，2,240个原始面板myCAF阳性；myCAF程序19条可评估已知供者记录中18条同向'],['支持层C9','3,754个细胞，2,686个原始面板myCAF阳性；矿化和iCAF程序有跨层探索性支持']],[86,425])
head('外部结局验证未获得关联支持')
rows=[]
for r in F.to_dict('records'):
 endpoint='转移标签 OR' if r['cohort']=='GSE21257' else '复发随访 HR'
 rows.append([r['cohort']+'\n'+str(r['n'])+'例',r['signature'],endpoint,f'{r["effect"]:.2f} ({r["CI_low"]:.2f}-{r["CI_high"]:.2f})',f'{r["P"]:.3f}',f'{r["FDR"]:.3f}'])
table(['队列','候选','终点与效应','效应值及95%CI','P','FDR'],rows,[90,37,86,162,68,68])
para('效应对应特征分数每升高1个队列标准差。GSE21257为34例有转移标签与19例无转移标签；GSE39055为36例有有效复发随访时间，其中18例发生复发。各终点的两个候选分别进行BH校正。','CNNote')
newpage();head('验证结果在已检查的敏感性分析中保持阴性')
story.append(Image(str(O/'figures/external_validation_forest.png'),width=505,height=212));md.append('![外部验证](../figures/external_validation_forest.png)\n')
para('两个队列中，C7和C9的95%置信区间均跨越无效应值1，P值及FDR均未达到0.05。C9在转移标签队列中的效应方向与复发随访队列不同；两者终点不同，不据此合并为一个转移效应。')
table(['GSE21257检查','C7 OR 及P值','C9 OR 及P值'],[[label,*(f'{S[(S.model==model)&(S.signature==sig)].iloc[0].effect:.2f}; P={S[(S.model==model)&(S.signature==sig)].iloc[0].P:.3f}' for sig in ['C7','C9'])] for model,label in [('exclude_baseline_metastasis','排除诊断时已转移14例，n=39'),('adjust_age_sex','调整年龄与性别，n=53'),('rank_score_sensitivity','改用基因百分位评分，n=53')]],[245,133,133])
head('组织表达共现仍提供程序线索')
para('两个候选在两项验证平台均达到100%基因覆盖。GSE21257中C7的14个、C9的15个基因与各自其余特征基因的平均表达均呈正相关；GSE39055分别为11/14和13/15。该结果显示这些基因在独立组织数据中仍有共表达结构，但不能替代临床结局关联，也不能确认其来自CAF细胞。')
head('现有配对和标签结果完整保留')
para('支持层仅1组患者满足骨、肺两侧各至少10个myCAF的集合比较条件；广义拓展层为2组。支持层C7在可比较患者的骨、肺两侧各仅3个myCAF，C9两侧均无该群myCAF，尚未形成候选群的配对程序验证。')
para('GSE162454的任何转移标签比较中，支持层C9为Yes 2例、No 3例，myCAF中位差+0.252、成骨+0.183、矿化+0.004，iCAF −0.176、apCAF −0.255。该小样本结果作为描述保留，不冒充确认的肺转移A/B1比较。')
newpage();head('客户原始目标与本次交付状态')
accept=[['A B1 B2及治疗手术标注','部分完成','52条记录主表已交付；缺失的肺转移时点、手术与取样关系继续标为未知'],['CAF与myCAF分群及三层探索','完成本轮计算','52个群、候选基因、逐细胞分层和证据表可复核；群名为操作性状态'],['成骨与矿化重点分析','完成本轮计算','客户完整面板、效应方向、供者一致性及外部候选验证均保留'],['同患者骨肺分析','完成可评估部分','4组配对线索均登记，细胞或临床条件不足者未强制比较'],['独立外部验证','组织层验证完成','53例转移标签与36例复发随访未支持明确关联；独立单细胞身份验证未完成'],['表面marker与分选panel','候选筛查完成','单/双marker共136种RNA组合及既有蛋白定位、抗体资料已提供；无最终可用panel'],['转移相关新CAF亚群','尚未达到最终目标','C7 C9为候选线索，不能写成已验证的新亚型']]
table(['客户目标','状态','本次交付及未完成项'],accept,[121,96,294]);pd.DataFrame(accept,columns=['客户目标','状态','本次交付及未完成项']).to_csv(O/'tables/customer_acceptance.tsv',sep='\t',index=False)
head('候选分选思路可反馈但不能直接作为最终方案')
para('FAP+THY1可保留为后续复核组合。支持层C7在16条有足量目标细胞的供者/部位记录中，RNA目标覆盖中位数31.2%，门内已知目标比例3.6%；C9在18条记录中分别为46.5%和11.6%。这两组数值不足以支持直接用于高纯度分选。')
para('候选门控思路为活细胞与单细胞门后，评估CD45、CD31、EPCAM排除门及FAP、THY1候选富集；同时保留血管壁与恶性间质细胞误入评估。最终抗体、荧光通道及阈值应结合消化条件、仪器和蛋白实测确定。')
para('建议先请客户确认是否接受C7和C9作为下一步探索对象，以及优先补充临床资料还是开展候选蛋白验证。本次材料适合作为阶段结果反馈，不建议按最终目标全部达成进行验收。')
newpage();head('分析方法')
para('候选特征在查看本轮外部结果前冻结：取支持层C7和C9中具有跨供者正向支持的描述性marker，按群内减群外中位表达效应排序，分别固定14和15个基因。基因名单和冻结时间见FROZEN_PROTOCOL.json及frozen_signatures.tsv。')
para('GSE21257采用已验收的作者标准化表达，按基因取探针中位数。GSE39055采用官方series matrix与GPL14951平台注释；完整验收29,377行平台记录，探针按基因取中位数。两队列保留作者的方差稳定化及标准化尺度。')
para('每个队列中先按基因跨样本标准化，再对固定基因集取均值并标准化为特征分数。GSE21257采用连续分数logistic回归，报告每标准差的OR及95%CI；另行检查基线转移排除、年龄性别调整和百分位评分。GSE39055采用Efron并列事件处理的Cox模型，报告复发结局HR及95%CI；1例随访时间为0者排除。')
para('多重检验按各队列、各终点的两个候选组成同类假设家族进行BH校正；敏感性分析按模型分别校正。统计样本量以患者为单位，单细胞分析中的细胞数不代替生物学重复。')
head('解释边界与未完成验证')
para('GSE21257的终点包含诊断时已有转移及随后转移，不等同于肺转移特异终点；GSE39055提供复发或末次随访时间，不能将其死亡状态套用该时间进行总生存分析。两队列为组织bulk表达，不能证明CAF比例、CAF身份或表面蛋白分选纯度。')
para('GSE21257此前曾用于其他固定通用程序的分析，本轮属于新候选的回顾性外部队列检验；GSE39055未参与候选发现。当前阴性结果不否定单细胞中已观察到的表达状态，但不支持将C7/C9升级为临床结局相关、已验证的新CAF亚群。')
para('分层中的L3代表内部跨层探索性支持，并非独立外部验证。患者身份未确认的样本按样本展示，不计入跨供者重复支持；跨来源供者记录不直接合并为独立患者数。广义拓展层含CNV冲突与身份未定细胞，各群风险已单列。')
head('交付内容与来源')
para('交付包包含本报告及可编辑Markdown、52样本临床表、冻结候选基因、全部外部效应与P/FDR、52群证据卡、完整程序/基因/表面组合表、图件和复现脚本。大体积表达对象与原始阶段结果保留在项目服务器。','CNNote')
para('官方来源：GSE21257 / PRJNA126279，https://www.ncbi.nlm.nih.gov/bioproject/PRJNA126279；GSE39055，https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39055；GPL14951，https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951。','CNNote')
def footer(can,doc):
 can.saveState();can.setFont('CN',8);can.setFillColor(colors.HexColor('#666666'));can.drawString(42,25,'骨肉瘤CAF阶段结果  2026年9月21日');can.drawRightString(A4[0]-42,25,str(doc.page));can.restoreState()
pdf=D/'骨肉瘤CAF阶段交付报告.pdf';doc=SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=42,leftMargin=42,topMargin=36,bottomMargin=42,title='骨肉瘤CAF阶段结果与外部验证',author='')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
(D/'骨肉瘤CAF阶段交付报告.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
message='您好，本次骨肉瘤CAF阶段分析结果已整理完成，包含全队列CAF/myCAF探索、成骨与矿化程序、分层证据及表面marker候选。\n\n我们补充检查了两个独立组织队列：当前C7/C9候选具有表达程序线索，但外部转移或复发结局尚未得到显著支持，现阶段也未形成可直接使用的高纯度分选panel。因此这版作为阶段结果先请您审阅，候选群、阴性结果及需要补充的信息都已列明。\n\n您可重点看C7和C9的候选方向，以及是否优先补临床分组资料或先做候选表面蛋白验证。我们再根据您的反馈确定后续重点。'
(D/'客户反馈消息草稿.txt').write_text(message,encoding='utf-8')
bundle=O/'client_package';bundle.mkdir(exist_ok=True)
for folder in ['report','external_validation','evidence_register','figures','methods','sample_metadata']:(bundle/folder).mkdir(exist_ok=True)
for f in D.glob('*'):
 if f.suffix in ['.pdf','.md'] and f.name.startswith('骨肉瘤CAF阶段交付报告'):shutil.copy2(f,bundle/'report'/f.name)
for f in (O/'tables').glob('*'):shutil.copy2(f,bundle/'external_validation'/f.name)
for f in (B/'CAF_evidence_20260921/tables').glob('*'):
 if 'analysis_labels' not in f.name and 'cluster_surface_donor_screen' not in f.name:shutil.copy2(f,bundle/'evidence_register'/f.name)
for f in (O/'figures').glob('*.png'):shutil.copy2(f,bundle/'figures'/f.name)
for f in [B/'CAF_expanded_20260921/figures/supported_joint_umap.png',B/'CAF_evidence_20260921/figures/program_evidence_levels.png']:shutil.copy2(f,bundle/'figures'/f.name)
for f in [O/'FROZEN_PROTOCOL.json',D/'VALIDATION_QA.json',B/'CAF_evidence_20260921/PROTOCOL.json']:shutil.copy2(f,bundle/'methods'/f.name)
for f in (O/'code').glob('*.py'):shutil.copy2(f,bundle/'methods'/f.name)
for f in (B/'CAF_evidence_20260921/code').glob('*.py'):shutil.copy2(f,bundle/'methods'/('phase18_'+f.name))
shutil.copy2(B/'CAF_client_steps_20260920/tables/01_clinical_sample_audit_updated.tsv',bundle/'sample_metadata/原40样本临床证据.tsv');shutil.copy2(O/'inputs/ScPCA_12_manifest.tsv',bundle/'sample_metadata/ScPCA原始元数据.tsv')
(bundle/'README.md').write_text('# 骨肉瘤CAF阶段性交付包\n\n先阅读report/骨肉瘤CAF阶段交付报告.pdf。外部验证已实施，当前未支持C7/C9与临床结局的明确关联；本包不代表最终分选panel已经验证。\n\nexternal_validation包含外部统计及52样本主表；evidence_register保留正向、阴性和不可评估结果；methods包含冻结方案与复现脚本。完整表达矩阵、逐细胞大表及门控明细保留于项目服务器的Phase17至Phase19目录，可按需提供。统计推断以患者为单位。\n',encoding='utf-8')
manifest={str(p.relative_to(bundle)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in bundle.rglob('*') if p.is_file() and p.name!='SHA256.json'}
(bundle/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zipfile_path=D/'CAF阶段性交付包_20260921.zip'
with zipfile.ZipFile(zipfile_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(bundle.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(bundle))
print('PDF',pdf,'PACKAGE',zipfile_path,'FILES',len(manifest),flush=True)
