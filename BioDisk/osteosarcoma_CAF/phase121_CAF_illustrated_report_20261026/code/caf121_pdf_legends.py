from pathlib import Path
import json,hashlib,math,datetime,shutil
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4,landscape
from reportlab.lib.units import mm
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.lib import colors
from PIL import Image
import fitz
W=Path('C:/Users/31937/caf121_report_work');D=Path('C:/Users/31937/Downloads/CAF_illustrated_report_20261026');D.mkdir(parents=True,exist_ok=True)
F=W/'figures';meta=json.loads((W/'FIGURES_COMPLETE.json').read_text());sig=meta['retained_signal_nodes'];decision=pd.read_csv(W/'resources/node_decisions_ZH.tsv',sep='\t').set_index('父群或子群');stats=pd.read_csv(W/'resources/focus_functional_audit.tsv',sep='\t');rv=pd.read_csv(W/'resources/target_collagen_candidate_review.tsv',sep='\t').set_index('node_id');lit=json.loads((W/'resources/LITERATURE_CROSSCHECK.json').read_text(encoding='utf-8'))
# Present every retained signal group on the same panel; no expression/statistical changes.
pr=pd.read_csv(W/'tables/marker_dotplot_values.tsv',sep='\t');genes=['DCN','LUM','PDGFRA','COL1A1','COL1A2','COL3A1','COL5A1','POSTN','MMP11','FNDC1','COMP','COL11A1','SERPINH1','PCOLCE','LOX','LOXL2','ACTA2','TAGLN','MYH11','RGS5','CD36','C7','RUNX2','ALPL','IBSP','PTPRC','PECAM1']
a=pr.pivot(index='node_id',columns='gene',values='mean').reindex(index=sig,columns=genes);b=pr.pivot(index='node_id',columns='gene',values='detection').reindex(index=sig,columns=genes);z=(a-a.mean())/a.std(ddof=0).replace(0,1)
plt.rcParams.update({'font.size':11,'pdf.fonttype':42});fig,ax=plt.subplots(figsize=(13,7));xx,yy=np.meshgrid(range(len(genes)),range(len(sig)));im=ax.scatter(xx.ravel(),yy.ravel(),s=b.to_numpy().ravel()*130,c=z.to_numpy().ravel(),cmap='RdBu_r',norm=Normalize(-2,2));ax.set_xticks(range(len(genes)),genes,rotation=65,ha='right',fontsize=10);ax.set_yticks(range(len(sig)),sig,fontsize=10);ax.set_ylim(len(sig)-.4,-.7);ax.set_axisbelow(True);ax.grid(alpha=.15);ax.set_title('All retained signal groups: predefined marker panel',pad=42)
handles=[ax.scatter([],[],s=v*130,c='gray') for v in [.1,.5,1]];ax.legend(handles,['10%','50%','100%'],title='RNA detection',ncol=3,loc='lower left',bbox_to_anchor=(0,1.01),frameon=False,fontsize=9);fig.colorbar(im,ax=ax,label='Gene-wise z-score of mean log1p(CP10k)',shrink=.7);fig.tight_layout();fig.savefig(F/'04b_all_signal_markers.png',dpi=300,bbox_inches='tight');fig.savefig(F/'04b_all_signal_markers.pdf',bbox_inches='tight');plt.close(fig)
pdfmetrics.registerFont(TTFont('YH','C:/Windows/Fonts/msyh.ttc'));pdfmetrics.registerFont(TTFont('YHB','C:/Windows/Fonts/msyhbd.ttc'))
pdfmetrics.registerFontFamily('YH',normal='YH',bold='YHB',italic='YH',boldItalic='YHB')
pw,ph=landscape(A4);left=16*mm;right=pw-left;usable=right-left;out=D/'骨肺CAF亚群_图文分析与Lung04推荐报告_图例补充版.pdf'
if out.exists() and not (W/'report_before_dotplot_weighting.pdf').exists():shutil.copyfile(out,W/'report_before_dotplot_weighting.pdf')
c=canvas.Canvas(str(out),pagesize=(pw,ph));c.setTitle('骨肺CAF亚群：基质相关信号与Lung_04优先推荐');c.setAuthor('CAF研究分析');page=0
style=ParagraphStyle('body',fontName='YH',fontSize=11,leading=17,textColor=colors.HexColor('#263b48'),wordWrap='CJK',spaceAfter=6)
small=ParagraphStyle('small',parent=style,fontSize=9,leading=13)
def para(text,y,size='body',width=None,x=None):
 p=Paragraph(text,style if size=='body' else small);w,h=p.wrap(width or usable,1000)
 assert y-h>=18*mm,(page,'paragraph overflow',y,h)
 p.drawOn(c,left if x is None else x,y-h);return y-h-6
def start(title,kicker='骨肺CAF | 固定6,633细胞队列'):
 global page
 if page:c.showPage()
 page+=1;c.setFillColor(colors.HexColor('#1c6470'));c.rect(0,ph-7*mm,pw,7*mm,fill=1,stroke=0);c.setFont('YH',9);c.setFillColor(colors.HexColor('#617580'));c.drawString(left,ph-14*mm,kicker);c.setFillColor(colors.HexColor('#173b4b'));c.setFont('YHB',21);c.drawString(left,ph-26*mm,title);c.setStrokeColor(colors.HexColor('#d1dfe3'));c.line(left,15*mm,right,15*mm);c.setFont('YH',8);c.setFillColor(colors.HexColor('#617580'));c.drawString(left,10*mm,'2026-10-26  |  Phase121  |  保留既有分群与统计结果');c.drawRightString(right,10*mm,str(page));c.bookmarkPage('p'+str(page));c.addOutlineEntry(title,'p'+str(page),level=0);return ph-34*mm
def table(data,y,widths,fs=10):
 st=ParagraphStyle('cell',parent=small,fontSize=fs,leading=fs*1.45)
 data=[[Paragraph(str(v),st) for v in row] for row in data];t=Table(data,colWidths=[v*mm for v in widths],repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5f0f2')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,-1),.3,colors.HexColor('#d6e2e6'))]));w,h=t.wrap(usable,1000);assert y-h>18*mm,(page,'table overflow',h);t.drawOn(c,left,y-h);return y-h-10
def figure(title,file,lead,caption):
 y=start(title);y=para(lead,y,'small');p=F/(file+'.png');iw,ih=Image.open(p).size
 cap=Paragraph(caption,small);_,ch=cap.wrap(usable,1000);bottom=22*mm+ch;maxh=y-bottom-5*mm;scale=min(usable/iw,maxh/ih);w,h=iw*scale,ih*scale;c.drawImage(str(p),left+(usable-w)/2,y-h,width=w,height=h,mask='auto');cap.drawOn(c,left,22*mm)
y=start('骨肺CAF亚群：基质相关信号与优先候选','研究报告 | UMAP · marker气泡图 · 患者支持 · 文献交叉')
y=para('<b>优先推荐 Lung_04，保留其他有信号群。</b> Lung_04在主分析中有7位有效患者支持，基质、胶原交联及FNDC1文献面板呈一致增强，是当前最适合推进后续验证的父群。',y)
y=table([['重点群','细胞数','主分析有效患者*','结果定位'],['Lung_04','233','7','优先推荐：基质/胶原相关，伴成骨程序'],['Bone_07','131','4','保留：骨侧基质与收缩相关线索'],['Bone_07_s1','35','3','保留：同父群内胶原相关方向线索'],['Lung_04_s2','36','1','保留：同父群内基质偏向线索'] ],y,[43,25,55,142],11)
y=para('*主分析有效患者：同来源、同病灶、同一患者的目标群及对照各至少5个细胞。父群与子群的细胞数不能相加。',y,'small')
y=para('本报告共保留13个有信号群。患者数量用于确定研究优先级，严格身份分析用于评估结果稳健性，不作为一票否决条件。Lung_04的严格身份层有5位有效患者，基质方向仍一致。',y)
y=para('阅读顺序：队列和推荐规则 → UMAP与marker → 患者构成和功能效应 → 文献对应 → 所有信号群 → 最终推荐。',y,'small')
y=start('队列与患者支持：人数决定证据权重')
y=table([['病灶组','来源','完整细胞','主分析细胞','来源内患者记录'],['骨','SCPCP000017 / SCPCP000023','2,788','2,788','12 + 4'],['肺','SCPCP000017 / SCPCP000023 / NCI','3,845','3,635','3 + 3 + 17'],['合计','固定原有队列','6,633','6,423','不同来源不直接视为全球去重人数']],y,[26,99,36,38,66],10)
y=para('210个已有软件双细胞阳性标记细胞保留在完整名单中，不进入本次主分析UMAP与气泡图。没有加入NCC或正常参考细胞。沿用骨组15个父群、肺组14个父群；仅生成可视化，不改变聚类标签。',y)
y=para('<b>推荐规则：</b>优先考虑至少6位有效患者、每位目标与对照各≥5个细胞、至少75%同向，并有目标程序原始P&lt;0.05的群；10位及以上可增加证据支持。每患者≥10细胞、身份敏感性和留一患者分析用于补充评估。',y)
y=para('6位是本项目的操作性下限，不是普适文献标准或充分功效保证。现有双侧精确检验在5位时最小P为0.0625，6位时为0.03125。文献强调独立患者重复，细胞数量不能替代患者数量[8,9]。',y,'small')
y=para('<b>本次口径更新：</b>按用户最新要求，用主分析患者支持确定优先推荐；Phase120较严格的身份层人数条件改为敏感性信息。历史记录保留，表达值、P值、FDR及分群不变。',y,'small')
figure('UMAP：骨、肺父群结构','01_parent_umap','骨组2,788个、肺组3,635个主分析细胞；图中数字对应父群编号后两位。','图1｜分别在保存的PCA/Harmony表示上计算UMAP，固定随机种子20261026，15近邻、min_dist=0.3。骨、肺嵌入独立计算，坐标及距离不能跨图直接比较；分离不等同于功能特异或患者重复。')
figure('UMAP：来源构成与重点父子群','02_sources_candidates','来源分布和重点候选显示在相同的病灶内坐标上。','图2｜左列按来源着色；右列灰色为其他父群，青色为目标父群中未单独着色的细胞，红/黄为指定子群。图例颜色类别互斥：Bone_07=66+35+30=131；Lung_04=197+36=233。来源集中不作为自动剔除依据。')
dotcap='气泡大小：RNA检出比例；颜色：逐基因在图示群间标准化的平均log1p(CP10k)。先在患者内平均，再对来源内患者记录等权汇总，避免少数细胞来源权重过大。预设功能/身份marker面板，不是新做的全基因差异排名；红色不代表统计显著，跨图颜色不能直接比较。'
figure('Marker气泡图：骨组全部父群','03_markers_Bone','并列展示成纤维、基质、胶原加工、收缩、成骨和竞争谱系标记。','图3｜'+dotcap)
figure('Marker气泡图：肺组全部父群','03_markers_Lung','Lung_04的基质相关表达应结合患者效应与成骨重叠背景解释。','图4｜'+dotcap)
figure('Marker气泡图：全部13个有信号群','04b_all_signal_markers','保留所有既有胶原候选及明确的基质表达线索，不因患者支持较少删除。','图5｜'+dotcap+' Lung_05_s2/s4等竞争谱系标记需结合身份记录核对；表达线索不直接证明纯CAF身份。')
figure('骨组：父群的逐患者细胞构成','05_patient_Bone','多个大父群主要由单一患者贡献；不能用细胞多替代跨患者支持。','图6｜每行一个父群、每列一位来源内患者；数字为细胞数，颜色为该患者占该父群细胞的比例。患者ID在来源内解释，不将不同来源相似ID自动合并。')
figure('肺组：父群的逐患者细胞构成','05_patient_Lung','Lung_04的患者构成较分散；其余有基质表达的群仍保留为探索线索。','图7｜读图方式同图6。该图是细胞贡献比例，不是患者患病率、亚群丰度组间检验或骨肺差异分析。')
figure('有效患者支持：Lung_04优先','06_patient_gate','由“出现过细胞”到“可稳定配对”有明显差别；红线为主分析建议下限6位。','图8｜灰色：两侧至少1细胞；蓝色：两侧至少5细胞；橙色：严格身份层且两侧至少5细胞。均为表中选定来源内有效配对患者。橙色用于稳健性说明，不作为推荐的一票否决。')
figure('功能效应：逐患者方向与均值','07_effects','主分析采用每侧至少5细胞的患者配对；每个小点为一位患者，菱形为均值。','图9｜父群与同病灶其余候选细胞比较，子群与同父群其他细胞比较；不同层级对照不同，不直接按效应大小排名。效应为log1p(CP10k)程序均分之差，不是倍数变化；无有效配对时没有点。')
figure('文献交叉：ECM与收缩、脂质、炎症分开解释','08_literature','Lung_04在FNDC1论文面板中增强；Bone_07更偏收缩/血管相关参考方向。','图10｜采用原文FNDC1、MYH11、CD36、C7基因面板[1]，显示min5患者均值差。Lung_04的FNDC1面板7/7同向，P=0.015625，FDR=0.06875；不能据此称原文亚群已被完整复现。NA为无有效配对。')
figure('Lung_04敏感性：配对减少主要来自对照不足','09_lung04_attrition','严格身份层从7位降到5位，剩余患者的基质方向仍一致。','图11｜两位患者的目标细胞仍较充足：6546为34→24，7978为14→13；但合格对照分别16→4、5→3，低于5细胞条件。这里是配对可评估性变化，不能解释为基质效应消失。')
y=start('Lung_04：核心效应与敏感性数值')
rows=[['分析层 / 程序','N','同向','平均效应','原始P','FDR']]
for mode,title in [('main_min5','主分析'),('identity_min5','身份敏感性')]:
 for ep,label in [('program__mCAF_matrix','基质'),('program__collagen_crosslinking','交联'),('audit__FNDC1_paper','FNDC1'),('audit__collagen_processing_nonCOL','非COL加工')]:
  q=stats[stats.node_id.eq('Lung_04')&stats.source.eq('NCI')&stats['mode'].eq(mode)&stats.endpoint.eq(ep)].iloc[0];rows.append([title+' / '+label,int(q.N),f'{int(q.positive_patients)}/{int(q.N)}',f'{q.mean_effect:+.3f}',f'{q.P:.5f}',f'{q.FDR:.5f}'])
y=table(rows,y,[83,20,27,40,47,48],10)
y=para('主分析结果属于名义显著信号，以上FDR均未达到0.05。身份层5位的若干程序5/5升高，P=0.0625受小样本精确检验分辨率限制。胶原加工面板不含COL家族，降低仅由共享胶原基因重复计证据的影响。',y,'small')
y=para('Lung_04仍伴成骨程序增强；既有结果不支持将其直接解释为实际矿化增强。优先推荐的是基质相关表达群，而不是已证实纤维沉积机制的细胞类型。',y,'small')
y=start('全部有信号群保留：推荐与探索分层')
notes={'Bone_00':'基质/胶原表达；单患者贡献高','Bone_07':'骨侧基质与收缩线索','Bone_07_s1':'胶原方向线索；加工面板不一致','Bone_07_s2':'次级基质线索；来源间方向不一致','Bone_13':'胶原/成骨重叠；单患者','Lung_04':'优先推荐：患者支持与功能一致性较好','Lung_04_s2':'父群内基质偏向；配对患者少','Lung_05_s2':'胶原线索；免疫/内皮标记需复核','Lung_05_s4':'胶原线索；免疫/内皮标记需复核','Lung_06_s2':'方向性线索；少患者、结构稳定性有限','Lung_07':'收缩+基质；单患者及血管壁线索','Lung_11':'基质表达；缺少足够同患者对照','Lung_12':'胶原加工/基质；单患者'}
rows=[['群','细胞','min5有效N*','定位与保留原因']]
for node in sig:
 q=decision.loc[node];rows.append([node,int(q['细胞']),int(q['最少5细胞有效患者']),notes[node]])
y=table(rows,y,[43,22,35,165],8.8)
y=para('*选定来源内的主分析有效配对N；不是总覆盖患者数。所有组保留，不以本表的N否定其表达特征。主文效应图展示10个重点群，另外3组在本表和13组marker图中完整保留。',y,'small')
y=start('方法与验收：图件不改变原有结果')
for text in ['<b>统计单位：</b>同来源、同病灶、同患者目标与对照均值差；每位患者等权。父群对全病灶其余候选，子群以同父群其余细胞为主要对照。样本/建库重复按原donor关系合并。','<b>功能与文献评分：</b>采用未整合的log1p(CP10k)表达和既有评分；FNDC1等论文面板按原文基因列表复核。气泡图为描述性表达展示，不能替代患者级检验。','<b>检验：</b>N≤16使用双侧精确符号翻转；更大N使用65,536次蒙特卡洛。N&lt;2或全零不补造P。P&lt;0.05用于名义线索；FDR在预设同类家族内计算并完整保留。','<b>身份敏感性：</b>沿用成纤维支持层；NCI限定作者diploid，并对明确免疫/内皮组合共表达作敏感性排除。该层不等于实验验证纯度，低CNV或diploid也不是CAF身份的充分证明。','<b>验收：</b>既有15,792项功能对比重新计算一致；UMAP展示6,423个主分析细胞；原细胞、分群及输入哈希未变。本次新增图件与PDF独立保存，不覆盖Phase119/120。','<b>解释范围：</b>未建立独立留出患者重复；文献相似性不是外部数据验证。RNA程序不直接证明蛋白分泌、纤维沉积、矿化、耐药或转移机制。']:
 y=para(text,y)
for chunk,title in [(lit[:5],'文献交叉证据（一）'),(lit[5:],'文献交叉证据（二）')]:
 y=start(title)
 for a in chunk:
  idx=int(a['id'][1:]);y=para(f'<b>[{idx}] {a["title"]}</b><br/>{a["finding"]} {a["boundary"]}<br/><font color="#1c6470">DOI: {a["doi"]}</font>',y,'small')
  y-=5
y=start('推荐结论：优先推进Lung_04，其余信号保留')
y=para('<b>推荐Lung_04作为本轮优先研究候选。</b> 在当前固定队列中，它同时具备较好的患者覆盖、基质/胶原程序一致性及FNDC1文献方向支持。',y)
y=table([['推荐依据','结果'],['主分析患者支持','7位有效患者，达到本项目建议下限6位；来源为NCI'],['功能方向','基质、胶原交联、纤维型胶原及FNDC1面板均7/7同向'],['非COL加工支持','6/7同向；避免仅依赖共享胶原基因'],['身份敏感性','5位可配对患者，基质与FNDC1仍5/5同向；配对损失主要来自对照不足'],['其他群','Bone_07及其他有信号父/子群全部保留，按患者支持和身份证据分层追踪']],y,[58,207],11)
y=para('<b>建议表述：</b>“在骨、肺病灶分别细分后，识别出多个基质相关表达群。其中Lung_04获得较充分的患者级一致性支持，优先推荐为后续验证的ECM/胶原相关CAF候选群。”',y)
y=para('后续验证重点：补充独立患者或来源，确认成纤维身份及基质蛋白表达，再检验胶原分泌/沉积。当前推荐不依赖宣称矿化、耐药或促转移机制已被证实。',y,'small')
c.save()
doc=fitz.open(out);assert len(doc)==19;full='\n'.join(p.get_text() for p in doc);assert 'Lung_04' in full and '0.06875' in full and '通过为0' not in full and '最终准入' not in full
for node in sig:assert node in full
assert '\ufffd' not in full
receipt=dict(status='PDF_CREATED_TEXT_QA_PASS',pages=len(doc),retained_signal_nodes=sig,priority='Lung_04',policy='Main patient support guides recommendation; identity analysis is sensitivity, not veto; historical statistics unchanged',pdf=str(out),sha256=hashlib.sha256(out.read_bytes()).hexdigest(),bytes=out.stat().st_size,visual_QA='pending',created_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(W/'PDF_BUILD_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8');(W/'REPORT_POLICY_UPDATE.json').write_text(json.dumps(dict(user_instruction='不使用通过为0或全部最终准入条件；保留有信号群，依据主分析患者支持优先推荐Lung_04',scope='Report interpretation and recommendation only; no alteration of cohort, labels, values, P or FDR',date='2026-10-26'),ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(receipt,ensure_ascii=False))
