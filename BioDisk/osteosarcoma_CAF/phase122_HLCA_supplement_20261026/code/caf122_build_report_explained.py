from pathlib import Path
import json,zipfile,hashlib
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,Table,TableStyle,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.enums import TA_CENTER
import pymupdf
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');D=B/('review_bundle_final' if (B/'review_bundle_final').exists() else 'review_bundle');F=B/'figures';F.mkdir(exist_ok=True)
for name,file in [('Song','simsun.ttc'),('Hei','simhei.ttf')]:pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
styles={'body':ParagraphStyle('body',fontName='Song',fontSize=10.5,leading=17,spaceAfter=8),'title':ParagraphStyle('title',fontName='Hei',fontSize=21,leading=30,spaceAfter=18),'h':ParagraphStyle('h',fontName='Hei',fontSize=15,leading=22,spaceAfter=12),'small':ParagraphStyle('small',fontName='Song',fontSize=9,leading=14,spaceAfter=6)}
for style in styles.values():style.wordWrap='CJK'
rec=json.loads((D/'ANALYSIS_COMPLETE.json').read_text(encoding='utf-8'));qa=json.loads((D/'SUPPLEMENT_QA.json').read_text(encoding='utf-8'));notes=json.loads((B/'REPORT_INTERPRETATION.json').read_text(encoding='utf-8'));assert qa['status']=='PASS'
idn=pd.read_csv(D/'tables/identity_patient_threshold_summary.tsv',sep='\t');bg=pd.read_csv(D/'tables/normal_IPF_background_comparison.tsv',sep='\t');rd=pd.read_csv(D/'tables/reference_donor_program_scores.tsv',sep='\t');qd=pd.read_csv(D/'tables/query_candidate_patient_program_scores.tsv',sep='\t');tt=pd.read_csv(D/'tables/IPF_within_study_patient_tests.tsv',sep='\t')
nodes=['Lung_04','Bone_07','Bone_07_s1','Lung_04_s2'];names=['mCAF_matrix','collagen_biosynthesis','collagen_crosslinking','fibrillar_I_III_V','collagen_processing_nonCOL','HLCA_IPF_reference'];labels=['Matrix program','Collagen synthesis','Crosslinking','Fibrillar collagen','Non-COL processing','HLCA IPF markers'];palette={'Fibroblast':'#238b8d','Pericyte':'#e3ae42','Smooth_muscle':'#a16ba2','Other_lineage':'#9ca3af'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
fig,ax=plt.subplots(figsize=(9,4.2));bottom=np.zeros(4)
for typ,col in palette.items():
 vals=idn[(idn.minimum_candidate_cells==5)&(idn.identity==typ)].set_index('node_id').reindex(nodes).equal_patient_fraction.to_numpy();ax.bar(nodes,vals,bottom=bottom,color=col,label=typ.replace('_',' '));bottom+=vals
ax.set_ylim(0,1);ax.set_ylabel('Mean patient-level cell fraction');ax.legend(loc='upper center',bbox_to_anchor=(.5,1.2),ncol=4,frameon=False);fig.tight_layout();fig.savefig(F/'identity.png');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(10,4.6),sharey=True)
scale=bg[(bg.query_min_cells==5)&(bg.reference_min_cells==20)&(bg.condition=='Healthy')&bg.reference_type.isin(['Alveolar fibroblasts','Adventitial fibroblasts'])&bg.program.isin(names)].difference.abs().max()
for ax,typ in zip(axs,['Alveolar fibroblasts','Adventitial fibroblasts']):
 b=bg[(bg.query_min_cells==5)&(bg.reference_min_cells==20)&(bg.condition=='Healthy')&(bg.reference_type==typ)];mat=b.pivot(index='node_id',columns='program',values='difference').reindex(index=nodes,columns=names).to_numpy();im=ax.imshow(mat,vmin=-scale,vmax=scale,cmap='RdBu_r',aspect='auto');ax.set_title(typ);ax.set_yticks(range(4),nodes);ax.set_xticks(range(6),labels,rotation=60,ha='right')
 for i in range(4):
  for j in range(6):
   if np.isfinite(mat[i,j]):ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',color='white' if abs(mat[i,j])>scale*.6 else 'black',fontsize=8)
fig.subplots_adjust(bottom=.4,right=.89,wspace=.08);fig.colorbar(im,cax=fig.add_axes([.91,.4,.02,.45]),label='Candidate - healthy score');fig.savefig(F/'background.png',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4.1));groups=[];xlabels=[];colorset=[]
for cond,col in [('Healthy','#7cae91'),('IPF','#d57b62')]:
 a=rd[(rd.cells>=20)&(rd.ann_finest_level=='Alveolar fibroblasts')&(rd.lung_condition==cond)];groups.append(a.HLCA_IPF_reference.to_numpy());xlabels.append(cond+'\n(n='+str(len(a))+')');colorset.append(col)
for node in nodes:
 a=qd[(qd.cells>=5)&(qd.node_id==node)];groups.append(a.HLCA_IPF_reference.to_numpy());xlabels.append(node+'\n(n='+str(len(a))+')');colorset.append('#547fba')
rng=np.random.default_rng(122)
for i,(v,col) in enumerate(zip(groups,colorset)):
 ax.scatter(i+rng.uniform(-.15,.15,len(v)),v,s=14,color=col,alpha=.7);ax.plot([i-.22,i+.22],[np.median(v)]*2,color='black',lw=2)
ax.set_xticks(range(len(groups)),xlabels,rotation=15);ax.set_ylabel('HLCA IPF 6-gene score\nmean log1p(CP10k)')
for name,col in [('Healthy','#7cae91'),('IPF','#d57b62'),('CAF candidates','#547fba')]:ax.scatter([],[],s=18,color=col,label=name)
ax.plot([],[],color='black',lw=2,label='Median');ax.legend(frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(.5,1.16),fontsize=8);fig.tight_layout();fig.savefig(F/'ipf.png');plt.close(fig)
marker=pd.read_csv(D/'tables/marker_plot_summary.tsv',sep='\t');ipfgenes=['CCL2','COL1A1','CTHRC1','MMP19','SERPINE1','HIF1A'];marker_genes=ipfgenes+sorted(set(marker.gene)-set(ipfgenes));groups=['Healthy | Alveolar fibroblasts','Healthy | Adventitial fibroblasts','IPF | Alveolar fibroblasts','IPF | Adventitial fibroblasts']+nodes;rowlabels=['Healthy alveolar','Healthy adventitial','IPF alveolar','IPF adventitial']+nodes
fig,ax=plt.subplots(figsize=(10,4.6));sub=marker[marker.group.isin(groups)];vmax=sub.mean_logCP10k.max()
for i,group in enumerate(groups):
 g=sub[sub.group==group].set_index('gene').reindex(marker_genes);sc=ax.scatter(np.arange(len(marker_genes)),np.full(len(marker_genes),i),s=130*g.fraction_detected,c=g.mean_logCP10k,vmin=0,vmax=vmax,cmap='YlOrRd',edgecolor='#777777',linewidth=.25)
ax.set_yticks(range(8),rowlabels);ax.set_xticks(range(len(marker_genes)),marker_genes,rotation=60,ha='right',fontsize=8);ax.invert_yaxis();ax.set_xlim(-.7,len(marker_genes)-.3);ax.axvline(5.5,color='#bbbbbb',ls='--',lw=.7);ax.axhline(3.5,color='#bbbbbb',lw=.7);fig.colorbar(sc,ax=ax,pad=.02,label='Mean log1p(CP10k)')
for fraction in [.1,.5,1]:ax.scatter([],[],s=130*fraction,c='white',edgecolor='#555555',label=f'{fraction:.0%}')
ax.legend(title='Detection fraction',loc='upper center',bbox_to_anchor=(.5,1.24),ncol=3,frameon=False,fontsize=8,title_fontsize=9);fig.tight_layout();fig.savefig(F/'markers.png',bbox_inches='tight');plt.close(fig)
out=Path('C:/Users/31937/Downloads/骨肺CAF四项候选_HLCA补充分析报告_说明完善版_20261027.pdf') if notes.get('final_ready') else B/'supplement_internal_draft.pdf';story=[]
def p(text,style='body'):story.append(Paragraph(text,styles[style]))
def heading(s):p(s,'h')
def figure(name,width=490):
 from PIL import Image as PILImage
 w,h=PILImage.open(F/name).size;story.append(Image(str(F/name),width=width,height=width*h/w))
def table(rows,widths):
 wrapped=[[Paragraph(str(x),styles['small']) for x in row] for row in rows];t=Table(wrapped,colWidths=widths,repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eff3')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#587080')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#d8dee4')),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7)]));story.append(t)

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
p('主要发现与候选推荐','title');p('身份参照 · 正常肺背景 · 肺纤维化状态比较','h');p('分析基准：原 6,633 个骨肺细胞；主分析保留 6,423 个。四项候选及原有层级保持不变。完成日期：2026-10-27。')
for text in notes['main']:p(text)
rows=[['候选','细胞数','候选患者数 ≥5 / ≥10细胞','定位']]
for n,c in zip(nodes,[233,131,35,36]):
 ns=rec['query_candidate_patient_counts'][n];rows.append([n,c,f"{ns['5']} / {ns['10']}",notes['roles'][n]])
table(rows,[90,55,150,195]);p('此处患者数用于候选群自身的外部参照比较；原报告的组内配对比较还要求同患者对照细胞充足，两者分母不同。父群与其子群存在包含关系。','small')
story.append(PageBreak());heading('一、身份参照：以患者为单位评估混合程度');figure('identity.png')
for text in notes['identity']:p(text)
p('图1｜HLCA健康核心图谱用于提供已知肺细胞身份参照。每位患者先计算候选细胞的身份比例，再对患者等权平均。纳入该候选至少5个细胞的患者；颜色图例给出全部身份类别。标签来自HLCA健康核心图谱50近邻投票。','small')
story.append(PageBreak());heading('二、正常肺背景：比较不同基质功能程序');figure('background.png')
for text in notes['background']:p(text)
p('图2｜HLCA健康肺成纤维细胞用于判断肺组织固有基质表达背景。数值为候选患者功能分数中位数减去健康参考的研究等权平均分数。单位为平均log1p(CP10k)；红色表示候选更高，蓝色表示参考更高。候选每位患者至少5细胞，参考每位供者每类至少20细胞。分别使用肺泡和外膜成纤维细胞作背景。','small')
story.append(PageBreak());heading('三、IPF功能参照：候选与纤维化程序的重叠');p('这一步承接前面的基质升高结果：进一步判断升高是否伴随更广泛的纤维化相关表达。IPF样本提供外部疾病状态参照，候选的研究价值仍结合原骨肺队列证据判断。');figure('ipf.png')
for text in notes['ipf']:p(text)
p('图3｜以HLCA中IPF（特发性肺纤维化）相关成纤维细胞作功能参照，使用原论文Fig.6h所列6基因：CCL2、COL1A1、CTHRC1、MMP19、SERPINE1、HIF1A。每个点为一位患者/供者的平均分数，黑线为组内中位数，括号为独立患者/供者数。健康与IPF点取肺泡成纤维细胞，参考每类至少20细胞；候选至少5细胞。此图为跨数据集描述比较。','small')
story.append(PageBreak());heading('四、单基因核查：识别COL1A1对IPF分数的贡献');figure('markers.png')
for text in notes['gene_evidence']:p(text)
nc=pd.read_csv(D/'tables/IPF_without_COL1A1_background.tsv',sep='\t');rows=[['候选','去COL1A1后健康百分位','去COL1A1后IPF百分位']]
for node in nodes:
 a=nc[nc.node_id==node].set_index('reference_group');rows.append([node,f"{a.loc['Healthy | Alveolar fibroblasts','median_reference_percentile']:.3f}",f"{a.loc['IPF | Alveolar fibroblasts','median_reference_percentile']:.3f}"])
table(rows,[130,180,180]);p('图4｜前6列为HLCA-IPF标志基因，其余为沿用程序中的基质与myCAF标志基因。颜色为平均log1p(CP10k)，点面积为表达细胞比例；先按供者汇总，参考再按研究等权平均。候选每患者至少5细胞，参考每供者每类至少20细胞。表中百分位使用肺泡成纤维细胞参照。','small')
story.append(PageBreak());heading('五、方法、敏感性检查与核验')
for text in notes['methods']:p(text)
heading('集中说明');p(notes['boundaries'])
heading('数据准确性核验');p(f"自动核验{qa['checks']}项，通过{qa['passed']}项。核验包含细胞标识及顺序、候选数量、原始整数计数、全库归一化分母、抽样重算功能分数、患者汇总唯一性及统计适用条件。完整患者结果、P值、FDR和文件哈希保存在服务器phase122目录。")
heading('参考资料');p('Sikkema L, et al. An integrated cell atlas of the lung in health and disease. Nature Medicine (2023). DOI: 10.1038/s41591-023-02327-2。https://pmc.ncbi.nlm.nih.gov/articles/PMC10287567/','small');p('HLCA模型与潜在空间：Zenodo记录7599104，v1.1；表达参照：CELLxGENE完整HLCA数据集9f222629-9e39-47d0-b83f-e08d610c7479，官方H5AD文件，下载后分块重算哈希。','small')
def footer(c,d):
 c.setFont('Song',8);c.setFillColor(colors.HexColor('#657382'));c.drawString(48,28,'骨肺CAF · HLCA补充分析 · 2026-10-27');c.drawRightString(547,28,str(d.page))
doc=SimpleDocTemplate(str(out),pagesize=(595,842),rightMargin=48,leftMargin=48,topMargin=45,bottomMargin=45);doc.build(story,onFirstPage=footer,onLaterPages=footer)
pdf=pymupdf.open(out);texts=[]
for i,page in enumerate(pdf):page.get_pixmap(matrix=pymupdf.Matrix(1.2,1.2)).save(B/f'explained_page_{i+1}.png');texts.append(page.get_text())
assert all(n in '\n'.join(texts) for n in nodes)
(B/'REPORT_EXPLAINED_RENDER_QA.json').write_text(json.dumps(dict(file=str(out),pages=len(pdf),sha256=hashlib.sha256(out.read_bytes()).hexdigest(),candidate_names_present=True),indent=2),encoding='utf-8');print(out,len(pdf))
