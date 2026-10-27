from pathlib import Path
import json,hashlib,shutil,html
import pandas as pd
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
import fitz
B=Path('C:/Users/31937/CAF_F2C3_subtype_review_20261010');I=B/'inputs';O=B/'output/pdf';TMP=B/'tmp/pdfs'
nompath=Path('C:/Users/31937/CAF_all_clusters_20261006/analysis/tables/nominal_upregulated_all_sources_with_sensitivity.tsv');shutil.copy2(nompath,I/nompath.name)
g=pd.read_csv(I/'SCPCP000017__F2_C3__corrected__Fibro.tsv.gz',sep='\t').set_index('symbol');allg=pd.read_csv(I/'SCPCP000017__F2_C3__corrected__all.tsv.gz',sep='\t');m=pd.read_csv(I/'marker_display.tsv',sep='\t');m=m[m.context=='SCPCP000017__F2_C3__corrected__all'];n=pd.read_csv(nompath,sep='\t');n=n[n.subgroup_id=='F2_C3'];assert len(n)==22 and len(g)==52521 and len(m)==12
markers={'COL5A1_MSC4':['COL5A1','COL16A1','EMILIN1','ADAMTS2'],'matrix_CAF6':['FAP','MMP11','COL1A1','COL11A1','POSTN','FN1'],'FNDC1_CAF6':['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3'],'atlas_mCAF4':['SPARC','THBS2','FN1','POSTN']};rows=[]
for panel,gs in markers.items():
 for gene in gs:
  r=g.loc[gene];rows.append({'literature_panel':panel,'gene':gene,'C3_detection':r.target_detection,'other_F2_Fibro_detection':r.reference_detection,'descriptive_log2_difference':r.mean_delta,'in_nominal22':gene in set(n.gene),'in_report_phenotype3':gene in ['CD248','POSTN','VCAN']})
pd.DataFrame(rows).to_csv(I/'literature_marker_expression_crosswalk.tsv',sep='\t',index=False);m.to_csv(I/'original_C3_top12_display.tsv',sep='\t',index=False)
refs=[
('Zheng et al. J Hematol Oncol. 2024.','A single-cell and spatially resolved atlas of human osteosarcomas.','10.1186/s13045-024-01598-7'),
('Zheng et al. J Cancer Res Clin Oncol. 2024;150:356.','Single-cell transcriptomic insights into chemotherapy-induced remodeling of the osteosarcoma tumor microenvironment.','10.1007/s00432-024-05787-2'),
('Zeng et al. Molecular Cancer. 2026;25:14.','Targeting OxLDL-mediated CD36+ CAF reprogramming potentiates PD-1 immunotherapy in osteosarcoma.','10.1186/s12943-025-02516-2'),
('He et al. Journal of Bone Oncology. 2024;46:100604.','Single-cell transcriptional profiling in osteosarcoma and the effect of neoadjuvant chemotherapy on the tumor microenvironment.','10.1016/j.jbo.2024.100604'),
('Liu et al. BMC Medicine. 2024.','Single-cell and spatial transcriptomics reveal metastasis mechanism and microenvironment remodeling of lymph node in osteosarcoma.','10.1186/s12916-024-03319-w'),
('Budhathoki et al. bioRxiv. 2026. 预印本v1。','Integrative Single-cell and Spatial Transcriptomic Analysis of Osteosarcoma Reveals Conserved and Distinct Ecosystems Across Sites and Species.','10.64898/2026.01.13.698472'),
('Huang et al. Theranostics. 2022;12:5877-5887.','Single-cell transcriptomics reveals the regulative roles of cancer associated fibroblasts in tumor immune microenvironment of recurrent osteosarcoma.','10.7150/thno.73714'),
('Hu et al. Cancer Drug Resistance. 2026;9:13.','Paired single cell analysis reveals chemotherapy resistance in osteosarcoma.','10.20517/cdr.2026.09'),
('Hu et al. J Immunother Cancer. 2024;12:e006991.','Cell membrane-anchored and tumor-targeted IL-12 T-cell therapy destroys cancer-associated fibroblasts and disrupts extracellular matrix in heterogenous osteosarcoma xenograft models.','10.1136/jitc-2023-006991'),
('Rossow et al. Oncogene. 2018;37:4921-4940.','LOX-catalyzed collagen stabilization is a proximal cause for intrinsic resistance to chemotherapy.','10.1038/s41388-018-0320-2')]
(I/'verified_reference_registry.json').write_text(json.dumps([{'number':i+1,'citation':a,'title':t,'doi':doi,'url':'https://doi.org/'+doi,'origin':'客户清单' if i+1 in [1,4,7,8,9,10] else '补充原始研究'} for i,(a,t,doi) in enumerate(refs)],ensure_ascii=False,indent=2),encoding='utf-8')
pdfmetrics.registerFont(TTFont('YaHei','C:/Windows/Fonts/msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('YaHeiB','C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0));pdfmetrics.registerFontFamily('YaHei',normal='YaHei',bold='YaHeiB',italic='YaHei',boldItalic='YaHeiB')
C=colors.HexColor;sty={
'body':ParagraphStyle('body',fontName='YaHei',fontSize=10.2,leading=16.9,spaceAfter=9,wordWrap='CJK',textColor=C('#263847')),
'small':ParagraphStyle('small',fontName='YaHei',fontSize=8.7,leading=13.2,spaceAfter=6,wordWrap='CJK',textColor=C('#536675')),
'h1':ParagraphStyle('h1',fontName='YaHeiB',fontSize=21,leading=29,spaceAfter=16,textColor=C('#143e54')),
'h2':ParagraphStyle('h2',fontName='YaHeiB',fontSize=15,leading=22,spaceAfter=12,textColor=C('#176c86')),
'h3':ParagraphStyle('h3',fontName='YaHeiB',fontSize=11.7,leading=18,spaceBefore=7,spaceAfter=7,textColor=C('#20596d')),
'cell':ParagraphStyle('cell',fontName='YaHei',fontSize=9.1,leading=14.3,wordWrap='CJK',textColor=C('#263847')),
'ref':ParagraphStyle('ref',fontName='YaHei',fontSize=8.6,leading=12.5,spaceAfter=10,wordWrap='CJK',textColor=C('#334957'))}
W,H=A4;fw=W-92;story=[]
def p(s,style='body'):return Paragraph(s,sty[style])
def add(s,style='body'):story.append(p(s,style))
def heading(s):add(s,'h2')
def table(head,rs,widths):
 data=[[p('<b>'+html.escape(str(v))+'</b>','cell') for v in head]]+[[p(str(v),'cell') for v in row] for row in rs];t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),C('#eaf2f5')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.6,C('#b8cdd6')),('LINEBELOW',(0,1),(-1,-1),.35,C('#dce5e9'))]));story.extend([t,Spacer(1,10)])
def box(s):
 t=Table([[p(s)]],colWidths=[fw]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),C('#edf5f7')),('BOX',(0,0),(-1,-1),.6,C('#c8dce4')),('LEFTPADDING',(0,0),(-1,-1),13),('RIGHTPADDING',(0,0),(-1,-1),13),('TOPPADDING',(0,0),(-1,-1),10),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,12)])
def page():story.append(PageBreak())
def cite(i):return f'[{i}]'
def markerline(gs):return '、'.join(gs)
def dval(gene):return f'{g.loc[gene,"mean_delta"]:+.2f}'
def det(gene):return f'{g.loc[gene,"target_detection"]*100:.1f}%'
# Page 1
add('F2-C3骨肉瘤亚群查新核查','h1');add('核查日期：2026年10月10日  |  对照对象：2026年10月6日报告','small')
box('<b>当前判断：</b>F2-C3与已报道的骨肉瘤胶原／基质相关细胞状态存在实质重叠，现有证据不足以认定为“首次发现的新亚群”。同时，尚未证明它与某个已报道亚群完全等同。更稳妥的定位是<b>已知基质相关框架内的候选细分状态</b>。')
heading('一、F2-C3在本项目中如何定义')
add('本报告只回答三个问题：本群如何定义、既往有哪些相近亚群、现有证据能否支持新亚群。以10月6日报告及其原始表达记录为本群依据，核查客户清单6篇文献，并补充4篇直接相关的原始研究／图谱。')
table(['定义层次','核实结果'],[
['聚类身份','F2-C3是本项目F2内部的工作聚类编号。冻结记录共38个细胞；主要来源SCPCP000017为34个细胞，患者级对照有效N=7。编号本身不代表跨研究通用亚群。'],
['表型描述','原报告以CD248、POSTN、VCAN及胶原程序描述C3。三者在主来源C3中的校正表达检出比例分别为85.3%、94.1%、94.1%；不是已验证的排他性鉴别组合。'],
['完整表达特征','已核验52,521条全转录组对比记录。原展示前12项为CD248、MMP9、OLFML3、VCAN、CFH、POSTN、P4HA2、COL16A1、PRSS23、COL5A3、FKBP11、AP001107.3；展示排名属于描述性结果。'],
['胶原／基质特征','III／V型胶原、PCOLCE、LOXL2及其他基质成员构成功能描述。22个名义上调基因来自预设机制基因范围，不能代替完整亚群身份签名。']], [85,fw-85])
add('<b>22个名义上调基因：</b>'+markerline(n.gene.tolist())+'。','small')
add('上述22项为原始P&lt;0.05、原家族FDR未通过0.05的探索性结果。P值描述表达比较证据，不是“新亚群概率”。下文依据基因组合及整体方向进行亚群对照。','small')
# Page2
page();heading('二、与骨肉瘤已报道亚群交叉对照')
add('最直接的两个对照：COL5A1+间质细胞与matrix CAF','h3')
add('<b>COL5A1+ MSC：最直接的胶原型先例。[2]</b> Zheng等2024年已在骨肉瘤中定义COL5A1+间质细胞，展示COL5A1、COL16A1、EMILIN1、ADAMTS2。四者均在C3中检出，且在原报告同F2、Fibro限定的表达对照中全部较高。这证明两者具有明确的共同表达特征，不能仅以“胶原高表达”作为新亚群依据。')
add('<b>Matrix CAF／mCAF：最直接的既有基质型框架。[1]</b> 补充图S2I展示FAP、MMP11、COL1A1、COL11A1、POSTN、FN1。C3具有其中多项特征；POSTN、COL1A1、COL11A1和FN1相对其他F2 Fibro较高，而FAP、MMP11在这一背景下没有同步升高。该差异可作为细分状态线索，尚不能据此建立独立新亚群。')
rs=[]
for gene in ['COL5A1','COL16A1','EMILIN1','ADAMTS2','POSTN','COL1A1','COL11A1','FN1','FAP','MMP11']:
 rs.append([gene,'COL5A1+ MSC' if gene in markers['COL5A1_MSC4'] else 'matrix CAF',det(gene),dval(gene)])
table(['文献展示基因','对应已报道群','C3检出比例','表达log2差'],rs,[88,137,115,fw-340])
add('表中取原报告同一来源的C3（34细胞）对其他F2内Fibro（70细胞）的描述性伪总体log2(CP10k+0.1)差，校正计数检出门槛为0.5。它不是文献原细胞的直接差异检验，也不是患者级P值。文献短表为展示标志，不等同作者完整分类器。','small')
#Page3
page();heading('二、交叉对照（续）：其他直接相关亚群')
add('FNDC1+ CAF：基质功能重叠，代表组合不完全吻合 [3]','h3')
add('Zeng等研究定义FNDC1、COMP、MMP2、COL11A1、COL10A1、P4HA3组成的FNDC1+ CAF。原C3相对其他F2 Fibro中，FNDC1、COMP、MMP2、P4HA3的描述性差为'+', '.join(dval(z) for z in ['FNDC1','COMP','MMP2','P4HA3'])+'；COL11A1为'+dval('COL11A1')+'，COL10A1仅3/34细胞检出。现有表达不支持将C3直接改称FNDC1+ CAF，但也不足以证明它属于全新类别。')
add('He等的SPP1、IGFBP4和COL10A成纤维细胞 [4]','h3')
add('该研究确实定义了三个成纤维细胞群，应全部纳入对照。C3中的SPP1和COL10A1检出分别为4/34、3/34；IGFBP4则为25/34，描述性log2差为+1.44。因此，不能只对照SPP1／COL10A群而忽略IGFBP4群。尚未取得足以完成整体同群判断的作者完整marker与细胞标签；客户表中的扩展marker不能直接当作已全部核实的标准签名。')
add('Liu等的myCAF、OstCAF及apCAF [5]','h3')
add('该文在骨肉瘤相关淋巴结中区分myCAF、OstCAF和apCAF，分别展示POSTN／TAGLN、ALPL／RUNX2以及抗原呈递相关分子。C3与其中的基质标志存在共享，但组织环境和亚群定义不同。名称或某个共同marker均不足以判定两群等同；该文也不能因未使用“F2-C3”名称就被排除。')
add('OsteoCAR图谱预印本：同来源的优先对照 [6]','h3')
add('已核查的作者公开代码列出mCAF相关展示签名SPARC、THBS2、FN1、POSTN。该图谱包含本项目主要来源SCPCP000017，必须进行来源去重；它是相关的既有图谱，不能直接计为独立复现。本次未取得C3与作者细胞标签的一一对应证据，因此不判定“已经重复”，也不据此宣称首次发现。')
add('本页对照结论','h3')
box('C3并非只与一个文献短名单存在联系：它同时落在胶原型／基质型表达框架中，并与IGFBP4等已报道群的部分特征相交。当前最需要区分的是<b>已有群体中的细分状态</b>与<b>可独立定义的新亚群</b>，不是比较聚类编号或名称是否相同。')
#Page4
page();heading('二、客户其余参考文献能否直接判断重复')
table(['客户文献','实际比较层次','对“是否新亚群”的作用'],[
['Huang等，2022 [7]','CAF与Div-CAF等较粗分层，并研究复发相关CAF。','证明骨肉瘤CAF及其异质性已有报道；未提供可与C3完整定义直接对应的精细亚群证据。'],
['Hu等，2026 [8]','总体CAF变化；所列九基因签名来自恶性细胞。','没有在相应比较中细分CAF。不能把恶性细胞签名当作C3亚群身份对照，也不能由“不重叠”证明C3新颖。'],
['Hu等，2024 [9]','总体CAF及ECM签名，主要为PDX研究。','未定义与C3对应的精细CAF群。它提供共同基质背景，不能直接裁定C3是否已经被报道。'],
['Rossow等，2018 [10]','非骨肉瘤的胶原／基质研究，关键过表达对象为4T1肿瘤细胞。','不属于骨肉瘤CAF亚群发现研究；保留为客户清单核查项，不作为本次亚群重复判据。']], [103,169,fw-272])
add('客户查新表中会影响亚群判断的三处口径','h3')
add('① CD248／POSTN／VCAN应称“原报告描述的表型特征”，不能直接写成已验证的排他性核心marker。<br/>② 22个名义上调基因是机制范围内的候选，不是完整的亚群marker集。<br/>③ Hu 2026的九基因属于恶性细胞；Rossow 2018并非骨肉瘤CAF亚群研究。两者不能用于证明“C3是新群”或“C3已被做过”。')
add('短名单交集只能辅助描述','h3')
ng=set(n.gene);rs=[]
for name,gs in markers.items():
 label={'COL5A1_MSC4':'COL5A1+ MSC四基因 [2]','matrix_CAF6':'matrix CAF六基因 [1]','FNDC1_CAF6':'FNDC1+ CAF六基因 [3]','atlas_mCAF4':'图谱mCAF四基因 [6]'}[name];shared=sorted(ng&set(gs));rs.append([label,str(len(shared)), '、'.join(shared) if shared else '无'])
table(['文献展示短表','与22项交集数','共有基因'],rs,[239,104,fw-343])
add('例如COL1A1和COL11A1不在这22项中，但均在C3表达并呈描述性正向差。故“交集为零／很小”不能解释成“不属于已知群”。因文献完整marker集及统一筛选背景不齐备，本报告不作超几何显著性检验，也不把交集数换算成新颖性分数。','small')
#Page5
page();heading('三、新颖性判定与可直接使用的结论')
table(['需要回答的问题','核查结论'],[
['本项目是否有可描述的C3状态？','有。原编号、细胞记录及完整表达特征可追溯，具有CD248／POSTN／VCAN及胶原相关表达组合。'],
['骨肉瘤是否已有相近亚群？','有。COL5A1+间质细胞和matrix CAF是最直接先例；FNDC1+ CAF、IGFBP4等群及同来源图谱也需纳入。'],
['是否足以认定“首次发现的新亚群”？','不足。共有表达框架明确，现有证据尚未建立区别于已知亚群的完整、稳定鉴别定义。'],
['是否可以判定“已被某篇文章完整做过”？','目前也不能。尚缺少本C3与作者参考细胞的整体表达及标签对应证据，不能仅凭功能相近或少量共有基因认定完全重复。']], [176,fw-176])
add('综合判断','h3')
box('<b>现有证据更支持将F2-C3定位为：已知胶原／基质相关细胞框架内的候选细分状态。</b><br/>“新亚群”目前属于未建立的主张；“与某个已知亚群完全等同”同样尚未建立。保留F2-C3工作编号，不把编号或新名称本身作为创新证据。')
add('建议对客户的回复','h3')
add('“已对10月6日报告的F2-C3与骨肉瘤相关文献做了交叉核查。既往已报道胶原型和基质型细胞亚群，其中COL5A1+间质细胞及matrix CAF与我们的C3最为接近；C3也有部分不同的表达组合。因此，目前不宜直接称为首次发现的新亚群，更适合定位为已有基质相关框架中的候选细分状态。同时，现有证据还不能说它与某篇文献中的亚群完全重复。”')
add('本轮核查完成范围','h3')
add('已完成原报告定义及全谱记录核验、客户6篇与补充4篇文献的逐篇对照，以及分层新颖性判断。文献依据为原始论文、相应补充图和作者公开代码；保留作者各自命名，不以myCAF／mCAF名称决定是否等同。')
add('判断边界：多篇文献仅提供展示marker，完整参考标签未统一取得；预印本作为相关先例单列。本次为截至2026年10月10日的定向查新，不把“未找到同名群”当作“首次发现”的证明，也不将本报告表述为已完成跨论文细胞的一一同一性验证。','small')
# Page6
page();heading('参考文献与核查依据')
for i,(author,title,doi) in enumerate(refs,1):
 add(f'<b>[{i}] {html.escape(author)}</b><br/>{html.escape(title)}<br/><link href="https://doi.org/{doi}" color="#176c86">DOI: {doi}</link>','ref')
add('资料对应','h3')
add('客户清单对应[1]、[4]、[7]-[10]；新增直接对照为[2]、[3]、[5]、[6]。[1]基因短表核至补充图S2I；[2]核至正文及S3F；[3]核至Results的CAF亚群定义；[4]核实原文三群命名，未将客户扩展列表当作完整签名；[5]核至CAF分群正文；[6]为v1预印本及作者公开代码。','small')
add('本群依据：2026年10月6日胶原与基质整合报告及同日功能分类合并稿；Phase62完整表达对比与展示marker；Phase69患者名义上调表。公开代码：<link href="https://github.com/kidcancerlab/Osteo_atlas" color="#176c86">kidcancerlab/Osteo_atlas</link>，已核查快照2829c082ff7aa0c60d1146b747f356aa55fb0f73。','small')

def footer(c,doc):
 c.saveState();c.setStrokeColor(C('#ccdce3'));c.setLineWidth(.6);c.line(46,H-40,W-46,H-40);c.setFont('YaHei',8);c.setFillColor(C('#607784'));c.drawString(46,H-29,'CAF研究 | F2-C3亚群查新');c.drawRightString(W-46,H-29,'2026-10-10');c.line(46,40,W-46,40);c.drawString(46,26,'定义核验 · 文献对照 · 新颖性判断');c.drawRightString(W-46,26,f'{doc.page} / 6');c.restoreState()
output=O/'F2-C3骨肉瘤亚群查新_三项核查报告_20261010.pdf';doc=SimpleDocTemplate(str(output),pagesize=A4,rightMargin=46,leftMargin=46,topMargin=58,bottomMargin=54,title='F2-C3骨肉瘤亚群查新：定义、文献对照与新颖性判断',author='CAF项目分析',allowSplitting=1);doc.build(story,onFirstPage=footer,onLaterPages=footer)
d=fitz.open(output);texts=[p.get_text() for p in d];(TMP/'pdf_extracted.txt').write_text('\n\n'.join(texts),encoding='utf-8');assert len(d)==6,('Unexpected page count',len(d));assert all(len(t)>200 for t in texts);assert sum(len(p.get_links()) for p in d)>=11;assert not any('\ufffd' in t for t in texts);assert 'COL5A1' in texts[1] and '新颖性判定' in texts[4];print('PDF',output,'PAGES',len(d),'BYTES',output.stat().st_size)
qa={'status':'text_and_structure_passed_pending_visual_review','pages':len(d),'page_text_lengths':[len(t) for t in texts],'links':sum(len(p.get_links()) for p in d),'full_expression_records_checked':len(g),'nominal_genes_checked':len(n),'original_display_markers_checked':len(m),'literature_records':len(refs),'scope':['C3 definition','published subtype comparison','novelty judgment'],'sha256':hashlib.sha256(output.read_bytes()).hexdigest()};(B/'PDF_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
