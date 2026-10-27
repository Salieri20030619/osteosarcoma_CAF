"""Build the final review package only after real cell-identity integration.
SYNTAX_CHECK: python -m py_compile build_delivery.py
RUN_COMMAND: python build_delivery.py --root <phase2> --phase3 <phase3>
EXPECTED_SUCCESS_MARKER: PHASE3_DELIVERY_COMPLETE
"""
import argparse,json,hashlib,zipfile,datetime
from pathlib import Path
import pandas as pd
SHORTS=['S0058','S0059','S0217','S0218','SC069','SC072','SC073']

def table(df):
    head='| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'
    return head+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--phase3',type=Path,required=True);a=ap.parse_args();p=a.phase3;idd=p/'identity'
assert (idd/'COMPLETE').exists()
delivery=p/'delivery';delivery.mkdir(exist_ok=True)
s=pd.read_csv(idd/'sample_identity_summary.tsv',sep='\t')
qc=pd.read_csv(idd/'tool_reference_QC.tsv',sep='\t')
cells=pd.read_csv(idd/'all_cell_identity_evidence.tsv.gz',sep='\t')
inputs=pd.read_csv(p/'review/input_summary.tsv',sep='\t')
coords=pd.read_csv(p/'resources/coordinate_audit.tsv',sep='\t')
changes=pd.read_csv(p/'review/annotation_changes.tsv',sep='\t')
manifest=pd.read_csv(a.root/'code/manifest.tsv',sep='\t')
manifest['dataset']='GSE270231';manifest['tissue']='Lung metastasis';manifest['clinical_group']='B2';manifest['source_status']='patient_direct_tissue_per_GEO_and_paper';manifest['geo_reported_assembly']='hg19';manifest['matched_feature_annotation']='GRCh38_Ensembl84';manifest['patient_mapping_status']=manifest.patient_id.map(lambda x:'pending' if x=='UNCONFIRMED' else 'known_within_dataset')
manifest.to_csv(delivery/'source_sample_manifest.tsv',sep='\t',index=False)
clinical=[]
geo_titles={};current_sample=None
for line in (p/'takeover_audit/GSE270231_family.soft').read_text(encoding='utf-8').splitlines():
    if line.startswith('^SAMPLE = '):current_sample=line.split(' = ',1)[1]
    elif line.startswith('!Sample_title = ') and current_sample:geo_titles[current_sample]=line.split(' = ',1)[1]
for r in manifest.itertuples():
    for field,raw,verified,status in [('patient_id',r.patient_id,r.patient_id,r.patient_mapping_status),('clinical_group',geo_titles.get(r.sample,'SOURCE_TITLE_NOT_CAPTURED'),'B2','source_supported'),('technology',r.technology,r.technology,'phase2_manifest_value_preserved_RNA_Peaks_discrepancy_retained'),('treatment_history','FIELD_NOT_PROVIDED','unknown','pending'),('bone_lung_sampling_interval','FIELD_NOT_PROVIDED','not_available','pending'),('genome_assembly','hg19','GRCh38 gene-coordinate annotation','GEO_conflict_retained_feature_match_33694_of_33694')]:
        clinical.append({'dataset':'GSE270231','sample_id':r.sample,'field':field,'original_value':raw,'reviewed_value':verified,'status':status,'source':'phase2 code/manifest.tsv; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+r.sample+'; coordinate_audit.tsv for gene annotation'})
pd.DataFrame(clinical).to_csv(delivery/'clinical_source_audit.tsv',sep='\t',index=False)
counts=cells.final_identity.value_counts();strict=int(cells.caf_strict.sum());expanded=int(cells.caf_extended.sum());discord=int((cells.cnv_concordance=='discordant_or_ambiguous').sum())
raw_conflicts=int(cells.original_copykat_conflict.sum())
selection=pd.read_csv(idd/'method_selection.tsv',sep='\t')
extra_completed=sum((p/m/s/'COMPLETE').exists() for m in ['copykat_balanced','scevan'] for s in SHORTS)
ss=s[['short','patient_id','retained','cnv_input','caf_strict','caf_extended','Identity_pending']].copy().fillna(0)
ss.columns=['样本','患者ID','保留细胞','CNV输入','严格CAF候选','扩展CAF候选','身份待定']
for c in ss.columns[2:]:ss[c]=ss[c].astype(int)
qbad=qc[~qc.reference_check_pass]
figcap='图中 UMAP1/UMAP2 为复用的分样本 UMAP 坐标，不是不同样本间可直接比较的坐标。每点代表一个细胞/核；颜色代表身份或对应工具结果。FibCand＝证据支持的非恶性成纤维细胞候选，MalOS＝支持恶性骨肉瘤身份，Mural＝周细胞/平滑肌，Other＝其他谱系，Unk＝身份待定；Low/High/Amb＝inferCNV低/高/模糊信号，Diploid/Aneuploid＝CopyKAT分类，NA＝工具未能评估，NotRun＝未进入该工具。CNV输入图中的HO表示留出参考。Marker图横轴为基因，纵轴为身份组；点面积为表达比例，颜色为平均log1p标准化表达。统计单位为描述性细胞/核，未作患者组间推断。CNV热图横轴为常染色体Chr，纵轴为展示细胞；每个输入角色与原簇最多抽取25个可评估细胞，所有统计仍使用全部细胞。Ref/Hold/Target分别为基线参考/留出参考/目标，Type为谱系。inferCNV按100基因窗口显示log2偏离，CopyKAT按20个相邻基因组bin显示中位偏离，色标在方法内固定，两个工具色值不能直接当作相同DNA拷贝数。展示细胞数和参数见CNV_plot_manifest.tsv；未能评估细胞不画入数值热图，但保留在身份证据表。'
figcap+=' SCEVAN的Normal/Tumor是工具分类，不等同于最终身份；BadQC表示该方法未通过参考质量检查。Check为未用于该版本基线的检查细胞。SCEVAN热图同样按100基因窗口展示。review_figures保留初次复核快照，最后5个参考排除见reference_exclusions.tsv，最终身份以identity图表为准。'
report=f'''# 骨肉瘤肺转移 CAF 身份鉴定阶段报告

生成时间：{datetime.datetime.now().astimezone().isoformat()}。

## 双方法证据支持{strict}个严格CAF候选

复用 GSE270231 七个样本的已有处理对象，共28,542个细胞/核。inferCNV和原始CopyKAT的14个基础任务完成标记全部存在；另有{extra_completed}个平衡参考或独立SCEVAN核查任务实际完成。最终严格 CAF 候选 **{strict}** 个，扩展候选 **{expanded}** 个；扩展集合包含严格集合。原始CopyKAT与inferCNV方向冲突 **{raw_conflicts}** 个细胞；采用经过参考质检的独立方法后，联合可评估但仍不一致或模糊的细胞 **{discord}** 个。原始、敏感性和替代方法结果均保留，不覆盖矛盾证据。

这里的 CAF 集合指肿瘤背景下、经谱系与非恶性证据筛选的成纤维细胞候选，不将其直接命名为转移特异“新亚群”。当前全部样本来自 B2，未开展 A/B1 临床风险比较。

最终上下文复核将S0058 C17中的30个初始扩展候选移出CAF集合：该簇具有CTAG2/PAGE2B/PAGE5等癌睾抗原程序，稀疏单细胞的低CNV信号不足以解除这一身份疑点。原始记录和细胞保留，严格集合数量未改变。

{table(ss)}

严格集合涉及样本：{', '.join(s.loc[s.caf_strict>0,'short']) or '无'}。三个 SC069/SC072/SC073 的患者映射仍未知，不能把7个样本当作7位独立患者，也不能将三个 UNCONFIRMED 合并成同一患者。

身份总数：{'; '.join(f'{k}={v}' for k,v in counts.items())}。

本轮用于身份判定的第二方法见下表，第一方法始终为inferCNV。选择依据是参考质量和明确的适用性问题，不是CAF数量。

{table(selection[['short','selected_source','reference_check_pass']])}

## 接管检查确认已有结果可复用

完整包大小372,577,149字节，SHA256为 baa3d99c454ec93834ff55720699b79f4651409a187fe44b4fb62cf3546850c0。Windows与服务器一致。七个正式 annotated.h5ad 可读取，细胞数逐一符合交接基准；每个对象33,694个RNA基因，counts为非负整数，与记录的UMI总数一致。标准化表达与log1p(每细胞1万UMI归一化)一致；正式h5ad中的复核注释与导出TSV一致。临时review_tmp文件不作为主输入。

已有质控、Scrublet、聚类和UMAP均复用。原始1,470个FibLike只是历史候选标签。此次依据已有簇marker增加{len(changes)}项谱系复核记录；例如S0217的C1/C8恢复为T/NK及髓系候选，S0058的C17因癌睾抗原信号转为待定间质，SC069的C5保留为身份待定的增殖间质。原标签未覆盖。

Scrublet原自动阈值保留：S0217约0.150，S0218约0.475；后者高于全部观测分数，解释了零预测双细胞。结合混合谱系与分数尾部标记问题细胞，没有强行删除固定比例。

## 基因坐标采用与实际特征完全匹配的 GRCh38 注释

七个对象的33,694个Ensembl ID和基因符号全部与GRCh38/Ensembl84吻合，其中1,551个ID未匹配GRCh37/Ensembl75。inferCNV使用官方Ensembl84 GTF生成的常染色体位置表；CopyKAT使用其hg38内置注释并按稳定Ensembl ID匹配，重复符号和歧义映射单独记录。没有使用HVG或整合矩阵进行CNV。

GEO提交信息中的hg19原值仍保留。上述证据确认的是表达特征注释及CNV坐标匹配；原FASTQ实际比对命令未取得，不宣称已经重比对或完整还原当初参考索引。

## 参考与目标集合

主分析输入{int(inputs.run_cells.sum())}个细胞，包含广泛间质目标、破骨细胞、其他间质混合信号以及同样本参考和留出检查细胞；并不限于原FibLike。全部原始28,542个保留细胞都进入身份证据总表，未运行CNV的细胞标记not_run。

参考细胞先以已复核免疫/内皮谱系和多基因检测筛选，再排除明显间质、破骨、上皮混合及低复杂度/高双细胞分数细胞；参考最低500基因和750UMI。每样本免疫/内皮候选分别最多400/160个，约20%留出，不参与基线。再次排除5个强混合谱系参考，将其保留为CNV目标。所有样本均使用同样本参考，未借用PDX或其他患者作为主基线。

## CNV方法与证据整合

inferCNV 1.26.0：cutoff=0.1，101基因平滑窗口，denoise=True，HMM=False，4线程，同样本免疫与内皮分组基线。细胞信号为去噪表达相对值的平均绝对log2偏离；阈值取0.05、留出参考99百分位、留出参考中位数加5倍未缩放MAD（MAD下限0.005）中的最大值。高信号还要求至少2条染色体绝对中位log2偏离≥0.1，或1条≥0.25。它是描述性CNV证据，不是临床诊断或患者推断检验。

CopyKAT 1.2.5，源码commit 01924c0999bb660064da8d3e1561615ac6efb76f：ngene.chr=5，win.size=25，KS.cut=0.1，LOW.DR=0.05、请求UP.DR=0.1，已知同样本正常参考，4线程。工具会在低基因覆盖时自动调整UP.DR，实际提示保留在各样本日志。外部随机种子为270231；此版本函数内部固定重设为1234。genome参数名hg20是该工具的人hg38注释选项，不表示另一个人基因组版本。

共同运行环境为复用的inferCNV专用R 4.5.3；默认R 4.5.2缺少所需工具。CopyKAT及补充依赖只安装到本项目Rlib。初始环境准备的失败记录保留，最终工具成功状态以各样本COMPLETE、cell_results.tsv和resource_usage.txt为准。

参考检查要求每工具至少10个可评估留出细胞，其中≥90%为对应低信号/非恶性分类；并分别检查免疫与内皮基线，每类至少10个可评估细胞且≥90%为正常分类。增加分谱系检查是因为S0058的总体参考结果掩盖了明显的内皮异常分类。未通过检查或失败的样本方法记录：{'; '.join(qbad.short.astype(str)+'/'+qbad.method.astype(str)) or '无'}。参考检查通过并不自动证明工具对全部间质谱系具有同等特异性。

S0058及SC072因持续的间质/内皮冲突，进行了等数量免疫/内皮参考的CopyKAT敏感性核查，并以SCEVAN 1.0.3作为第三种独立方法评估。SCEVAN使用官方preprocessingMtx和classifyTumorCells分类步骤，GRCh38/EnsDb.Hsapiens.v86内置坐标，ngenes_chr=5，perc_genes=0.1，beta_vega=0.5，4线程，保留同样本平衡参考及原留出检查细胞。findConfident=False、FIXED_NORMAL_CELLS=False，不用成纤维签名预设正常身份，不强制其余细胞为肿瘤，不执行不必要的后续亚克隆分析。内部预处理/分类种子分别为123和1，参考抽样种子270231。逐样本完成/失败状态见对应目录。

严格CAF候选要求：复核FibLike谱系、至少2个成纤维模块基因、inferCNV与质检选用的独立方法低信号/非恶性一致、无强竞争谱系/混合/癌睾抗原证据、至少500基因，并通过参考检查。若所在簇联合可评估细胞中≥20%呈双工具高信号，则不纳入严格集合，以防肿瘤簇中的低信号细胞漏判。扩展集合用于敏感性分析，允许部分谱系边界或单工具不可评估，但不纳入选用方法出现相反高CNV证据的细胞。被判定不适用的原CopyKAT分类不作为最终否决证据，其冲突仍在原始字段和问题列中保留。

若既有前50位簇marker包含至少2个预定义癌睾抗原基因，该簇整体排除出严格及扩展CAF集合；细胞保留原工具结果，证据不足者仍为身份待定。这一上下文规则防止把疑似肿瘤程序簇中的低复杂度或掉零细胞当作成纤维细胞。

双工具高信号还需符合间质肿瘤背景，才能支持恶性骨肉瘤身份。周细胞/平滑肌、破骨、免疫、内皮和上皮证据单独保留。工具运行失败或未能评估不视为CNV阴性。不一致细胞不强制归类。

## 主要不一致证据与后续边界

SC072显示明显原始工具不一致：CopyKAT对多种间质/上皮细胞给出非整倍体分类，而inferCNV没有相应大范围信号。S0058原CopyKAT还将大量复核内皮及参考细胞判为非整倍体。已检查热图并开展针对性核查，不能仅凭其中一个结果确认该群为恶性或CAF。完整原始值、替代方法值、最终采用方法及未解决问题均在细胞表中；簇级一致性见cluster_CNV_concordance.tsv，参考情况见tool_reference_QC.tsv。

当前数据支持的是B2内细胞身份筛选，不能回答未发生肺转移患者与将发生/已发生肺转移患者的CAF差异。CNV为表达推断，低信号仍需结合谱系；未明确的细胞不反向否定已经由多项证据支持的细胞集合。新CAF状态、成骨/矿化机制与活细胞分选panel尚未完成。

## 图注与基因证据

{figcap} Marker图括号内为该身份组细胞数，未出现的身份组不绘制空行；Expr %为检测比例图例。FibCand颜色表示身份候选集合，并不等同于严格集合；严格与扩展集合的具体细胞位置见各样本CAF_sets图。

每张图的具体坐标、分组、统计单位和主要结果见figure_captions.tsv。每样本身份marker图实际展示DCN、LUM、PDGFRA、FAP、基质/收缩、成骨/矿化及排除谱系标志；数值见identity_marker_expression.tsv。原始完整簇marker与复核marker比例表保留用于后续状态识别。CD36未预设为CAF阳性标志，ACTA2/COL1A1等不作为活细胞表面分选标志。

## 下一步

先依据本轮严格集合的样本/患者分布开展CAF再聚类，并以扩展集合进行必要敏感性分析；检查患者主导、技术差异，再分析肌成纤维、成骨、矿化和基质重塑程序与真实候选基因。临床分组不足时只做B2内探索，不用细胞数扩大患者样本量。

待补充：SC069/SC072/SC073患者映射；骨原发样本的肺转移史及相对采样时间；骨/肺准确间隔和取材前治疗；流式仪器配置及抗体可得性。ScPCA临床补充与主线并行，PDX和PDX衍生模型不进入人患者CAF主分析。

## 文件与复现

细胞级总表：identity/all_cell_identity_evidence.tsv.gz。严格/扩展清单：identity/caf_strict_cells.tsv、caf_extended_cells.tsv。大型表达对象单独位于identity/objects。原始七个h5ad保留在OS_CAF_phase2/results，未覆盖。

代码、环境、参数、来源、图表和日志摘要随轻量包回传。所有脚本写明SYNTAX_CHECK、RUN_COMMAND、EXPECTED_SUCCESS_MARKER。原始CNV大矩阵/RDS、GTF、完整运行日志保留在服务器。Methods和解释以实际参数、源文件及完整日志为准。

参考入口：[GEO](https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270231)、[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11739783/)、[inferCNV](https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV)、[CopyKAT](https://github.com/navinlabcode/copykat)、[SCEVAN](https://github.com/AntonioDeFalco/SCEVAN)、[Ensembl84 GTF](https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz)。
'''
(delivery/'阶段报告.md').write_text(report,encoding='utf-8')
readme='''# 执行与传输

本次不是仅准备代码：以14个工具COMPLETE和identity/COMPLETE确认实际计算状态。

服务器项目：/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF
阶段目录：/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914

原始阶段包上传命令（本次已完成）：
```powershell
scp "C:\\Users\\31937\\Downloads\\OS_CAF_phase2_bundle.tar.gz" sgs10:/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/incoming/
```

轻量回传与表达对象下载：
```powershell
scp sgs10:/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/delivery/OS_CAF_phase3_review.zip "C:\\Users\\31937\\Downloads\\OS_CAF_phase3_review.zip"
scp sgs10:/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/identity/objects/caf_strict.h5ad "C:\\Users\\31937\\Downloads\\caf_strict.h5ad"
scp sgs10:/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/identity/objects/caf_extended.h5ad "C:\\Users\\31937\\Downloads\\caf_extended.h5ad"
```
若某集合为空，清单仍保存但不会产生对应h5ad，以sample_identity_summary.tsv为准。

单图下载查看示例：
```powershell
scp sgs10:/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914/identity/figures/S0058_identity_UMAP.png "C:\\Users\\31937\\Downloads\\S0058_identity_UMAP.png"
```

检查状态：
```bash
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
cat "$ROOT/identity/sample_identity_summary.tsv"
cat "$ROOT/logs/cohort_infercnv_0.log"
cat "$ROOT/logs/cohort_copykat_1.log"
cat "$ROOT/logs/cohort_copykat_2.log"
```

复现顺序：takeover_audit.py核验对象；review_inputs.py复用counts与UMAP生成证据；refine_reference.py移除强混合参考；coordinate_audit.py核实注释版本；run_cnv_job.sh运行单个工具样本（已完成会跳过）；integrate_identity.py在全部任务完成后生成身份表及集合；build_delivery.py生成回传包。review_inputs和坐标下载/环境准备无需重复运行已有完成结果。不得运行phase2/run_server.sh start覆盖已有结果。

完成等待使用wait_and_integrate.sh，一次服务器阻塞循环每900秒检查进程存在；运行中日志用于判断是否需要干预。每个工具任务有job.pid、job.started、job.exit_code、resource_usage.txt、status.json和COMPLETE/FAILED。
'''
(delivery/'运行与传输.md').write_text(readme,encoding='utf-8')
items=[]
def add(f,arc):
    if f.is_file():items.append((f,arc))
for f in (p/'code').iterdir():add(f,'code/'+f.name)
for f in delivery.iterdir():
    if f.suffix in ['.md','.tsv']:add(f,f.name)
for f in idd.rglob('*'):
    if f.is_file() and 'objects' not in f.parts:add(f,'identity/'+f.relative_to(idd).as_posix())
for f in (p/'review').iterdir():
    if f.is_file():add(f,'review/'+f.name)
for f in (p/'figures').glob('*.png'):add(f,'review_figures/'+f.name)
for method in ['infercnv','copykat','copykat_balanced','scevan']:
    for short in SHORTS:
        d=p/method/short
        for name in ['cell_results.tsv','parameters.json','status.json','resource_usage.txt','sessionInfo.txt','COMPLETE','FAILED','job.exit_code','job.started','reference_cells.tsv','balanced_reference_cells.tsv']:
            add(d/name,method+'/'+short+'/'+name)
        # Original vendor plots remain on the server; curated numeric heatmaps are in identity/figures.
for name in ['coordinate_audit.tsv','coordinate_provenance.json','reference_download.json','copykat_commit.json','SCEVAN_commit.json','yaGST_commit.json']:
    add(p/'resources'/name,'provenance/'+name)
for name in ['manifest.tsv','review_map.tsv','phase2.py','finalize_review.py']:add(a.root/'code'/name,'phase2_code/'+name)
for name in ['sample_summary.tsv','lineage_counts_reviewed.tsv','cnv_candidates_reviewed.tsv']:add(a.root/'results'/name,'phase2_summary/'+name)
for row in manifest.itertuples():
    for name in ['cluster_markers.tsv','marker_evidence_reviewed.tsv','cluster_identity_evidence.tsv','complete.json','qc_thresholds.json']:
        add(a.root/'results'/row.sample/name,'phase2_summary/'+row.short+'/'+name)
for f in (p/'takeover_audit').glob('*'):add(f,'takeover_audit/'+f.name)
logdir=delivery/'log_summaries';logdir.mkdir(exist_ok=True)
for f in (p/'logs').iterdir():
    if not f.is_file():continue
    if f.stat().st_size<250000:add(f,'logs/'+f.name)
    else:
        tail=logdir/(f.name+'.tail.txt');lines=f.read_text(errors='replace').splitlines();tail.write_text(f'Full log remains on server: {f}\nSHA256: {hashlib.sha256(f.read_bytes()).hexdigest()}\nTotal lines: {len(lines)}; last 300 follow\n'+'\n'.join(lines[-300:]),encoding='utf-8');add(tail,'logs/'+tail.name)
filemanifest=[{'path':arc,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f,arc in items]
mf=delivery/'package_manifest.json';mf.write_text(json.dumps(filemanifest,ensure_ascii=False,indent=2),encoding='utf-8');items.append((mf,mf.name))
archive=delivery/'OS_CAF_phase3_review.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for f,arc in items:z.write(f,arc)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
(delivery/'COMPLETE').write_text('PHASE3_DELIVERY_COMPLETE\n')
print('PHASE3_DELIVERY_COMPLETE',strict,expanded,archive.stat().st_size)
