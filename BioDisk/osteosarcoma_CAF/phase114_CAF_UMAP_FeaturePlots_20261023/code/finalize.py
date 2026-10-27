from pathlib import Path
import json,hashlib,zipfile,datetime,shutil
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase114_CAF_UMAP_FeaturePlots_20261023'
def dump(p,x): p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
qa=json.loads((O/'resources/DATA_QA.json').read_text());assert qa['status']=='PASS'
manifest=json.loads((O/'resources/FIGURE_MANIFEST.json').read_text());assert len(manifest)==21
for m in manifest:
 for f in m['files']:assert sha(O/'figures'/f['name'])==f['sha256']
caption='''CAF 候选 UMAP 与六基因 FeaturePlot（Phase114，2026-10-23）

本批共21张独立图：3张Fig1g候选UMAP，18张Fig1j基因FeaturePlot。每张提供PNG及无损TIFF，4440×3720像素，600 DPI。PNG和TIFF像素已逐一核对一致。

图注
Fig1g：在既有全细胞UMAP中显示当前候选池的细胞位置。SCPCP000017包含3269个细胞、15个来源内患者记录；SCPCP000023包含1834个细胞、7个来源内患者记录；NCI包含1530个作者注释CAF细胞、17个来源内患者记录。原始两个来源沿用冻结的F1-F5候选类别。NCI沿用作者CAF注释，未强行分配F1-F5。点代表细胞；颜色代表既有类别。
Fig1j：在对应来源相同的坐标、细胞名单和坐标范围上展示PNKD、TOR1AIP1、GRID2、TAF10、BEST1、RGS2-AS1的RNA表达。颜色为log1p(CP10k)，即log1p(原始RNA计数/该细胞完整文库计数×10000)。同一基因在三个来源统一使用0到全体显示细胞最大值的色标，不截断极值。零值先绘制，阳性细胞按表达从低到高叠加。图中RNA detected为原始计数大于0的细胞数，分母为该来源全部显示的候选细胞。

方法与范围
细胞名单及六基因表达直接继承已核对的Phase111 Fig1i细胞表，对应Phase103主分析扩展候选池。6633个细胞均通过唯一标识符一对一匹配，未丢失细胞。每个细胞的log1p(CP10k)由原始计数和完整文库计数重新核对。
原始来源坐标来自Phase53各来源冻结的全细胞UMAP。NCI坐标来自OsteoCAR patient_mets.qs2作者原始对象的既有UMAP（该作者对象共221271个细胞），仅提取本批1530个NCI作者CAF的位置；不拼接逐样本独立UMAP，不重新拟合坐标，不聚类或改变标签。

解释边界
三个来源的UMAP为不同坐标空间，位置不能跨来源直接比较；这些图是全细胞嵌入上的候选子集，并非重新构建的CAF专属嵌入。F1-F5为既有扩展候选分类，不代表已确认纯CAF。39为来源内患者记录数之和，不构成跨研究独立性的新增验证。NCI原始来源是单核转录组，统一色标不消除平台差异。
本批仅作RNA分布描述，无新增P值或FDR，不以细胞数充当生物学重复。嵌入未为六基因另行排除特征，表达与空间聚集不是独立机制验证。RNA检出不等于蛋白表达或流式阳性。NCI沿用作者合并对象基因覆盖记录，原始逐文库合并前特征掩码不可用。

文件索引
PNG/、TIFF/：各21张单图，文件名包含来源和基因。
本说明是图注及方法附件。细胞级数据、代码、输入哈希和验收证据保留在服务器phase114_CAF_UMAP_FeaturePlots_20261023目录。
'''
cap=O/'resources/Figure_captions_ZH.txt';cap.write_text(caption,encoding='utf-8')
visual=dict(status='PASS',reviewed_figures=21,method='All 21 previews inspected on three source contact sheets; individual NCI UMAP and SCPCP000023 BEST1 preview inspected for legends and colorbar; revised NCI legend rechecked',correction='NCI author CAF legend wrapped to prevent right-edge clipping',remaining_required_issues=[])
dump(O/'resources/VISUAL_QA.json',visual)
archive=O/'delivery/CAF_Fig1g_Fig1j_21figures_600DPI.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for m in manifest:
  for f in m['files']:
   p=O/'figures'/f['name'];z.write(p,('PNG/' if p.suffix=='.png' else 'TIFF/')+p.name)
 z.write(cap,'Figure_captions_ZH.txt')
with zipfile.ZipFile(archive) as z:assert z.testzip() is None;assert len(z.infolist())==43
status=dict(status='completed_verified_server',phase=O.name,completed_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),figures=21,dpi=600,cells=6633,sources=['SCPCP000017','SCPCP000023','NCI'],genes=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1'],numeric_QA='PASS',visual_QA='PASS',server_archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,archive_members=43,zip_crc='PASS',original_results_modified=False,local_delivery_verified=False,build_exit_code=0,all_figure_adaptation_complete=False)
dump(O/'FINAL_STATUS.json',status);dump(O/'STATUS.json',status)
print(json.dumps(status,ensure_ascii=False))
