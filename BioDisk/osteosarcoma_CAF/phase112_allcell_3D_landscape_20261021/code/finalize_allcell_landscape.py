from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile,shutil
import pandas as pd
P=Path('C:/Users/31937/CAF_paper_figures_20261021/allcell_3d')
O=Path('C:/Users/31937/Downloads/CAF_allcell_3D_Fig1b_1d_20261021')
mf=P/'resources/FIGURE_MANIFEST.json';m=json.loads(mf.read_text(encoding='utf-8'))
browser=json.loads((P/'resources/INTERACTIVE_BROWSER_QA.json').read_text(encoding='utf-8'));assert browser['status']=='PASS'
assert m['static_figures']==12 and len(m['interactive_files'])==2
files=[]
for f in m['figures']:
 for e in f['exports']:
  p=Path(e['file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==e['sha256'];files.append(p)
for f in m['interactive_files']:
 p=Path(f['file']);assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256'];files.append(p)
methods=(P/'图注与方法_草稿.txt').read_text(encoding='utf-8')
diag=pd.read_csv(P/'embedding_diagnostic.tsv',sep='\t')
methods+='\n本次邻域诊断结果（固定2000细胞子集，k=15）：\n'+diag.to_csv(sep='\t',index=False)
methods+='\n交付图目：\n'
for i,f in enumerate(m['figures'],1):methods+=f"{i:02d}. {f['name']}\n"
methods+='\n另附两个来源各自的离线交互HTML，已以无界面Edge核查全部细胞数、四个部位筛选及相机旋转。静态图共12张，每张PNG/TIFF均为600 DPI，已逐张检查文字、图例、坐标、颜色、计数与版式。\n'
caption=O/'全细胞3D_组成图_图注与方法.txt';caption.write_text(methods,encoding='utf-8');files.append(caption)
readme=O/'先读我_图件索引.txt'
readme.write_text('''全细胞景观与组成图（2026-10-21）

第一批范围：SCPCP000017、SCPCP000023，合计345219个上游called细胞。采用全细胞分母和冻结复核注释，不是此前6633个CAF候选池。NCI不包含在本套。

Fig1b：每个来源各有总体、骨、肺3张三维UMAP；同一来源共用坐标、范围和视角。另附每个来源原始二维UMAP对照，共8张。
Fig1d：每个来源各有患者等权的部位汇总组成图、逐患者组成图，共4张。
每张PNG与TIFF为同图两种格式，全部600 DPI。

HTML：双击对应来源的interactive_3D.html可离线旋转与缩放；上方按钮可筛选总体、骨、肺、未知/混合部位。包含全部细胞，首次加载可能需要等待。HTML属于交互附件，无DPI概念。

骨/肺是解剖病灶分组，不是PD-1治疗响应分组。未知部位、未定谱系和混合类型均保留。各来源存在同时贡献骨肺样本的患者，组成图为描述性展示，不把细胞或患者-部位记录当作额外独立患者。

请同时阅读“全细胞3D_组成图_图注与方法.txt”。数据、脚本、坐标和核查记录归档于sgs10：
/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase112_allcell_3D_landscape_20261021
''',encoding='utf-8');files.append(readme)
now=datetime.now(timezone.utc).isoformat();m.update(status='PASS',visual_QA='All 12 static previews inspected individually: axes/legends/counts readable, no clipping; interactive all/rotated screenshots inspected',interactive_browser_QA=browser,completed_at=now)
mf.write_text(json.dumps(m,ensure_ascii=False,indent=2),encoding='utf-8')
archive=O/'CAF_Fig1b三维UMAP_Fig1d组成图_600DPI.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,p.name)
with zipfile.ZipFile(archive) as z:assert z.testzip() is None and len(z.namelist())==28
receipt=dict(status='completed_verified',phase='phase112_allcell_3D_landscape_20261021',completed_at=now,sources=['SCPCP000017','SCPCP000023'],all_cells=345219,static_figures=12,static_dpi=600,interactive_html=2,all_cells_displayed=True,composition_denominator='all upstream called cells per patient/site',annotations='Phase55 frozen reviewed expression lineages',clinical_groups=['Bone','Lung','Unknown_or_mixed'],paired_patients_preserved=True,statistical_inference='none; descriptive composition',original_analysis_files_modified=False,local_delivery=str(O),archive=str(archive),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_bytes=archive.stat().st_size,zip_members=28,zip_crc='PASS',visual_QA='PASS',interactive_QA='PASS')
(P/'resources/FINAL_STATUS.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(caption,P/'resources'/caption.name);shutil.copy2(readme,P/'resources'/readme.name)
print(json.dumps(receipt,ensure_ascii=False))
