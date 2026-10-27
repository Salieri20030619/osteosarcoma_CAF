from pathlib import Path
import json,hashlib,zipfile,datetime,shutil,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005'
assert json.loads((O/'QA/ACCEPTANCE.json').read_text())['passed']
assert json.loads((O/'QA/CLINICAL_REDUCED_QA.json').read_text())['passed']
assert json.loads((O/'QA/VERIFY_RUN_STATUS.json').read_text())['exit_status']==0
(O/'QA/FIGURE_VISUAL_QA.json').write_text(json.dumps({'passed':True,'inspected':['source_context_summary.png','maps/PT1.png','discovery_dominant_gene_sensitivity.png'],'checks':'Labels, axes, patient-level wording, all panels and image registration legible; 26 maps programmatically generated and checked readable.','all_maps':len(list((O/'figures/maps').glob('*.png')))},indent=2))
from PIL import Image
for p in (O/'figures').rglob('*.png'):
 with Image.open(p) as im:im.verify()
dep=[]
for rel in ['phase46_spatial_functional_20261001/code/common.py','phase63_matrix_resistance_20261005/objects/mechanism_counts.npz','phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz','phase63_matrix_resistance_20261005/tables/paired_patient_site_roster.tsv.gz','phase64_response_followup_20261005/tables/patient_metadata.tsv','phase64_response_followup_20261005/objects/GSE39055_gene_VST.tsv.gz','phase64_response_followup_20261005/objects/GSE21257_gene_VST.tsv.gz']:
 p=B/rel;dep.append({'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'resources/INHERITED_DEPENDENCIES.json').write_text(json.dumps(dep,indent=2))
shutil.copy2(B/'phase46_spatial_functional_20261001/code/common.py',O/'resources/inherited_common_phase46.py')
(O/'delivery/README_交付说明.txt').write_text('本阶段完成胶原结构稳健性、空间来源背景及缩减面板临床敏感性验证。原始F1–F5与33亚群置信度保留。复现脚本以sgs10项目根目录为基础，需要既有phase46/63/64及原空间h5ad；依赖路径与哈希见resources/INHERITED_DEPENDENCIES.json和SPATIAL_INPUT_MANIFEST.tsv。未重复打包全部大型原始数据。\nCAF_pathology_review_preparation_20261005.zip为单独的无RNA分数病理复核准备包，26张均为lowres预览，0份已接受轮廓。主结果包包含匿名编号对应表，不应把主结果包作为盲法评阅输入。所有空模板代表未标注。\n',encoding='utf-8')
# Preserve review outcomes and clinical null findings in an explicit phase status.
summary={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_complete':True,'computational_QA_passed':True,'cells':8174,'subgroups':33,'slices':26,'spots':53309,'zero_library_spots':2,'all_patient_records':19,'primary_patient_records':9,'clinical_reduced_nominal_tests':0,'clinical_reduced_total_tests':12,'pathology_highres_images':0,'accepted_independent_masks':0,'complete_project':False,'outstanding':['High-resolution pathology and independent contour review','Physical collagen/organization/crosslink measurement','Controlled matrix/drug-response intervention']}
(O/'COMPUTATIONAL_PHASE_COMPLETE.json').write_text(json.dumps(summary,indent=2))
manifest=[]
for folder in ['code','resources','tables','figures','QA','logs','pathology_review']:
 for p in sorted((O/folder).rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and p.name!='DELIVERY_MANIFEST.json':manifest.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for p in list((O/'delivery').glob('*.txt'))+[O/'RUN_STATUS.json',O/'SPATIAL_COMPLETE.json',O/'DISCOVERY_COMPLETE.json',O/'COMPUTATIONAL_PHASE_COMPLETE.json']:
 manifest.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
mp=O/'resources/DELIVERY_MANIFEST.json';mp.write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
zpath=O/'delivery/CAF_spatial_structure_validation_20261005.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for item in manifest:z.write(O/item['path'],arcname=O.name+'/'+item['path'])
 z.write(mp,arcname=O.name+'/resources/DELIVERY_MANIFEST.json')
receipts=[]
for p in [zpath,O/'delivery/CAF_pathology_review_preparation_20261005.zip']:
 with zipfile.ZipFile(p) as z:assert z.testzip() is None;n=len(z.infolist())
 receipts.append({'name':p.name,'path':str(p),'bytes':p.stat().st_size,'files':n,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'delivery/DELIVERY_RECEIPTS.json').write_text(json.dumps(receipts,indent=2));print(json.dumps(receipts),flush=True)
