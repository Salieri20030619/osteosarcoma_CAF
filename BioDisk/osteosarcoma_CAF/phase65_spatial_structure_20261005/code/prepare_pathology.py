from pathlib import Path
import pandas as pd,json,gzip,io,hashlib,zipfile
from PIL import Image
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase65_spatial_structure_20261005';D=O/'pathology_review'
for sub in ['images','annotations']:(D/sub).mkdir(parents=True,exist_ok=True)
mf=pd.read_csv(O/'resources/SPATIAL_INPUT_MANIFEST.tsv',sep='\t');scores=pd.read_csv(O/'tables/spatial_gene_and_program_scores.tsv.gz',sep='\t',usecols=['sample','spot_id','x','y','patient_id','site_group'],low_memory=False)
links=[];public=[];coords=[]
for i,r in enumerate(mf.sort_values('sample').itertuples(),1):
 code=f'H{i:03d}';src=Path(r.image_path);hi=Path(str(src).replace('_lowres_','_hires_'));chosen=hi if hi.exists() else src;resolution='hires' if hi.exists() else 'lowres'
 data=gzip.decompress(chosen.read_bytes());im=Image.open(io.BytesIO(data));out=D/'images'/f'{code}.png';out.write_bytes(data)
 prefix=src.name.split('_tissue_')[0]
 scale_file=src.parent/(prefix+'_scalefactors_json.json.gz');sc=json.loads(gzip.decompress(scale_file.read_bytes()))
 scale=float(sc['tissue_'+resolution+'_scalef'])
 public.append({'image_id':code,'image_file':f'images/{code}.png','width_px':im.width,'height_px':im.height,'available_resolution':resolution,'annotation_status':'unannotated'})
 t=scores[scores['sample']==r.sample]
 links.append({'image_id':code,'sample':r.sample,'patient_id':t.patient_id.iloc[0],'site_group':t.site_group.iloc[0],'source_image':str(chosen),'fullres_to_image_scale':scale,'image_sha256':hashlib.sha256(data).hexdigest()})
 q=t[['spot_id','x','y']].copy();q['image_id']=code;q['image_x']=q.x*scale;q['image_y']=q.y*scale
 assert q.image_x.between(0,im.width,inclusive='left').all() and q.image_y.between(0,im.height,inclusive='left').all()
 coords.append(q)
 (D/'annotations'/f'{code}.geojson').write_text(json.dumps({'type':'FeatureCollection','image_id':code,'coordinate_system':'image_pixels_top_left_y_down','features':[]},indent=2))
pd.DataFrame(public).to_csv(D/'image_index.tsv',sep='\t',index=False)
pd.DataFrame(links).to_csv(O/'resources/PATHOLOGY_PRIVATE_KEY.tsv',sep='\t',index=False)
pd.concat(coords,ignore_index=True).to_csv(O/'tables/pathology_spot_coordinate_registration.tsv.gz',sep='\t',index=False)
instructions="""独立病理复核准备包（未标注）

本包仅含原始H&E图像的匿名编号副本、尺寸目录及空轮廓模板，不含胶原或细胞程序分数。请独立于RNA结果评阅；所有空模板表示“尚未标注”，不是阴性，也不是正常组织。

当前图像分辨率见image_index.tsv。本批可用图像以lowres预览为主；如果分辨率不足以可靠区分活肿瘤、骨样基质、坏死和血管，应标记“无法判断”，并索取对应原始高分辨率图或切片。不得仅凭低分辨率概览生成被视为金标准的轮廓。图像本身不提供胶原特异蛋白或纤维染色证据。

建议区域类别：viable_tumor（可辨识活肿瘤）、necrosis（坏死）、osteoid_cartilage（骨样/软骨基质）、vascular_region（可辨识血管区域）、non_tumor_stroma（非肿瘤间质）、uncertain（不确定）。不要求所有类别均出现，不强迫模糊区域归类。
标注时记录评阅者、日期、置信度、是否使用高分辨率配准图及无法判断的原因。建议由第二评阅者核对关键肿瘤边界，记录分歧。

每个annotations/Hxxx.geojson是空FeatureCollection。多边形坐标基于对应PNG：左上为(0,0)，x向右、y向下；首尾坐标闭合。每个Feature的properties包含label、reviewer、reviewed_at、confidence及comment。可使用支持图像坐标的工具标注后转换；不要填写经纬度或地理投影。
回收后需先检查匿名编号、图像哈希、坐标范围和图像配准。仅在真实人工轮廓通过复核后，才可进入肿瘤边缘距离、周界覆盖与连续性统计。当前包包含0份已接受的独立肿瘤轮廓。
"""
(D/'README_病理复核说明.txt').write_text(instructions,encoding='utf-8')
q={'images':len(public),'lowres':sum(x['available_resolution']=='lowres' for x in public),'hires':sum(x['available_resolution']=='hires' for x in public),'accepted_masks':0,'blank_templates':26,'all_coordinates_in_bounds':True,'anonymous_key_outside_blinded_package':True,'purpose':'Pathology review preparation only; no inferred independent contour'}
(O/'QA/PATHOLOGY_PREPARATION_QA.json').write_text(json.dumps(q,indent=2))
zpath=O/'delivery/CAF_pathology_review_preparation_20261005.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for p in D.rglob('*'):
  if p.is_file():z.write(p,arcname=str(p.relative_to(O)))
print(json.dumps(q),flush=True)
