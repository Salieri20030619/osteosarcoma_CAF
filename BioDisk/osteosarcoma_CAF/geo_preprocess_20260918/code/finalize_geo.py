from pathlib import Path
import json,datetime,html,zipfile
import pandas as pd
from PIL import Image,ImageOps,ImageDraw
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_preprocess_20260918')
assert (P/'review/OUTPUTS_VERIFIED').exists()
s=pd.read_csv(P/'review/geo_processing_summary.tsv',sep='\t')
n=pd.read_csv(P/'review/cnv_readiness.tsv',sep='\t')
s=s.merge(n,on='sample',suffixes=('','_review'),validate='one_to_one')
p=pd.read_csv(P/'review/cluster_review_packets.tsv',sep='\t',dtype={'cluster':str}).fillna('')
s['all_CAF_calls_pending_CNV']=True
s.to_csv(P/'delivery/GEO_sample_acceptance_and_QC.tsv',sep='\t',index=False)
p.to_csv(P/'delivery/GEO_cluster_lineage_review.tsv',sep='\t',index=False)
for i in range(0,len(s),4):
 batch=s.iloc[i:i+4];canvas=Image.new('RGB',(2200,650*len(batch)), 'white');d=ImageDraw.Draw(canvas)
 for j,r in enumerate(batch.to_dict('records')):
  d.text((10,j*650+5),r['sample'],fill='black')
  for k,fn in enumerate(['QC.png','UMAP.png']):
   im=Image.open(P/'results'/r['sample']/fn).convert('RGB');im.thumbnail((1080,600));canvas.paste(im,(k*1100,j*650+30))
 canvas.save(P/'review'/('QA_QC_UMAP_'+str(i//4+1)+'.png'))
for i in range(0,len(s),4):
 batch=s.iloc[i:i+4];canvas=Image.new('RGB',(2600,1000*len(batch)), 'white');d=ImageDraw.Draw(canvas)
 for j,r in enumerate(batch.to_dict('records')):
  d.text((10,j*1000+5),r['sample'],fill='black')
  im=Image.open(P/'results'/r['sample']/'Markers.png').convert('RGB');im.thumbnail((2550,950));canvas.paste(im,(0,j*1000+30))
 canvas.save(P/'review'/('QA_markers_'+str(i//4+1)+'.png'))
print(s[['sample','input_cells','retained','doublets','proposed_fibroblast_cells','eligible_immune_reference','eligible_endothelial_reference']].to_string(index=False))
print('GEO_DELIVERY_TABLES_AND_QA_READY',flush=True)
