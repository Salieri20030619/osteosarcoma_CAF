from pathlib import Path
import json,gzip,io
import pandas as pd,numpy as np,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005';z=pd.read_csv(O/'tables/spatial_spot_scores.tsv.gz',sep='\t',low_memory=False);figdir=O/'figures/spatial_maps'
files=[]
for sm,t in z.groupby('sample',sort=True):
 a=ad.read_h5ad(P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad',backed='r');im=Image.open(io.BytesIO(gzip.open(a.uns['provenance']['lowres_image_path'],'rb').read())).convert('RGB');scale=float(a.uns['spatial_scalefactors']['tissue_lowres_scalef']);xy=t[['x','y']].to_numpy()*scale;a.file.close()
 fig,axs=plt.subplots(2,3,figsize=(13,9.6));names=['cancer','fibroblast','CD8_RNA','full_program__collagen_biosynthesis','full_program__collagen_crosslinking','full_program__ECM_degradation'];titles=['Tumor signed IC','Fibroblast signed IC','CD8 RNA','Collagen biosynthesis RNA','Collagen crosslinking RNA','ECM degradation RNA']
 for ax,name,title in zip(axs.ravel(),names,titles):
  ax.imshow(im);sc=ax.scatter(xy[:,0],xy[:,1],c=t[name],s=6,cmap='viridis',alpha=.75);ax.set_xlim(0,im.width);ax.set_ylim(im.height,0);ax.set_aspect('equal');ax.set_axis_off();ax.set_title(title,pad=12,fontsize=11);fig.colorbar(sc,ax=ax,shrink=.65,pad=.025)
 fig.suptitle(sm+' | '+str(t.patient_id.iloc[0])+' | '+str(t.site_group.iloc[0]),fontsize=14,y=.985);fig.text(.015,.025,'RNA programs and signed components over tissue image. No independently annotated tumor boundary or physical collagen measurement.',fontsize=8.5);fig.subplots_adjust(left=.015,right=.985,bottom=.085,top=.91,hspace=.28,wspace=.18);path=figdir/f'{sm}.png';fig.savefig(path,dpi=150);plt.close(fig);files.append(path)
for start in range(0,len(files),8):
 fig,axs=plt.subplots(4,2,figsize=(13,19))
 for ax in axs.ravel():ax.set_axis_off()
 for ax,path in zip(axs.ravel(),files[start:start+8]):ax.imshow(Image.open(path));ax.set_title(path.stem,fontsize=12,pad=2)
 fig.subplots_adjust(left=.005,right=.995,top=.98,bottom=.005,wspace=.02,hspace=.05);fig.savefig(O/'QA'/f'maps_contact_{start//8+1}.png',dpi=120);plt.close(fig)
print('SPATIAL_FIGURES_RENDERED',len(files))
