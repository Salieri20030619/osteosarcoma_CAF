import sys
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase46_spatial_functional_20261001/code')
from common import *
import gzip,io,anndata as ad
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
Q=P/'phase48_clinical_programs_delivery_20261001';H=Q/'figures/histology';H.mkdir(exist_ok=True);z=pd.read_csv(O/'tables/all_spot_RNA_and_components.tsv.gz',sep='\t');rows=[]
for sm,t in z.groupby('sample',sort=True):
 a=ad.read_h5ad(P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad',backed='r');f=Path(a.uns['provenance']['lowres_image_path']);sc=float(a.uns['spatial_scalefactors']['tissue_lowres_scalef']);diam=float(a.uns['spatial_scalefactors']['spot_diameter_fullres'])*sc;raw=gzip.open(f,'rb').read();im=Image.open(io.BytesIO(raw));im.load();width,height=im.size;xy=a.obsm['spatial']*sc;aligned=(xy[:,0]>=0)&(xy[:,0]<width)&(xy[:,1]>=0)&(xy[:,1]<height);assert aligned.all();assert list(a.obs_names)==t.spot_id.tolist();assert np.allclose(xy,t[['x','y']].to_numpy()*sc)
 rows.append(dict(sample=sm,image_path=str(f),image_SHA256=sha(f),decoded_image_SHA256=hashlib.sha256(raw).hexdigest(),width=width,height=height,lowres_scalefactor=sc,lowres_spot_diameter_pixels=diam,spots=len(t),coordinates_in_image=int(aligned.sum()),image_decode_success=True,coordinate_alignment_accepted=True,tumor_mask_accepted=False,geometry_endpoint='unassessable_no_independent_tumor_boundary',reason='low-resolution tissue image available; not independently annotated cancer/stroma/CD8 cell masks'))
 fig,axs=plt.subplots(1,4,figsize=(13,4));arr=np.asarray(im.convert('RGB'))
 for ax in axs:ax.imshow(arr,origin='upper');ax.set_xlim(0,width);ax.set_ylim(height,0);ax.set_aspect('equal');ax.set_axis_off()
 axs[0].set_title('Tissue lowres')
 for ax,name,label in zip(axs[1:],['F3_SELECTED_REMODELING','CD8_RNA','cancer'],['F3 remodelling RNA','CD8 RNA','Cancer signed IC']):
  v=ax.scatter(xy[:,0],xy[:,1],c=t[name],s=6,alpha=.65,cmap='viridis',rasterized=True);ax.set_title(label);fig.colorbar(v,ax=ax,shrink=.55,pad=.02)
 fig.suptitle(sm+'  '+t.patient_id.iloc[0]+'  low-resolution tissue with RNA overlay');fig.tight_layout();fig.savefig(H/f'{sm}_tissue_RNA_overlay.png',dpi=180,bbox_inches='tight');plt.close(fig)
 if sm in ['PT1','MR1']:im.convert('RGB').save(H/f'{sm}_original_lowres.png')
 a.file.close()
frame=pd.DataFrame(rows);frame.to_csv(Q/'tables/histology_image_coordinate_audit.tsv',sep='\t',index=False);qa={'images_available':len(rows),'images_decoded_and_aligned':int(frame.coordinate_alignment_accepted.sum()),'spots_in_bounds':int(frame.coordinates_in_image.sum()),'independent_tumor_masks_accepted':0,'prior_statement_corrected':'only embedded AnnData images were checked; downloaded gzipped PNG provenance paths contain actual images','geometry_endpoint_status':'unassessable due independent tumor boundary/region labels absent, not due image absence','existing_RNA_statistics_unchanged':True};dump(Q/'delivery/HISTOLOGY_IMAGE_AUDIT.json',qa);(Q/'HISTOLOGY_AUDIT_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
