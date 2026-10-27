from pathlib import Path
import anndata as ad,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase9_joint_20260919');a=ad.read_h5ad(O/'objects/joint_CAF_integrated.h5ad');ob=a.obs.copy();ob['site_display']=ob.site.str.lower().replace({'femur':'bone'});clusters=sorted(ob.joint_cluster.astype(str).unique(),key=int)
plt.rcParams.update({'font.size':10})
fig,axs=plt.subplots(2,3,figsize=(17,10))
for row,rep in enumerate(['X_umap_uncorrected','X_umap_joint_main']):
 for col,key in enumerate(['dataset','site_display','joint_cluster']):
  ax=axs[row,col];vals=ob[key].astype(str);cats=sorted(vals.unique(),key=int) if key=='joint_cluster' else sorted(vals.unique());palette=plt.get_cmap('tab20')
  for j,v in enumerate(cats):
   m=(vals==v).to_numpy();xy=a.obsm[rep][m];ax.scatter(xy[:,0],xy[:,1],s=3,c=[palette(j)],label=('J'+v if key=='joint_cluster' else v),rasterized=True,alpha=.8)
  ax.set_title(('Before correction' if row==0 else 'Harmony: study + technology')+' | '+key.replace('_display',''));ax.set_xticks([]);ax.set_yticks([]);ax.legend(markerscale=3,fontsize=8,loc='best',framealpha=.8)
fig.tight_layout();fig.savefig(O/'figures/01_joint_umap.png',dpi=180);plt.close(fig)
x=a.layers['counts'];sy=a.var.symbol.astype(str).to_numpy();lib=ob.full_library_counts.to_numpy()
markers=['DCN','LUM','PDGFRA','FAP','THY1','ACTA2','TAGLN','MYL9','CXCL12','IL6','CCL2','HLA-DRA','CD74','COL1A1','COL3A1','POSTN','CTHRC1','MGP','ALPL','RUNX2','S100A4','PTPRC','PECAM1','EPCAM','RGS5','CSPG4'];markers=[g for g in markers if g in sy];means=[];dets=[]
for cl in clusters:
 m=ob.joint_cluster.astype(str).to_numpy()==cl;v=[];d=[]
 for g in markers:
  raw=np.asarray(x[m][:,np.flatnonzero(sy==g)].sum(1)).ravel();v.append(np.log1p(raw/lib[m]*1e4).mean());d.append((raw>0).mean())
 means.append(v);dets.append(d)
means=np.array(means);dets=np.array(dets);norm=(means-means.mean(0))/np.maximum(means.std(0),1e-9)
fig,ax=plt.subplots(figsize=(15,7));xx,yy=np.meshgrid(np.arange(len(markers)),np.arange(12));im=ax.scatter(xx.ravel(),yy.ravel(),s=dets.ravel()*180,c=norm.ravel(),cmap='RdBu_r',vmin=-2,vmax=2);ax.set_xticks(range(len(markers)),markers,rotation=60,ha='right');ax.set_yticks(range(12),['J'+c for c in clusters]);ax.invert_yaxis();ax.set_title('Raw RNA markers | dot area = detection; color = standardized mean logCP10k');fig.colorbar(im,ax=ax,label='Z score')
handles=[ax.scatter([],[],s=v*180,c='#777777',label=str(int(v*100))+'%') for v in [.1,.5,1]];ax.legend(handles=handles,title='Detection fraction',loc='upper center',bbox_to_anchor=(.5,-.20),ncol=3,frameon=False);fig.subplots_adjust(bottom=.27,right=.9);fig.savefig(O/'figures/05_marker_panel.png',dpi=180,bbox_inches='tight');plt.close(fig)
print('FIGURE_DISPLAY_FIXES_COMPLETE')