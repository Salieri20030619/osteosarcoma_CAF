from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='2'
import json,hashlib,datetime,zipfile,traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=R/'phase114_CAF_UMAP_FeaturePlots_20261023'
genes=['PNKD','TOR1AIP1','GRID2','TAF10','BEST1','RGS2-AS1']
sources=['SCPCP000017','SCPCP000023','NCI']
def dump(name,x): (O/name).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as r:
  for b in iter(lambda:r.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 f=R/'phase111_reference_figure_adaptation_20261021/tables/Fig1i_cell_values.tsv.gz'
 inputs=[f,O/'tables/author_umap_coordinates.tsv',R/'patient_expansion_download_investigation_20261013_v1/downloads/OsteoCAR/patient_mets.qs2']
 c=pd.read_csv(f,sep='\t');assert len(c)==6633 and c.cell_id.is_unique
 assert set(c.source)==set(sources)
 parts=[];coverage=[]
 for src in sources:
  q=c[c.source.eq(src)].copy()
  if src=='NCI':
   a=pd.read_csv(inputs[1],sep='\t').rename(columns={'umap_1':'UMAP1','umap_2':'UMAP2'})
   origin='Author patient_mets full-object UMAP; NCI author CAF subset'
  else:
   p=R/'phase53_de_novo_reannotation_20261002/objects'/src/'embeddings.npz';inputs.append(p)
   with np.load(p) as z:a=pd.DataFrame({'cell_id':z['cell_id'],'UMAP1':z['UMAP'][:,0],'UMAP2':z['UMAP'][:,1]})
   origin='Frozen Phase53 source all-cell UMAP; current candidate subset'
  assert a.cell_id.is_unique
  q=q.merge(a,on='cell_id',how='left',validate='one_to_one')
  assert np.isfinite(q[['UMAP1','UMAP2']].to_numpy()).all()
  for g in genes:
   assert (q[g+'__raw']>=0).all()
   assert np.allclose(q[g+'__log'],np.log1p(q[g+'__raw']/q.full_library_counts*1e4),rtol=0,atol=1e-12)
  coverage.append(dict(source=src,cells=len(q),source_patient_records=q.donor.nunique(),matched=len(q),coordinate_origin=origin,classes=q['class'].value_counts().to_dict()))
  parts.append(q)
 data=pd.concat(parts,ignore_index=True);assert set(data.cell_id)==set(c.cell_id)
 data.to_csv(O/'tables/plot_cells.tsv.gz',sep='\t',index=False,encoding='utf-8')
 hashes=[dict(path=str(p),sha256=sha(p)) for p in inputs]
 dump('resources/INPUT_HASHES.json',hashes)
 protocol=dict(cells=6633,genes=genes,coverage=coverage,scope='Phase103 expanded original candidates plus NCI author CAF; frozen membership and classes',coordinate_policy='Subset existing embeddings without refitting. Different source panels have independent coordinate systems. Author NCI coordinates originate from the broader 221271-cell patient_mets object.',expression='log1p(raw RNA counts / full cell library counts * 10000)',color_scale='Per gene, common zero-to-maximum scale across all three displayed sources, with no clipping',draw_order='Zeros first, positive cells ascending expression',inference='Descriptive cell maps; no new statistical test or protein/flow validation',identity='F1-F5 are inherited candidate classes, not purified CAF identity. NCI uses author CAF label, not F1-F5 transfer.',embedding_gene_exclusion='No new feature selection. Six displayed genes were not newly excluded from inherited embeddings; these maps are not independent evidence of separation.',patient_counts='Source-local patient records; not proof of cross-study independence',NCI_feature_coverage='Inherited author merged-object feature coverage; per-library premerge feature masks unavailable',no_reclustering=True,original_results_modified=False)
 dump('resources/PROTOCOL.json',protocol)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':14,'axes.labelsize':10,'savefig.facecolor':'white'})
 palette={'F1':'#2878B5','F2':'#E5A326','F3':'#409A71','F4':'#B55A88','F5':'#8265AE','Author CAF':'#258B8B'}
 cmap=LinearSegmentedColormap.from_list('rna',['#dddddd','#fee0b6','#ef8a62','#b2182b','#67001f'])
 manifest=[]
 for src in sources:
  q=data[data.source.eq(src)].copy();x=q.UMAP1.to_numpy();y=q.UMAP2.to_numpy()
  dx=max(np.ptp(x),1);dy=max(np.ptp(y),1)
  lims=((x.min()-.07*dx,x.max()+.07*dx),(y.min()-.07*dy,y.max()+.07*dy))
  for gene in [None]+genes:
   fig,ax=plt.subplots(figsize=(7.4,6.2));fig.subplots_adjust(left=.12,right=.78,bottom=.17,top=.81)
   title='CAF candidate classes' if gene is None else gene+' RNA expression'
   fig.text(.12,.94,title,fontsize=16,weight='bold')
   fig.text(.12,.895,src+' | '+format(len(q),',')+' cells | '+str(q.donor.nunique())+' source-local patient records',fontsize=10)
   fig.text(.12,.855,'Author CAF subset' if src=='NCI' else 'Expanded candidate pool | frozen F1-F5 labels',fontsize=9,color='#555555')
   if gene is None:
    shuffled=q.sample(frac=1,random_state=20261023)
    ax.scatter(shuffled.UMAP1,shuffled.UMAP2,c=shuffled['class'].map(palette),s=10,linewidths=0,alpha=.85)
    for cl in sorted(q['class'].unique()):ax.scatter([],[],c=palette[cl],s=26,label=cl+('\n' if cl=='Author CAF' else ' ')+'(n='+str((q['class']==cl).sum())+')')
    ax.legend(loc='upper left',bbox_to_anchor=(1.02,1.02),frameon=False,fontsize=8,handletextpad=.3)
    name='Fig1g_'+src+'_CAF_candidate_UMAP';extra={}
   else:
    vals=q[gene+'__log'].to_numpy();order=np.argsort(vals,kind='stable');vmax=float(data[gene+'__log'].max());assert vmax>0
    m=ax.scatter(x[order],y[order],c=vals[order],cmap=cmap,vmin=0,vmax=vmax,s=11,linewidths=0)
    cax=fig.add_axes([.825,.35,.025,.30]);cb=fig.colorbar(m,cax=cax);cb.set_label('log1p(CP10k)',fontsize=9);cb.ax.tick_params(labelsize=8)
    fig.text(.81,.73,'RNA detected\n'+str(int((vals>0).sum()))+' / '+str(len(q))+' cells',fontsize=8,linespacing=1.7)
    fig.text(.81,.22,'Same gene:\nshared color scale\nacross sources',fontsize=8,linespacing=1.5,color='#555555')
    name='Fig1j_'+src+'_'+gene+'_FeaturePlot';extra=dict(gene=gene,detected=int((vals>0).sum()),vmin=0,vmax=vmax)
   ax.set(xlim=lims[0],ylim=lims[1],xlabel='UMAP 1',ylabel='UMAP 2');ax.set_aspect('equal',adjustable='box')
   ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
   fig.text(.12,.065,'Existing all-cell embedding; current candidate cells shown.\nSource coordinates are independent. Descriptive RNA view.',fontsize=8,color='#555555',linespacing=1.5)
   png=O/'figures'/(name+'.png');fig.savefig(png,dpi=600)
   fig.savefig(O/'previews'/(name+'.png'),dpi=140);plt.close(fig)
   tif=png.with_suffix('.tiff')
   with Image.open(png) as im:
    im.convert('RGB').save(tif,compression='tiff_lzw',dpi=(600,600))
    dims=im.size;assert dims==(4440,3720) and abs(im.info['dpi'][0]-600)<.1
    with Image.open(tif) as t:assert np.array_equal(np.asarray(im.convert('RGB')),np.asarray(t)) and t.size==dims
   manifest.append(dict(name=name,source=src,cells=len(q),pixels=list(dims),dpi=600,PNG_TIFF_pixels_identical=True,files=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in [png,tif]],**extra))
   print('DONE',name,flush=True)
 assert len(manifest)==21
 for h in hashes:assert sha(h['path'])==h['sha256']
 dump('resources/FIGURE_MANIFEST.json',manifest)
 dump('resources/DATA_QA.json',dict(status='PASS',matched_cells=6633,expected_cells=6633,no_duplicate_cell_ids=True,expression_normalization_reproduced=True,source_inputs_unchanged=True,all_21_PNG_TIFF_pairs_pixel_identical=True,coverage=coverage,visual_QA='pending'))
 dump('STATUS.json',dict(status='numeric_QA_complete_visual_QA_pending',figures=21,updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 # Legible overview sheets for visual review; originals remain separate files.
 for src in sources:
  files=sorted((O/'previews').glob('*'+src+'*.png'))
  sheet=Image.new('RGB',(1554,2604),'white')
  for i,f in enumerate(files):
   with Image.open(f) as im:
    im=im.convert('RGB');im.thumbnail((777,651));sheet.paste(im,((i%2)*777,(i//2)*651))
  sheet.save(O/'previews'/('CONTACT_'+src+'.jpg'),quality=95)
if __name__=='__main__':
 try:main()
 except Exception:
  dump('FAILURE.json',dict(traceback=traceback.format_exc()));raise
