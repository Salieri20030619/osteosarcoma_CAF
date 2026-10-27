from pathlib import Path
import json,hashlib,math
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from PIL import Image
import plotly.graph_objects as go
P=Path('C:/Users/31937/CAF_paper_figures_20261021/allcell_3d');P.mkdir(exist_ok=True)
O=Path('C:/Users/31937/Downloads/CAF_allcell_3D_Fig1b_1d_20261021');O.mkdir(exist_ok=True)
for s in ['previews','resources']:(P/s).mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7,'savefig.facecolor':'white'})
palette={'T cell':'#E34A33','Macrophage':'#E3A13B','Monocyte / DC unresolved':'#BC7D29','pDC-like':'#A35F18','Mixed B / plasma':'#B887CB','Plasma-like':'#8B4DAA','Mast':'#C05588','Endothelial':'#2A9D8F','Lymphatic endothelium':'#7CC9BA','Fibroblast-like':'#79A737','Mural / pericyte':'#397E71','Osteogenic-like':'#6C65A7','Chondrocyte-like':'#8F9ACB','Mesenchymal unresolved':'#B7B774','Myogenic':'#A1785C','Osteoclast-like':'#A74555','Epithelial':'#3689C7','Renal epithelial-like':'#65B8E1','Erythroid':'#85454A','Mixed endothelial / epithelial':'#6D9B9B','Unresolved lineage':'#B7BEC5'}
order=list(palette);manifest=[];interactive=[]
sites=['Bone','Lung','Unknown_or_mixed'];pretty={'All':'All cells','Bone':'Bone','Lung':'Lung','Unknown_or_mixed':'Unknown / mixed site'}
def save(fig,name,kind,source,scope):
 exports=[]
 for ext in ['png','tiff']:
  path=O/(name+'.'+ext);fig.savefig(path,dpi=600,**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
  with Image.open(path) as im:
   im.load();dpi=list(map(float,im.info['dpi']));assert max(abs(v-600) for v in dpi)<.05
   assert im.size==tuple(round(v*600) for v in fig.get_size_inches())
   exports.append(dict(file=str(path),pixels=im.size,dpi=dpi,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 fig.savefig(P/'previews'/(name+'.png'),dpi=110);plt.close(fig)
 manifest.append(dict(name=name,kind=kind,source=source,scope=scope,exports=exports));print('SAVED',name,flush=True)
def handles(labels):return [Line2D([0],[0],marker='o',ls='',color=palette[t],markersize=6,label=t) for t in labels]
means=pd.read_csv(P/'site_patient_balanced_composition.tsv',sep='\t');comp=pd.read_csv(P/'patient_site_composition.tsv',sep='\t');cohort=pd.read_csv(P/'cohort_counts.tsv',sep='\t')
for source in ['SCPCP000017','SCPCP000023']:
 d=pd.read_csv(P/(source+'_plot_cells.tsv.gz'),sep='\t');assert d.cell_id.is_unique and np.isfinite(d[['UMAP1','UMAP2','UMAP3']]).all().all()
 labels=[t for t in order if t in set(d.cell_type)];assert set(d.cell_type)<=set(palette)
 xyz=d[['UMAP1','UMAP2','UMAP3']].to_numpy();lo=xyz.min(0);hi=xyz.max(0);pad=(hi-lo)*.05
 for site in ['All','Bone','Lung']:
  q=d if site=='All' else d[d.site.eq(site)]
  # All cells shown, deterministic drawing order; no display subsampling.
  ix=np.random.default_rng(20261021).permutation(len(q));q=q.iloc[ix]
  fig=plt.figure(figsize=(11.8,8.4),dpi=120);ax=fig.add_axes([.015,.13,.70,.73],projection='3d')
  ax.scatter(q.UMAP1,q.UMAP2,q.UMAP3,c=q.cell_type.map(palette),s=.65,alpha=.8,depthshade=False,linewidths=0,rasterized=True)
  ax.view_init(elev=20,azim=-55);ax.set_box_aspect((1,1,.85))
  for j,(setter,label) in enumerate([(ax.set_xlim,'UMAP 1'),(ax.set_ylim,'UMAP 2'),(ax.set_zlim,'UMAP 3')]):setter(lo[j]-pad[j],hi[j]+pad[j])
  ax.set_xlabel('UMAP 1',labelpad=6);ax.set_ylabel('UMAP 2',labelpad=6);ax.set_zlabel('UMAP 3',labelpad=5)
  for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
   axis.set_major_locator(MaxNLocator(4));axis.pane.fill=False;axis._axinfo['grid'].update(color='#E3E7EB',linewidth=.5)
  ax.tick_params(labelsize=8,pad=1)
  fig.text(.06,.94,f'{source} | {pretty[site]} | 3D UMAP',fontsize=17,fontweight='bold')
  fig.text(.06,.895,f'{len(q):,} cells | {q.donor_record.nunique()} patient records | fixed source coordinates and view',fontsize=11,color='#53616B')
  fig.legend(handles=handles(labels),loc='upper left',bbox_to_anchor=(.735,.84),frameon=False,fontsize=9,labelspacing=.72,handletextpad=.6)
  fig.text(.06,.07,'Frozen reviewed expression lineages; unresolved and mixed populations retained.',fontsize=10,color='#53616B')
  fig.text(.06,.035,'All upstream called cells shown. UMAP axes represent expression structure, not tissue location.',fontsize=10,color='#53616B')
  save(fig,f'Fig1b_{source}_{site}_3D_UMAP','UMAP3D',source,site)
 # Original 2D coordinates are delivered as a matched-context comparison, not projected 3D.
 fig=plt.figure(figsize=(11.8,8.4),dpi=120);ax=fig.add_axes([.08,.18,.61,.63]);q=d.iloc[np.random.default_rng(20261021).permutation(len(d))]
 ax.scatter(q.legacy_UMAP1,q.legacy_UMAP2,c=q.cell_type.map(palette),s=.65,alpha=.8,linewidths=0,rasterized=True)
 ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False);ax.set_aspect('equal',adjustable='datalim')
 fig.text(.06,.94,f'{source} | All cells | original 2D UMAP',fontsize=17,fontweight='bold')
 fig.text(.06,.895,f'{len(d):,} cells | same frozen graph and annotation as the 3D view',fontsize=11,color='#53616B')
 fig.legend(handles=handles(labels),loc='upper left',bbox_to_anchor=(.735,.84),frameon=False,fontsize=9,labelspacing=.72)
 fig.text(.06,.065,'Phase53 two-dimensional embedding retained; this is not a 2D projection of the new 3D embedding.',fontsize=9,color='#53616B')
 fig.text(.06,.035,'Frozen reviewed expression lineages. Separate source embeddings cannot be compared by absolute coordinates.',fontsize=9,color='#53616B')
 save(fig,f'Fig1b_{source}_All_original_2D_UMAP','UMAP2D_reference',source,'All')
 # Interactive self-contained HTML: each cell appears once; site filtering changes trace visibility.
 traces=[];tracekeys=[]
 for site in sites:
  for celltype in labels:
   z=d[d.site.eq(site)&d.cell_type.eq(celltype)]
   if z.empty:continue
   traces.append(go.Scatter3d(x=z.UMAP1.round(5),y=z.UMAP2.round(5),z=z.UMAP3.round(5),mode='markers',name=celltype,legendgroup=celltype,marker=dict(size=1.5,color=palette[celltype],opacity=.8),text=z.donor_record.astype(str)+' | '+pretty[site],hovertemplate=celltype+'<br>%{text}<extra></extra>'))
   tracekeys.append((site,celltype))
 buttons=[]
 for site in ['All']+sites:
  visible=[site=='All' or a==site for a,b in tracekeys];seen=set();show=[]
  for vis,(ss,tt) in zip(visible,tracekeys):show.append(vis and tt not in seen);seen.add(tt) if vis else None
  if site=='All':
   for trace,sl in zip(traces,show):trace.showlegend=sl
  buttons.append(dict(label=pretty[site],method='update',args=[dict(visible=visible,showlegend=show),dict(title=dict(text=source+' | '+pretty[site]+' | all-cell 3D UMAP',x=.02,y=.985,xanchor='left',yanchor='top'))]))
 pl=go.Figure(traces);pl.update_layout(title=dict(text=source+' | All cells | all-cell 3D UMAP',x=.02,y=.985,xanchor='left',yanchor='top'),scene=dict(xaxis=dict(title='UMAP 1',range=[lo[0]-pad[0],hi[0]+pad[0]]),yaxis=dict(title='UMAP 2',range=[lo[1]-pad[1],hi[1]+pad[1]]),zaxis=dict(title='UMAP 3',range=[lo[2]-pad[2],hi[2]+pad[2]]),aspectmode='cube',camera=dict(eye=dict(x=1.6,y=-1.9,z=.9))),updatemenus=[dict(buttons=buttons,x=0,y=1.06,direction='right',type='buttons')],margin=dict(l=0,r=0,t=160,b=50),legend=dict(font=dict(size=11)),annotations=[dict(text='All upstream called cells; frozen reviewed lineages. Rotate / zoom; filter sites above. Coordinates are source-specific.',xref='paper',yref='paper',x=0,y=-.05,showarrow=False)],uirevision='fixed_scene')
 html=O/f'Fig1b_{source}_interactive_3D.html';pl.write_html(html,include_plotlyjs=True,full_html=True);assert sum(len(t.x) for t in traces)==len(d)
 interactive.append(dict(file=str(html),cells=len(d),traces=len(traces),sha256=hashlib.sha256(html.read_bytes()).hexdigest(),all_cells_once=True));print('INTERACTIVE_SAVED',source,flush=True)
 # Patient-balanced site composition, three bars, preserving unknown-site context.
 fig=plt.figure(figsize=(11.8,8.4),dpi=120);ax=fig.add_axes([.10,.20,.55,.58]);bottom=np.zeros(3)
 for label in labels:
  rows=means[means.dataset.eq(source)&means.cell_type.eq(label)].set_index('site');vals=rows.loc[sites,'mean_patient_percent'].to_numpy()
  ax.bar(range(3),vals,bottom=bottom,width=.57,color=palette[label],edgecolor='white',lw=.4);bottom+=vals
 assert np.allclose(bottom,100)
 for i,site in enumerate(sites):
  row=cohort[cohort.dataset.eq(source)&cohort.site.eq(site)].iloc[0]
  ax.text(i,102,f'n={int(row.patients)}\n{int(row.cells):,} cells',ha='center',fontsize=10,linespacing=1.4)
 ax.set_ylim(0,100);ax.set_yticks([0,25,50,75,100]);ax.set_xticks(range(3),['Bone','Lung','Unknown / mixed']);ax.set_ylabel('Mean patient cell-type proportion (%)',fontsize=11);ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='x',length=0,pad=10)
 fig.text(.06,.94,f'{source} | Cell-type composition',fontsize=17,fontweight='bold')
 fig.text(.06,.895,'All-cell denominator; equal patient weighting within each anatomical site',fontsize=11,color='#53616B')
 fig.legend(handles=[Patch(facecolor=palette[t],label=t) for t in labels],loc='upper left',bbox_to_anchor=(.70,.83),frameon=False,fontsize=9,labelspacing=.72)
 fig.text(.06,.085,'Each bar averages patient-specific proportions; it is not the proportion after pooling all cells.',fontsize=10,color='#53616B')
 fig.text(.06,.045,'Shared patients may contribute to both sites. Descriptive composition; no independent-group significance claim.',fontsize=9,color='#53616B')
 save(fig,f'Fig1d_{source}_patient_balanced_site_composition','site_composition',source,'All')
 # Per-patient/site bars expose heterogeneity and all denominators.
 q=comp[comp.dataset.eq(source)].copy();meta=q[['site','donor','denominator']].drop_duplicates();meta['site_order']=meta.site.map({s:i for i,s in enumerate(sites)});meta=meta.sort_values(['site_order','donor'])
 fig=plt.figure(figsize=(13.4,8.8),dpi=120);ax=fig.add_axes([.07,.32,.64,.44]);bottom=np.zeros(len(meta))
 for label in labels:
  lookup=q[q.cell_type.eq(label)].set_index(['site','donor']).fraction_percent
  vals=np.array([lookup.loc[(r.site,r.donor)] for r in meta.itertuples()]);ax.bar(range(len(meta)),vals,bottom=bottom,width=.82,color=palette[label],edgecolor='white',lw=.15);bottom+=vals
 assert np.allclose(bottom,100)
 for i,row in enumerate(meta.itertuples()):ax.text(i,101,str(row.denominator),ha='center',fontsize=6.5,rotation=90)
 for site,gg in meta.groupby('site_order',sort=True):
  ix=np.flatnonzero(meta.site_order.eq(site));ax.text(float(ix.mean()),1.22,pretty[sites[site]],transform=ax.get_xaxis_transform(),ha='center',fontsize=10)
  if ix[0]>0:ax.axvline(ix[0]-.5,color='#53616B',lw=.8)
 ax.set_ylim(0,100);ax.set_xlim(-.6,len(meta)-.4);ax.set_yticks([0,25,50,75,100]);ax.set_ylabel('Cell-type proportion per patient/site (%)',fontsize=10)
 ax.set_xticks(range(len(meta)),[x.split('|')[-1] for x in meta.donor],rotation=65,ha='right',fontsize=7);ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='x',length=0)
 fig.text(.06,.94,f'{source} | Patient-level cell-type composition',fontsize=17,fontweight='bold')
 fig.text(.06,.895,'One bar per patient/site; all-cell counts shown above bars',fontsize=11,color='#53616B')
 fig.legend(handles=[Patch(facecolor=palette[t],label=t) for t in labels],loc='upper left',bbox_to_anchor=(.75,.83),frameon=False,fontsize=9,labelspacing=.72)
 fig.text(.06,.09,'All upstream called cells; frozen reviewed lineages including mixed and unresolved populations.',fontsize=10,color='#53616B')
 fig.text(.06,.05,'Repeated patient IDs across sites are retained and documented; cell counts are not biological sample sizes.',fontsize=9,color='#53616B')
 save(fig,f'Fig1d_{source}_per_patient_composition','patient_composition',source,'All')
out=dict(status='export_pass_visual_pending',static_figures=len(manifest),interactive_files=interactive,figures=manifest,all_static_600dpi=True,all_3D_cells_plotted=True,no_display_downsampling=True,fixed_camera_and_axis_per_source=True,clinical_group_coordinates_not_refitted=True)
(P/'resources/FIGURE_MANIFEST.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print('PLOTTING_COMPLETE',len(manifest),flush=True)
