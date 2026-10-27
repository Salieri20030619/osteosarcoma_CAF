from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,json,hashlib,datetime,zipfile
from scipy import sparse
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import umap
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=R/'phase119_Bone_Lung_functional_subgroups_20261024';A=R/'phase120_patient_supported_CAF_review_20261026';O=R/'phase121_CAF_illustrated_report_20261026'
for x in ['code','figures','tables','resources','logs','delivery']:(O/x).mkdir(parents=True,exist_ok=True)
def read(n,base=I):return pd.read_csv(base/'tables'/n,sep='\t',low_memory=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def savefig(fig,name):
 fig.savefig(O/'figures'/(name+'.png'),dpi=250,bbox_inches='tight',facecolor='white');fig.savefig(O/'figures'/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
m=read('cohort_registry.tsv.gz');l=read('hierarchical_cell_labels.tsv.gz');nodes=read('hierarchical_nodes.tsv');rv=read('target_collagen_candidate_review.tsv');d=read('functional_audit.tsv',A);e=read('patient_effects.tsv.gz',A);sp=read('node_donor_support.tsv',A);dec=read('node_decisions_ZH.tsv',A)
assert m.cell_id.equals(l.cell_id);fit=~m.doublet;assert fit.sum()==6423
X=sparse.load_npz(I/'objects/log1p_cp10k.npz');genes=read('genes.tsv').symbol;gi={g:i for i,g in enumerate(genes)}
files=[I/'tables/cohort_registry.tsv.gz',I/'tables/hierarchical_cell_labels.tsv.gz',I/'objects/log1p_cp10k.npz',A/'tables/functional_audit.tsv',A/'PROTOCOL_V1.json'];before={str(p):sha(p) for p in files};dump('INPUT_HASHES.json',before)
core=['Bone_00','Bone_07','Bone_07_s1','Bone_07_s2','Bone_13','Lung_04','Lung_04_s2','Lung_07','Lung_11','Lung_12']
signal=sorted(set(core)|set(rv.loc[rv.collagen_evidence.isin(['nominal_collagen_candidate','consistent_direction_candidate']),'node_id']))
dump('PLAN_AND_SCOPE.json',dict(user_request='PDF with UMAP, marker dotplots and additional evidence plots; preserve signal groups, prioritize Lung_04',core_nodes=core,retained_signal_nodes=signal,no_reclustering=True,scientific_plotting='Authorized for this report; prior global figure backlog not resumed',umap='Visualization only on saved fitted PCA/Harmony representations; separate Bone/Lung, n_neighbors15,min_dist0.3,seed20261026',dotplot='Patient mean/detection first, average within source then equal weight sources represented in node; gene-wise z across displayed groups for colour; detection fraction for bubble size; predefined panel, not de novo DE markers',inference='Inherited Phase120 patient-level results; no changed thresholds or final-admission claims'))
coords=np.full((len(m),2),np.nan);pal={}
fig,axs=plt.subplots(1,2,figsize=(13,5.4))
for ax,site in zip(axs,['Bone','Lung']):
 z=np.load(I/'objects'/('main_'+site+'_root_representation.npz'));tr=z['transform_rows'];keep=~m.loc[tr,'doublet'].to_numpy();ix=tr[keep];rep=z['embedding'][keep]
 xy=umap.UMAP(n_neighbors=15,min_dist=.3,random_state=20261026,n_jobs=2).fit_transform(rep);coords[ix]=xy
 groups=sorted(l.loc[ix,'parent'].unique());cmap=plt.get_cmap('tab20');pal.update({g:matplotlib.colors.to_hex(cmap(i)) for i,g in enumerate(groups)})
 for g in groups:
  b=l.loc[ix,'parent'].eq(g).to_numpy();ax.scatter(xy[b,0],xy[b,1],s=4,c=pal[g],alpha=.8,rasterized=True);center=np.median(xy[b],axis=0);ax.text(*center,g.split('_')[1],fontsize=9,ha='center',va='center',bbox=dict(boxstyle='round,pad=.15',fc='white',ec='none',alpha=.75))
 ax.set(title=f'{site}: {len(ix):,} cells / {len(groups)} parents',xlabel='UMAP 1',ylabel='UMAP 2');ax.set_xticks([]);ax.set_yticks([])
fig.suptitle('Fixed parent labels on site-specific UMAP',fontsize=15);fig.tight_layout();savefig(fig,'01_parent_umap')
cc=m[['cell_id','site','source','donor','doublet']].copy();cc['parent']=l.parent;cc['leaf']=l.leaf;cc['UMAP1']=coords[:,0];cc['UMAP2']=coords[:,1];cc.to_csv(O/'tables/umap_coordinates.tsv.gz',sep='\t',index=False)
fig,axs=plt.subplots(2,2,figsize=(12,9))
for row,site in enumerate(['Bone','Lung']):
 ix=m.index[fit&m.site.eq(site)];xy=coords[ix]
 for src,color in zip(['SCPCP000017','SCPCP000023','NCI'],['#4e79a7','#f28e2b','#59a14f']):
  b=m.loc[ix,'source'].eq(src).to_numpy()
  if b.any():axs[row,0].scatter(xy[b,0],xy[b,1],s=5,c=color,label=src,alpha=.75,rasterized=True)
 axs[row,0].legend(fontsize=8,frameon=False);axs[row,0].set_title(site+' / source')
 ax=axs[row,1];ax.scatter(xy[:,0],xy[:,1],s=4,c='#d9dfe3',rasterized=True)
 targets=['Bone_07','Bone_07_s1','Bone_07_s2'] if site=='Bone' else ['Lung_04','Lung_04_s2']
 for g,color in zip(targets,['#1c7f89','#d34b42','#dfb34b']):
  b=l.loc[ix,'parent' if '_s' not in g else 'leaf'].eq(g).to_numpy();ax.scatter(xy[b,0],xy[b,1],s=10,c=color,label=f'{g} ({b.sum()})',rasterized=True)
 ax.legend(fontsize=8,frameon=False);ax.set_title(site+' / highlighted candidates')
 for ax in axs[row]:ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
fig.tight_layout();savefig(fig,'02_sources_candidates')
markers=['DCN','LUM','PDGFRA','COL1A1','COL1A2','COL3A1','COL5A1','POSTN','MMP11','FNDC1','COMP','COL11A1','SERPINH1','PCOLCE','LOX','LOXL2','ACTA2','TAGLN','MYH11','RGS5','CD36','C7','RUNX2','ALPL','IBSP','PTPRC','PECAM1'];assert all(g in gi for g in markers)
def profiles(groups):
 rows=[]
 for node in groups:
  col='leaf' if '_s' in node else 'parent';ix=m.index[fit&l[col].eq(node)]
  z=X[ix][:,[gi[g] for g in markers]].toarray();meta=m.loc[ix,['source','donor']].reset_index(drop=True)
  for src in meta.source.unique():
   av=[];dt=[]
   for _,idx in meta[meta.source.eq(src)].groupby('donor').groups.items():av.append(z[list(idx)].mean(0));dt.append((z[list(idx)]>0).mean(0))
   for j,g in enumerate(markers):rows.append(dict(node_id=node,source=src,gene=g,mean=np.mean(av,axis=0)[j],detection=np.mean(dt,axis=0)[j]))
 return pd.DataFrame(rows).groupby(['node_id','gene'])[['mean','detection']].mean().reset_index()
allgroups=nodes.node_id.tolist();pr=profiles(allgroups);pr.to_csv(O/'tables/marker_dotplot_values.tsv',sep='\t',index=False)
def dotplot(groups,title,name,figsize=(13,7)):
 a=pr[pr.node_id.isin(groups)].pivot(index='node_id',columns='gene',values='mean').reindex(index=groups,columns=markers);b=pr.pivot(index='node_id',columns='gene',values='detection').reindex(index=groups,columns=markers);std=a.std(axis=0,ddof=0).replace(0,1);z=(a-a.mean(0))/std
 fig,ax=plt.subplots(figsize=figsize);xx,yy=np.meshgrid(np.arange(len(markers)),np.arange(len(groups)));im=ax.scatter(xx.ravel(),yy.ravel(),s=b.to_numpy().ravel()*120,c=z.to_numpy().ravel(),cmap='RdBu_r',norm=Normalize(-2,2),edgecolors='none')
 ax.set_xticks(np.arange(len(markers)),markers,rotation=65,ha='right',fontsize=9);ax.set_yticks(np.arange(len(groups)),groups);ax.invert_yaxis();ax.set_ylim(len(groups)-.4,-.7);ax.set_title(title,pad=45);ax.grid(alpha=.14);ax.set_axisbelow(True)
 handles=[plt.scatter([],[],s=p*120,c='#777777') for p in [.1,.5,1]];ax.legend(handles,['10%','50%','100%'],title='RNA detection',ncol=3,loc='lower left',bbox_to_anchor=(0,1.01),frameon=False,fontsize=8);fig.colorbar(im,ax=ax,pad=.02,label='Gene-wise z-score of mean log1p(CP10k)',shrink=.65);fig.tight_layout();savefig(fig,name)
for site in ['Bone','Lung']:dotplot(nodes.loc[nodes.site.eq(site)&nodes.level.eq('parent'),'node_id'].tolist(),site+' parent marker panel','03_markers_'+site,(13,7))
dotplot(core,'Retained signal groups: marker panel','04_markers_candidates',(13,6.4))
for site in ['Bone','Lung']:
 q=sp[sp.site.eq(site)&sp.level.eq('parent')].copy();q['patient']=q.source+' / '+q.donor.str.split('|').str[-1].str.replace('NCI_PATIENT_','NCI_',regex=False)
 mat=q.pivot_table(index='node_id',columns='patient',values='cells',fill_value=0);frac=mat.div(mat.sum(1),axis=0)
 fig,ax=plt.subplots(figsize=(12,7));im=ax.imshow(frac,cmap='Blues',vmin=0,vmax=1,aspect='auto');ax.set_yticks(range(len(mat)),mat.index);ax.set_xticks(range(len(mat.columns)),[x.split(' / ')[-1] for x in mat.columns],rotation=65,ha='right',fontsize=8)
 for i in range(len(mat)):
  for j in range(len(mat.columns)):
   if mat.iloc[i,j]>0:ax.text(j,i,str(int(mat.iloc[i,j])),ha='center',va='center',fontsize=6.5,color='white' if frac.iloc[i,j]>.55 else '#333333')
 ax.set_title(site+' parent / patient contribution (numbers = cells)');fig.colorbar(im,ax=ax,label='Fraction of cells in parent',shrink=.7);fig.tight_layout();savefig(fig,'05_patient_'+site)
summary=dec.set_index('父群或子群').loc[core];fig,ax=plt.subplots(figsize=(11,6));y=np.arange(len(core));labels=['>=1 cell/side','>=5 cells/side','Identity + >=5 cells/side'];cols=['最少1细胞有效患者','最少5细胞有效患者','身份敏感性有效患者']
for j,(col,lab,color) in enumerate(zip(cols,labels,['#9cabb8','#278a94','#dc9e43'])):
 v=summary[col].to_numpy();ax.barh(y+(j-1)*.23,v,height=.21,label=lab,color=color)
 for yy,vv in zip(y+(j-1)*.23,v):ax.text(vv+.08,yy,str(vv),va='center',fontsize=8)
ax.axvline(6,color='#b23b3b',ls='--',lw=1.5,label='Minimum = 6');ax.set_yticks(y,core);ax.invert_yaxis();ax.set_xlim(0,12);ax.set_xlabel('Evaluable donor pairs in selected source');ax.legend(loc='lower right',fontsize=8);ax.set_title('Patient support, not cell abundance');fig.tight_layout();savefig(fig,'06_patient_gate')
ends=['program__mCAF_matrix','program__collagen_crosslinking','program__fibrillar_I_III_V','audit__FNDC1_paper'];titles=['Matrix program','Collagen crosslinking','Fibrillar I/III/V','FNDC1 paper panel']
fig,axs=plt.subplots(2,2,figsize=(12,8));fr=[]
for ax,ep,title in zip(axs.flat,ends,titles):
 for i,node in enumerate(core):
  source=summary.loc[node,'评估来源'];v=e[e.node_id.eq(node)&e.source.eq(source)&e['mode'].eq('main_min5')&e.endpoint.eq(ep)].effect.to_numpy();row=d[d.node_id.eq(node)&d.source.eq(source)&d['mode'].eq('main_min5')&d.endpoint.eq(ep)].iloc[0];fr.append(row.to_dict())
  ax.scatter(v,np.full(len(v),i),s=15,c='#89b6c3',alpha=.8);ax.scatter([row.mean_effect],[i],marker='D',s=32,c='#b8443d' if node=='Lung_04' else '#263d50',zorder=3)
 ax.axvline(0,c='gray',lw=.8);ax.set_yticks(range(len(core)),core,fontsize=8);ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('Patient mean difference (log1p CP10k)');ax.grid(axis='x',alpha=.15)
fig.suptitle('Within-source patient effects | dots: donors; diamond: mean',fontsize=13);fig.tight_layout();savefig(fig,'07_effects');pd.DataFrame(fr).to_csv(O/'tables/core_plot_statistics.tsv',sep='\t',index=False)
papers=['audit__FNDC1_paper','audit__MYH11_paper','audit__CD36_paper','audit__C7_paper'];mat=[]
for node in core:
 src=summary.loc[node,'评估来源'];q=d[d.node_id.eq(node)&d.source.eq(src)&d['mode'].eq('main_min5')].set_index('endpoint');mat.append([q.loc[ep,'mean_effect'] for ep in papers])
mat=np.array(mat);fig,ax=plt.subplots(figsize=(8,6.5));mx=max(abs(mat[np.isfinite(mat)]).max(),.1);im=ax.imshow(mat,cmap='RdBu_r',vmin=-mx,vmax=mx,aspect='auto');ax.set_yticks(range(len(core)),core);ax.set_xticks(range(4),['FNDC1 / matrix','MYH11 / vascular','CD36 / lipid','C7 / inflammatory'],rotation=20,ha='right')
for i in range(len(core)):
 for j in range(4):ax.text(j,i,'NA' if not np.isfinite(mat[i,j]) else f'{mat[i,j]:+.2f}',ha='center',va='center',fontsize=9,color='white' if abs(mat[i,j])>.65*mx else '#222222')
ax.set_title('Paper-derived programs: min5 patient contrasts');fig.colorbar(im,ax=ax,label='Mean log-expression difference',shrink=.8);fig.tight_layout();savefig(fig,'08_literature')
fig,axs=plt.subplots(1,2,figsize=(12,4.5));lung=sp[sp.node_id.eq('Lung_04')&sp.source.eq('NCI')].sort_values('cells',ascending=False);x=np.arange(len(lung));axs[0].bar(x,lung.cells,color='#a6bdca',label='Main cells');axs[0].bar(x,lung.identity_cells,color='#278a94',label='Identity sensitivity');axs[0].set_xticks(x,lung.donor.str.replace('NCI_PATIENT_','',regex=False),rotation=55);axs[0].set_ylabel('Target cells');axs[0].set_title('Lung_04: donor coverage');axs[0].legend(fontsize=8)
q=pd.read_csv(A/'tables/Lung_04_identity_attrition.tsv',sep='\t');q=q[q.source.eq('NCI')&q.donor.isin(['NCI_PATIENT_6546','NCI_PATIENT_7978'])&~q.target];xx=np.arange(len(q));axs[1].bar(xx-.18,q.cells,width=.36,color='#a6bdca',label='Main control');axs[1].bar(xx+.18,q.identity,width=.36,color='#278a94',label='Identity control');axs[1].axhline(5,color='#b23b3b',ls='--',label='Minimum = 5');axs[1].set_xticks(xx,q.donor.str.replace('NCI_PATIENT_','',regex=False));axs[1].set_title('Two lost pairs: control-cell attrition');axs[1].set_ylabel('Control cells');axs[1].legend(fontsize=8);fig.tight_layout();savefig(fig,'09_lung04_attrition')
dec[dec['父群或子群'].isin(signal)].to_csv(O/'tables/retained_signal_groups.tsv',sep='\t',index=False)
for src in [A/'resources/LITERATURE_CROSSCHECK.json',A/'tables/focus_functional_audit.tsv',A/'tables/node_decisions_ZH.tsv',I/'tables/target_collagen_candidate_review.tsv']:
 import shutil
 shutil.copyfile(src,O/'resources'/src.name)
for p,h in before.items():assert sha(Path(p))==h
figs=sorted((O/'figures').glob('*.png'));assert len(figs)==11
dump('FIGURES_COMPLETE.json',dict(status='PASS',figures=len(figs),cells=6633,plotted_main_cells=6423,umap_cells=int(np.isfinite(coords[:,0]).sum()),labels_unchanged=True,input_hashes_unchanged=True,retained_signal_nodes=signal,core_nodes=core,figure_hashes={p.name:sha(p) for p in figs}))
with zipfile.ZipFile(O/'delivery/report_inputs.zip','w',zipfile.ZIP_DEFLATED) as z:
 for directory in ['figures','tables','resources']:
  for p in (O/directory).glob('*'):
   if p.is_file() and p.suffix!='.pdf':z.write(p,str(p.relative_to(O)))
 z.write(O/'FIGURES_COMPLETE.json','FIGURES_COMPLETE.json')
print(json.dumps(dict(status='PASS',figures=len(figs),retained=len(signal),archive_sha256=sha(O/'delivery/report_inputs.zip'))),flush=True)
