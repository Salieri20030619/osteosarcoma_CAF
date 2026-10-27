from pathlib import Path
import json,hashlib,shutil
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase121_CAF_illustrated_report_20261026';A=R/'phase120_patient_supported_CAF_review_20261026'
q=pd.read_csv(O/'tables/umap_coordinates.tsv.gz',sep='\t');q=q[~q.doublet];assert len(q)==6423
backup=O/'resources/before_legend_revision';backup.mkdir(exist_ok=True)
for n in ['01_parent_umap','02_sources_candidates','07_effects']:
 for ext in ['png','pdf']:
  p=O/'figures'/(n+'.'+ext)
  if not (backup/p.name).exists():shutil.copyfile(p,backup/p.name)
if not (backup/'CAF_illustrated_report.pdf').exists():shutil.copyfile(O/'delivery/CAF_illustrated_report.pdf',backup/'CAF_illustrated_report.pdf')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(O/'figures'/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
 plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(14,6.3))
for ax,site in zip(axs,['Bone','Lung']):
 z=q[q.site.eq(site)];groups=sorted(z.parent.unique());cmap=plt.get_cmap('tab20');handles=[]
 for i,g in enumerate(groups):
  b=z[z.parent.eq(g)];color=cmap(i);ax.scatter(b.UMAP1,b.UMAP2,s=4,c=[color],alpha=.8,rasterized=True);center=b[['UMAP1','UMAP2']].median();ax.text(*center,g.split('_')[1],fontsize=9,ha='center',va='center',bbox=dict(boxstyle='round,pad=.15',fc='white',ec='none',alpha=.75));handles.append(Line2D([],[],marker='o',ls='',color=color,label=g,markersize=6))
 ax.set(title=f'{site}: {len(z):,} cells / {len(groups)} parents',xlabel='UMAP 1',ylabel='UMAP 2');ax.set_xticks([]);ax.set_yticks([]);ax.legend(handles=handles,title='Parent cluster',loc='upper center',bbox_to_anchor=(.5,-.10),ncol=4,frameon=False,fontsize=8,columnspacing=.8,handletextpad=.4)
fig.suptitle('Fixed parent labels on site-specific UMAP',fontsize=15);fig.subplots_adjust(bottom=.27,wspace=.18,top=.88);save(fig,'01_parent_umap')
fig,axs=plt.subplots(2,2,figsize=(13,9))
for row,site in enumerate(['Bone','Lung']):
 z=q[q.site.eq(site)];ax=axs[row,0];handles=[]
 for src,color in zip(['SCPCP000017','SCPCP000023','NCI'],['#4e79a7','#f28e2b','#59a14f']):
  b=z[z.source.eq(src)]
  if len(b):ax.scatter(b.UMAP1,b.UMAP2,s=5,c=color,alpha=.75,rasterized=True);handles.append(Line2D([],[],marker='o',ls='',color=color,label=src,markersize=6))
 ax.legend(handles=handles,title='Source',fontsize=8,frameon=False,loc='upper right');ax.set_title(site+' / source')
 ax=axs[row,1];parent='Bone_07' if site=='Bone' else 'Lung_04';children=['Bone_07_s1','Bone_07_s2'] if site=='Bone' else ['Lung_04_s2'];colors=['#d34b42','#dfb34b'];handles=[]
 for mask,color,label in [(~z.parent.eq(parent),'#d9dfe3','Other parents'),(z.parent.eq(parent)&~z.leaf.isin(children),'#1c7f89',parent+' / remaining cells')]+[(z.leaf.eq(ch),co,ch) for ch,co in zip(children,colors)]:
  b=z[mask];ax.scatter(b.UMAP1,b.UMAP2,s=4 if label=='Other parents' else 10,c=color,rasterized=True);handles.append(Line2D([],[],marker='o',ls='',color=color,label=label+' (n='+str(len(b))+')',markersize=6))
 ax.legend(handles=handles,title='Displayed cell groups',fontsize=7.5,frameon=False,loc='upper right');ax.set_title(site+' / highlighted candidates')
 for ax in axs[row]:ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
fig.tight_layout();save(fig,'02_sources_candidates')
core=json.loads((O/'PLAN_AND_SCOPE.json').read_text())['core_nodes'];dec=pd.read_csv(A/'tables/node_decisions_ZH.tsv',sep='\t').set_index('父群或子群');e=pd.read_csv(A/'tables/patient_effects.tsv.gz',sep='\t');stats=pd.read_csv(O/'tables/core_plot_statistics.tsv',sep='\t')
ends=['program__mCAF_matrix','program__collagen_crosslinking','program__fibrillar_I_III_V','audit__FNDC1_paper'];fig,axs=plt.subplots(2,2,figsize=(13,9))
for ax,ep,title in zip(axs.flat,ends,['Matrix program','Collagen crosslinking','Fibrillar I/III/V','FNDC1 paper panel']):
 for i,node in enumerate(core):
  src=dec.loc[node,'评估来源'];v=e[e.node_id.eq(node)&e.source.eq(src)&e['mode'].eq('main_min5')&e.endpoint.eq(ep)].effect.to_numpy();mu=stats[stats.node_id.eq(node)&stats.endpoint.eq(ep)].mean_effect.iloc[0];ax.scatter(v,np.full(len(v),i),s=18,c='#89b6c3',alpha=.8);ax.scatter([mu],[i],marker='D',s=35,c='#b8443d' if node=='Lung_04' else '#263d50',zorder=3)
 ax.axvline(0,c='gray',lw=.8);ax.set_yticks(range(len(core)),core,fontsize=9);ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('Patient mean difference (log1p CP10k)');ax.grid(axis='x',alpha=.15)
handles=[Line2D([],[],marker='o',ls='',c='#89b6c3',label='Individual donor effect'),Line2D([],[],marker='D',ls='',c='#263d50',label='Other group mean'),Line2D([],[],marker='D',ls='',c='#b8443d',label='Lung_04 mean'),Line2D([],[],c='gray',label='Zero difference')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.965),ncol=4,frameon=False,fontsize=10);fig.suptitle('Within-source patient effects',fontsize=14);fig.tight_layout(rect=[0,0,1,.91]);save(fig,'07_effects')
receipt=dict(status='PASS',changed_figures=['01_parent_umap','02_sources_candidates','07_effects'],unchanged='UMAP coordinates, cells, cluster labels, expression and statistical values',changes=['Complete colour-to-parent legends for Bone and Lung UMAP','Grey background and nonoverlapping highlighted groups explicitly labelled','Donor dots, other means, Lung_04 mean and zero line explicitly labelled'],other_figures='Marker plots have size legends and colour bars; patient heatmaps and literature heatmap have labelled colour bars; support/attrition bar plots have categorical legends',hashes={n:hashlib.sha256((O/'figures'/(n+'.png')).read_bytes()).hexdigest() for n in ['01_parent_umap','02_sources_candidates','07_effects']})
(O/'resources/LEGEND_REVISION.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt))
