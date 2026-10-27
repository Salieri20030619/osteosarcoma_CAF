from pathlib import Path
import pandas as pd,numpy as np,json,shutil,hashlib
from scipy import sparse
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=R/'phase119_Bone_Lung_functional_subgroups_20261024';O=R/'phase121_CAF_illustrated_report_20261026'
m=pd.read_csv(I/'tables/cohort_registry.tsv.gz',sep='\t');l=pd.read_csv(I/'tables/hierarchical_cell_labels.tsv.gz',sep='\t');nodes=pd.read_csv(I/'tables/hierarchical_nodes.tsv',sep='\t');old=pd.read_csv(O/'tables/marker_dotplot_values.tsv',sep='\t');genes=pd.read_csv(I/'tables/genes.tsv',sep='\t').symbol;gi={g:i for i,g in enumerate(genes)};markers=['DCN','LUM','PDGFRA','COL1A1','COL1A2','COL3A1','COL5A1','POSTN','MMP11','FNDC1','COMP','COL11A1','SERPINH1','PCOLCE','LOX','LOXL2','ACTA2','TAGLN','MYH11','RGS5','CD36','C7','RUNX2','ALPL','IBSP','PTPRC','PECAM1'];X=sparse.load_npz(I/'objects/log1p_cp10k.npz')[:,[gi[g] for g in markers]].toarray();records=[]
backup=O/'resources/source_balanced_dotplots_v1';backup.mkdir(exist_ok=True)
for f in [O/'tables/marker_dotplot_values.tsv',O/'figures/03_markers_Bone.png',O/'figures/03_markers_Lung.png']:
 if not (backup/f.name).exists():shutil.copyfile(f,backup/f.name)
for node in nodes.node_id:
 col='leaf' if '_s' in node else 'parent';ix=m.index[~m.doublet&l[col].eq(node)];a=[];b=[]
 for _,group in m.loc[ix].groupby(['source','donor']):a.append(X[group.index].mean(0));b.append((X[group.index]>0).mean(0))
 for j,g in enumerate(markers):records.append(dict(node_id=node,gene=g,mean=np.mean(a,axis=0)[j],detection=np.mean(b,axis=0)[j],donor_records=len(a)))
pr=pd.DataFrame(records);pr.to_csv(O/'tables/marker_dotplot_values.tsv',sep='\t',index=False)
for site in ['Bone','Lung']:
 groups=nodes[nodes.site.eq(site)&nodes.level.eq('parent')].node_id.tolist();a=pr.pivot(index='node_id',columns='gene',values='mean').reindex(index=groups,columns=markers);b=pr.pivot(index='node_id',columns='gene',values='detection').reindex(index=groups,columns=markers);z=(a-a.mean())/a.std(ddof=0).replace(0,1)
 fig,ax=plt.subplots(figsize=(13,7));xx,yy=np.meshgrid(range(len(markers)),range(len(groups)));im=ax.scatter(xx.ravel(),yy.ravel(),s=b.to_numpy().ravel()*130,c=z.to_numpy().ravel(),cmap='RdBu_r',norm=Normalize(-2,2));ax.set_xticks(range(len(markers)),markers,rotation=65,ha='right',fontsize=10);ax.set_yticks(range(len(groups)),groups,fontsize=10);ax.set_ylim(len(groups)-.4,-.7);ax.set_title(site+' parent marker panel',pad=45);ax.grid(alpha=.14);ax.set_axisbelow(True);ax.legend([ax.scatter([],[],s=v*130,c='gray') for v in [.1,.5,1]],['10%','50%','100%'],title='RNA detection',ncol=3,loc='lower left',bbox_to_anchor=(0,1.01),frameon=False,fontsize=9);fig.colorbar(im,ax=ax,label='Gene-wise z-score of mean log1p(CP10k)',shrink=.65);fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(O/'figures'/('03_markers_'+site+'.'+ext),dpi=300,bbox_inches='tight')
 plt.close(fig)
note=dict(status='PASS',reason='Use equal donor-record weights for descriptive dotplots; avoid one-cell secondary source receiving half the weight of NCI232 cells in Lung_04',scope='Visualization only; no change to source-stratified effects, P, FDR, labels, or recommendation',method='Within donor mean/detection then equal mean of source-local donor records; no claim of cross-source global patient deduplication')
(O/'resources/DOTPLOT_WEIGHTING_UPDATE.json').write_text(json.dumps(note,indent=2),encoding='utf-8')
print(json.dumps(note))
