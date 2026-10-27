from pathlib import Path
import json,hashlib
import anndata as ad,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase17_all_cohort_expanded_20260921')
palette=list(plt.get_cmap('tab20').colors)+list(plt.get_cmap('tab20b').colors)
for layer in ['supported','expanded']:
 a=ad.read_h5ad(O/'objects'/f'{layer}_integrated.h5ad',backed='r');xy=a.obsm['X_umap'].copy();obs=a.obs.copy();a.file.close()
 cols=['cluster','dataset','myCAF_customer_positive' if layer=='supported' else 'CNV_state']
 fig,axs=plt.subplots(1,3,figsize=(18,5.5))
 for ax,col in zip(axs,cols):
  values=obs[col].astype(str);names=sorted(values.unique(),key=int) if col=='cluster' else sorted(values.unique())
  for i,name in enumerate(names):
   mask=values.eq(name).to_numpy();ax.scatter(xy[mask,0],xy[mask,1],s=1,alpha=.65,color=palette[i],label=name,rasterized=True)
   if col=='cluster':
    xx,yy=np.median(xy[mask],axis=0);ax.text(xx,yy,name,fontsize=7,ha='center',va='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.75,pad=1))
  ax.set_title({'cluster':'Exploratory clusters','dataset':'Data source','myCAF_customer_positive':'Raw myCAF panel','CNV_state':'CNV evidence'}[col]);ax.set_xticks([]);ax.set_yticks([])
  ax.legend(markerscale=5,fontsize=6,loc='upper left',bbox_to_anchor=(1.01,1),ncol=2 if col=='cluster' else 1,frameon=False)
 fig.tight_layout();fig.savefig(O/'figures'/f'{layer}_joint_umap.png',dpi=180);plt.close(fig)
report=O/'delivery/全队列拓展探索结果.md'
t=report.read_text();anchor='## 全队列联合分析已形成探索性表达群'
summary='## 首轮优先查看支持层群7和群9，并单列血管壁富集群3\n\n支持层群7共2,509个细胞，原始面板myCAF为2,240个，血管壁标志比例7.3%；可与严格子集比较的3条供者记录中2条同方向，myCAF程序群内减群外中位效应为+0.548/+0.285。它提供了可继续复核的收缩程序候选线索，但方向并非全部一致。\n\n群9共3,754个细胞，其中2,391个来自严格层，原始面板myCAF为2,686个，血管壁标志比例0.7%；然而8条可比较记录只有4条同方向，说明不能将此群的myCAF程序直接视为稳定的新状态。群3有6,379个细胞，原始面板myCAF为4,807个，但71.7%带血管壁程序，应作为血管壁/收缩身份边界群单独复核。\n\n'
if summary not in t:report.write_text(t.replace(anchor,summary+anchor))
h={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['tables','figures','delivery','code'] for p in (O/folder).glob('*') if p.is_file() and p.name!='SHA256.json'}
(O/'delivery/SHA256.json').write_text(json.dumps(h,indent=2))
print('FIGURES_AND_RESULT_INTERPRETATION_COMPLETE')
