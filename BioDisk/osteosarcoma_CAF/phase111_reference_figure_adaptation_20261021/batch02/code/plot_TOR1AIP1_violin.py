from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

P=Path('C:/Users/31937/CAF_paper_figures_20261021')
O=Path('C:/Users/31937/Downloads/CAF_reference_figures_600DPI_20261021')
O.mkdir(exist_ok=True,parents=True)
name='Fig1h_TOR1AIP1_patient_expression_violin'
f=P/'Fig1i_cell_values.tsv.gz';d=pd.read_csv(f,sep='\t',encoding='utf-8')
assert len(d)==6633 and d.cell_id.is_unique and d.donor.nunique()==39
assert np.isfinite(d['TOR1AIP1__log']).all() and (d['TOR1AIP1__raw']>=0).all()
err=float(np.max(np.abs(np.log1p(d['TOR1AIP1__raw']/d.full_library_counts*1e4)-d['TOR1AIP1__log'])))
assert err<1e-12
v=d.groupby(['source','class','donor'],sort=True).agg(cells=('cell_id','size'),mean_log1pCP10k=('TOR1AIP1__log','mean'),raw_TOR1AIP1_counts=('TOR1AIP1__raw','sum')).reset_index()
v['gene']='TOR1AIP1'
v.to_csv(P/(name+'_patient_values.tsv'),sep='\t',index=False,encoding='utf-8')
sources=['SCPCP000017','SCPCP000023','NCI']
groups=[(s,c) for s in sources for c in sorted(v.loc[v.source.eq(s),'class'].unique())]
prior=pd.read_csv(P/'Fig1i_six_candidates_dotplot.tsv',sep='\t',encoding='utf-8')
summary=[]
for s,c in groups:
 q=v[(v.source==s)&(v['class']==c)];x=q.mean_log1pCP10k.to_numpy()
 old=prior[(prior.source==s)&(prior.group==c)&prior.gene.eq('TOR1AIP1')].iloc[0]
 assert len(q)==old.donors and q.cells.sum()==old.cells and abs(x.mean()-old.mean_log1pCP10k)<1e-12
 summary.append(dict(source=s,group=c,patients=len(q),cells=int(q.cells.sum()),mean=float(x.mean()),median=float(np.median(x)),minimum=float(x.min()),maximum=float(x.max()),zero_patients=int((x==0).sum()),density_drawn=bool(len(q)>=3 and np.ptp(x)>1e-12)))
pd.DataFrame(summary).to_csv(P/(name+'_summary.tsv'),sep='\t',index=False,encoding='utf-8')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.8,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(9.6,6.6),dpi=140)
ax=fig.add_axes([.10,.27,.86,.50])
colors={'F1':'#477DA8','F2':'#5A9B86','F3':'#9D82BC','F4':'#D27A79','F5':'#C3A353','Author CAF':'#84929E'}
rng=np.random.default_rng(20261021)
for i,((src,cl),st) in enumerate(zip(groups,summary)):
 q=v[(v.source==src)&(v['class']==cl)];x=q.mean_log1pCP10k.to_numpy();color=colors[cl]
 if st['density_drawn']:
  kde=gaussian_kde(x,bw_method='scott');grid=np.linspace(0,x.max(),350)
  density=kde(grid)+kde(-grid);width=.34*density/density.max()
  ax.fill_betweenx(grid,i-width,i+width,facecolor=color,edgecolor=color,alpha=.36,lw=.85,zorder=2)
 jitter=rng.uniform(-.14,.14,len(x)) if len(x)>1 else np.zeros(1)
 ax.scatter(i+jitter,x,s=25,c=color,edgecolors='white',linewidths=.5,zorder=4)
 median=float(np.median(x));ax.plot([i-.20,i+.20],[median,median],color='#25313A',lw=1.4,zorder=5)
 ax.text(i,1.025,f"n={len(x)}",transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=9)
ax.axvline(4.5,color='#BBC2C8',lw=.8);ax.axvline(8.5,color='#BBC2C8',lw=.8)
ax.set_xlim(-.6,9.6);ax.set_ylim(-0.07,2.12)
ax.set_xticks(range(10),[c if c!='Author CAF' else 'Author CAF' for _,c in groups],fontsize=10)
ax.set_yticks(np.arange(0,2.05,0.5));ax.set_ylabel('Patient mean TOR1AIP1 expression\nlog1p(CP10k)',labelpad=12,fontsize=11)
ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='x',length=0,pad=10)
ax.grid(axis='y',color='#E8EBEE',lw=.65,zorder=0)
for left,right,src in [(0,4,'SCPCP000017'),(5,8,'SCPCP000023'),(9,9,'NCI')]:
 ax.text((left+right)/2,-.18,src,transform=ax.get_xaxis_transform(),ha='center',fontsize=10,color='#53616B')
fig.text(.10,.93,'TOR1AIP1 expression across CAF candidate groups',fontsize=16,fontweight='bold')
fig.text(.10,.875,'Patient-level distributions within each source and existing class',fontsize=11,color='#53616B')
handles=[Line2D([0],[0],marker='o',ls='',markerfacecolor='#84929E',markeredgecolor='white',markersize=6,label='One patient'),Line2D([0],[0],color='#25313A',lw=1.4,label='Median')]
fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.09,.092),frameon=False,ncol=2,fontsize=10)
fig.text(.10,.065,'Density shown only for groups with >=3 patients and nonzero variation.',fontsize=9,color='#53616B')
fig.text(.10,.027,'NCI: author-annotated CAFs; original F classes: expanded candidates. Descriptive, no P values.',fontsize=9,color='#53616B')
for ext in ['png','tiff']:
 kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
 fig.savefig(O/(name+'.'+ext),dpi=600,**kw)
fig.savefig(P/(name+'_preview.png'),dpi=140);plt.close(fig)
exports=[]
for ext in ['png','tiff']:
 path=O/(name+'.'+ext)
 with Image.open(path) as im:
  im.load();assert im.size==(5760,3960);dpi=tuple(float(z) for z in im.info['dpi']);assert max(abs(z-600) for z in dpi)<.05
  exports.append(dict(file=str(path),pixels=im.size,dpi=dpi,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
qa=dict(status='data_and_export_PASS_visual_pending',gene='TOR1AIP1',cell_count=len(d),source_resolved_patient_records=d.donor.nunique(),patient_class_records=len(v),groups=10,input_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),normalization_max_error=err,previous_dotplot_group_means_reproduced=True,density_method='Scott Gaussian KDE with zero-boundary reflection, displayed within observed range; equal-width normalization; only >=3 patients with variable expression',zero_values_retained=True,inference='descriptive only; no inferential testing across overlapping patient groups',exports=exports)
(P/(name+'_QA.json')).write_text(json.dumps(qa,indent=2),encoding='utf-8')
caption='''TOR1AIP1在CAF候选分类中的患者级表达分布（参考原文Fig.1h图型）

每个点代表一名患者在对应来源/分类内的平均TOR1AIP1表达：先按每个细胞的完整转录组UMI归一化为CP10k，再取自然对数log1p，最后求患者内均值。小提琴表示患者均值的分布，各患者等权；深色短横线为中位数。n为该分类的患者数，不是细胞数。同一患者可贡献多个分类，故各行/组n不能直接相加。

使用与上一张六候选气泡图相同的6633个细胞、39个来源内患者记录。源数据为Phase103扩展候选池和NCI作者CAF。来源保持分开，原始F1–F5分类未重新定义；NCI使用作者CAF注释，不与某个F分类等同。SCPCP000023未纳入F3细胞，不显示空组。

小提琴基于患者均值而非把单细胞当生物学重复。使用Scott带宽高斯核密度，零边界反射并限定在观察范围；各小提琴最大宽度相同，宽度不表示样本量。只有至少3名患者且表达有变异时绘制密度。其余组保留全部实际点和中位数，不绘制人为密度。横向抖动只用于显示重叠点，不改变表达值。

本图为描述性结果，不进行跨组显著性检验，不产生P值。零值保留。少数患者、少量细胞组不能据此作总体推断。TOR1AIP1 RNA不能解释为流式或蛋白阳性。原始来源保留扩展候选池，不声称全部细胞为已确认纯CAF；NCI逐文库合并前基因测量掩码缺失的边界沿用上一张。来源、技术及病灶差异均可能影响表达，不能由本图单独认定亚群特异性或功能。

核对：患者级组均值、患者数、细胞数均与上一张气泡图一致；由原始计数重建的log1p(CP10k)通过数值核对。数据、代码和验收记录保存于sgs10的phase111_reference_figure_adaptation_20261021。
各基因图采用覆盖自身观测范围的纵轴；比较基因时需读取刻度，不直接比较图中高度。
输出：同一张单图的PNG与LZW TIFF，5760×3960像素，600 DPI，9.6×6.6英寸。
'''
(O/'Fig1h_TOR1AIP1_图注与口径.txt').write_text(caption,encoding='utf-8')
print(json.dumps(qa));print(pd.DataFrame(summary).to_string(index=False))
