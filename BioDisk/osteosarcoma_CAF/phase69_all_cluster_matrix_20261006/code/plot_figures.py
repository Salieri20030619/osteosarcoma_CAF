
from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
W=Path(r'C:\Users\31937\CAF_all_clusters_20261006');D=W/'analysis/tables';F=W/'figures';F.mkdir(exist_ok=True)
font_manager.fontManager.addfont(r'C:\Windows\Fonts\msyh.ttc')
plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':10,'savefig.dpi':200})
s=pd.read_csv(D/'subgroup_summary.tsv',sep='\t');order=s.subgroup_id.tolist()
p=pd.read_csv(D/'all_source_program_statistics.tsv',sep='\t')
g=pd.read_csv(D/'inherited_all_gene_statistics.tsv.gz',sep='\t')
q=pd.read_csv(D/'additional_nonprimary_source_signals.tsv',sep='\t')
short=lambda x:x.replace('_original','-').replace('_','-')
def save(fig,name):
 fig.savefig(F/(name+'.png'),bbox_inches='tight',facecolor='white');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
programs=['collagen_biosynthesis','collagen_assembly','collagen_crosslinking','collagen_components','ECM_organization','fibrillar_I_III_V','collagen_IV']
names=['胶原合成','胶原装配','交联相关','胶原组分','基质组织','I/III/V 型','IV 型']
a=p[(p.matrix=='corrected')&(p.scope=='all')&(p.minimum_cells==1)&(p.version=='full')&p.is_primary_source&p.program.isin(programs)]
a=a.set_index(['subgroup_id','program'])
values=np.array([[a.loc[(st,pr),'mean_effect_excluded'] for pr in programs] for st in order])
fig,ax=plt.subplots(figsize=(9.3,10.7))
cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
im=ax.imshow(values,cmap=cmap,vmin=-2,vmax=2,aspect='auto')
labels=[]
for st in order:
 n=int(a.loc[(st,programs[0]),'N_excluded']);ds=s.set_index('subgroup_id').loc[st,'primary_source'].replace('SCPCP0000','S')
 labels.append(f'{short(st)}  |  {ds}  |  N={n}')
ax.set_yticks(range(33),labels,fontsize=9);ax.set_xticks(range(7),names,rotation=30,ha='right')
for i,st in enumerate(order):
 for j,pr in enumerate(programs):
  r=a.loc[(st,pr)];txt=f'{values[i,j]:+.2f}'
  if pd.notna(r.P_excluded) and r.P_excluded<.05:txt+='*'
  ax.text(j,i,txt,ha='center',va='center',fontsize=8,color='white' if abs(values[i,j])>1.2 else 'black')
ax.set_title('全部 33 簇：主要来源内的患者级程序效应',loc='left',pad=13,fontweight='bold')
cb=fig.colorbar(im,ax=ax,shrink=.45,pad=.02);cb.set_label('平均程序效应（继承原程序评分尺度）')
fig.text(.02,.015,'S17/S18/S23 = SCPCP000017/18/23；* 原始 P < 0.05，不代表 FDR 显著。N = 有效患者数。',fontsize=9)
fig.subplots_adjust(left=.30,bottom=.095,right=.92,top=.95)
save(fig,'01_all33_programs')
# Source heterogeneity for the seven distinct subgroup-gene pairs added outside primary source.
pairs=list(dict.fromkeys(zip(q.subgroup_id,q.gene)))
sources=['GSE152048','GSE162454','GSE250015','GSE270231','GSE318840','SCPCP000017','SCPCP000018','SCPCP000023']
mat=np.full((len(pairs),len(sources)),np.nan);ann={}
for i,(st,ge) in enumerate(pairs):
 for j,ds in enumerate(sources):
  x=g[(g.subgroup_id==st)&(g.gene==ge)&(g.dataset==ds)&(g.matrix=='corrected')]
  if len(x) and x.iloc[0].N_excluded>0:
   r=x.iloc[0];mat[i,j]=r.mean_effect_excluded;ann[i,j]=f'{r.mean_effect_excluded:+.2f}'+('*' if r.P_excluded<.05 else '')+f'\nN={int(r.N_excluded)}'
fig,ax=plt.subplots(figsize=(10,4.4));im=ax.imshow(mat,cmap=cmap,vmin=-4,vmax=4,aspect='auto')
for i in range(len(pairs)):
 for j in range(len(sources)):ax.text(j,i,ann.get((i,j),'—'),ha='center',va='center',fontsize=8,color='white' if abs(mat[i,j])>2.7 else 'black')
ax.set_yticks(range(len(pairs)),[short(st)+' / '+ge for st,ge in pairs],fontsize=9)
ax.set_xticks(range(len(sources)),[x.replace('SCPCP0000','S') for x in sources],rotation=30,ha='right',fontsize=8)
ax.set_title('新增候选的来源差异：同一基因不强求跨来源显著',loc='left',fontweight='bold')
fig.colorbar(im,ax=ax,pad=.015,shrink=.8,label='患者平均 log2 表达比')
fig.tight_layout();save(fig,'02_new_signals_sources')
# Gene expression display at fixed primary source; no significance encoding.
e=pd.read_csv(D/'all_source_expression.tsv',sep='\t')
panel=['COL1A1','COL1A2','COL3A1','COL5A1','COL6A3','COL4A1','ADAMTS14','PCOLCE','LOXL2','FN1','SPARC','BGN','SDC2','TNC','VCAN','ACAN']
a=e[(e.matrix=='corrected')&e.is_primary_source&e.gene.isin(panel)].set_index(['subgroup_id','gene'])
fig,ax=plt.subplots(figsize=(10.2,10.5))
for i,st in enumerate(order):
 for j,ge in enumerate(panel):
  r=a.loc[(st,ge)]
  if r.assayed_cells:
   sc=ax.scatter(j,i,s=6+r.detection_ge1_pct*.95,c=[r.mean_log1p_CP10k],cmap='viridis',vmin=0,vmax=6,edgecolors='none')
ax.set_xlim(-.7,len(panel)-.3);ax.set_ylim(32.6,-.6);ax.set_yticks(range(33),[short(st) for st in order],fontsize=10)
ax.set_xticks(range(len(panel)),panel,rotation=60,ha='right',fontsize=9)
ax.set_axisbelow(True);ax.grid(color='#eeeeee',linewidth=.7)
ax.set_title('全部 33 簇：胶原与基质相关基因的表达分布',loc='left',pad=15,fontweight='bold')
cb=fig.colorbar(sc,ax=ax,pad=.02,shrink=.45);cb.set_label('平均 log1p(CP10k)')
for pct in [10,50,100]:ax.scatter([],[],s=6+pct*.95,c='#667788',label=f'{pct}%')
ax.legend(title='校正计数 ≥ 1',loc='upper left',bbox_to_anchor=(1.03,1),frameon=False,fontsize=9)
fig.text(.02,.012,'固定各簇原主要来源；点大小为检测比例，颜色为表达量。跨来源色深不用于推断簇间差异。',fontsize=9)
fig.subplots_adjust(left=.15,right=.84,bottom=.115,top=.95)
save(fig,'03_all33_expression')
# Sensitivity of all nine new source-specific rows.
cols=['raw_mean_effect','mean_effect_excluded','min5_mean_effect','fibro_mean_effect']
values=q[cols].to_numpy()
fig,ax=plt.subplots(figsize=(9.2,4.6));im=ax.imshow(values,cmap=cmap,vmin=-3,vmax=3,aspect='auto')
for i,r in q.iterrows():
 for j,c in enumerate(cols):
  n=r.N_excluded if j<2 else r.min5_N if j==2 else r.fibro_N
  value=r[c];txt=f'{value:+.2f}\nN={int(n)}' if n and pd.notna(value) else '不可评估'
  ax.text(j,i,txt,ha='center',va='center',fontsize=9,color='white' if pd.notna(value) and abs(value)>2 else 'black')
ax.set_yticks(range(len(q)),[f'{short(r.subgroup_id)} / {r.gene} / '+r.dataset.replace('SCPCP0000','S') for r in q.itertuples()],fontsize=9)
ax.set_xticks(range(4),['原始表达','校正表达','每侧 ≥ 5 细胞','仅 Fibro'],fontsize=10)
ax.set_title('新增信号的敏感性：方向与有效患者数',loc='left',fontweight='bold')
fig.colorbar(im,ax=ax,pad=.02,shrink=.8,label='患者平均 log2 表达比')
fig.tight_layout();save(fig,'04_new_signals_sensitivity')
(W/'FIGURE_QA.json').write_text(json.dumps({'figures':4,'all33_present':len(order)==33,'new_rows':len(q),'source_gene_pairs':len(pairs),'new_P_computed':False},indent=2),encoding='utf-8')
print('FIGURES_COMPLETE')
