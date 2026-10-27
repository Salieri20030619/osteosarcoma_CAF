from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from threadpoolctl import threadpool_limits
threadpool_limits(4)
p=Path.cwd();root=p/'report_figures_20260917';out=root/'figures';out.mkdir(exist_ok=True);rng=np.random.default_rng(20260917)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white'})
cols=['#327DA8','#DE8650','#6BA484','#8F79A9','#C6A255','#5FA6A5']
E=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
def save(fig,name,caption):
 fig.savefig(out/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white');fig.savefig(out/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig);captions[name]=caption;print('FIGURE',name,flush=True)
def tidy(ax):ax.set(xticks=[],yticks=[],xlabel='UMAP 1',ylabel='UMAP 2');ax.set_aspect('equal',adjustable='datalim')
def catplot(ax,xy,lab,palette,title,size=3):
 for k,c in palette.items():
  ix=np.asarray(lab)==k
  if ix.any():ax.scatter(xy[ix,0],xy[ix,1],s=size,c=c,label=k,rasterized=True,linewidths=0,alpha=.85)
 tidy(ax);ax.set_title(title,pad=12)
def group(e):
 z=pd.Series('Other / unresolved',index=e.index)
 z.loc[e.lineage_reviewed.isin(['T_NK','B_Plasma','Myeloid','DC','Mast'])]='Immune'
 z.loc[e.lineage_reviewed=='Endo']='Endothelial'
 z.loc[e.lineage_reviewed=='OC']='Osteoclast'
 z.loc[e.final_identity=='Mural_SMC']='Mural'
 z.loc[e.final_identity=='Malignant_OS_candidate']='Malignant candidate'
 z.loc[e.caf_strict]='CAF candidate';return z
pal={'Other / unresolved':'#d5d8dc','Immune':'#6BA484','Endothelial':'#8F79A9','Osteoclast':'#C6A255','Mural':'#5FA6A5','CAF candidate':'#327DA8','Malignant candidate':'#D26D5B'}
captions={}
# Representative chromosome heatmap; sampling independent of CNV amplitude.
sample='GSE152048__BC22';e=E[E['sample']==sample].copy();prof=pd.read_csv(root/'data'/f'{sample}_profiles.tsv.gz',sep='\t').set_index('tool_id');bins=pd.read_csv(root/'data'/f'{sample}_bins.tsv',sep='\t');e=e[e.tool_id.isin(prof.index)].copy();e['plot_group']='Other / unresolved'
e.loc[e.cnv_role.str.contains('immune',na=False),'plot_group']='Immune reference';e.loc[e.cnv_role.str.contains('endo',na=False),'plot_group']='Endothelial reference';e.loc[e.caf_strict,'plot_group']='CAF candidate';e.loc[e.final_identity=='Malignant_OS_candidate','plot_group']='Malignant candidate'
groups=['Immune reference','Endothelial reference','CAF candidate','Malignant candidate','Other / unresolved'];colors=['#6BA484','#8F79A9','#327DA8','#D26D5B','#b6bdc5'];selected=[];sizes=[];total=[]
for name in groups:
 ids=e.index[e.plot_group==name].to_numpy();total.append(len(ids));chosen=rng.choice(ids,min(180,len(ids)),replace=False);chosen=sorted(chosen,key=lambda k:e.loc[k,'infercnv_signal']);selected.extend(chosen);sizes.append(len(chosen))
mat=prof.loc[e.loc[selected,'tool_id']].to_numpy(dtype=float);assert mat.shape[1]==len(bins)
fig=plt.figure(figsize=(12.5,7.4));gs=fig.add_gridspec(1,3,width_ratios=[.18,10,.3],wspace=.12);axg=fig.add_subplot(gs[0]);ax=fig.add_subplot(gs[1]);cax=fig.add_subplot(gs[2]);v=np.repeat(np.arange(len(groups)),sizes)
axg.imshow(v[:,None],aspect='auto',cmap=ListedColormap(colors),interpolation='nearest');axg.set_xticks([]);axg.set_yticks(np.cumsum(sizes)-np.array(sizes)/2,labels=[f'{g}\nshown {n} / total {N}' for g,n,N in zip(groups,sizes,total)],fontsize=9)
im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-.35,vmax=.35,interpolation='nearest');chrom=bins.chr.astype(str).str.replace('chr','',regex=False);ticks=[];names=[]
for cc in chrom.drop_duplicates():
 idx=np.where(chrom==cc)[0];ticks.append(idx.mean());names.append(cc)
 if idx[0]>0:ax.axvline(idx[0]-.5,color='#666666',lw=.45)
for edge in np.cumsum(sizes)[:-1]:ax.axhline(edge-.5,color='#333333',lw=.65)
ax.set(yticks=[],xticks=ticks,xticklabels=names,xlabel='Chromosome (100-gene windows within each chromosome)');ax.tick_params(axis='x',labelsize=8);fig.colorbar(im,cax=cax,extend='both').set_label('Median log2 relative expression\n(inferCNV-derived signal)');fig.suptitle('BC22 | Chromosome-scale inferred CNV profiles',x=.58,y=.96,fontsize=15);fig.text(.24,.04,'Red: increased signal   Blue: decreased signal   |   Deterministic stratified display sample',fontsize=10);fig.subplots_adjust(left=.25,right=.93,bottom=.13,top=.87)
e.loc[selected,['tool_id','plot_group','infercnv_signal']].to_csv(root/'data/heatmap_display_cells.tsv',sep='\t');save(fig,'Fig01_BC22_CNV_heatmap','BC22骨病灶的inferCNV衍生染色体窗口热图。每类随机展示至多180个细胞，类内按信号排序，完整类别总数标在左侧；染色体内每100个基因聚合为窗口。颜色为去噪表达的log2相对信号，并非DNA实测拷贝数。身份由既有多方法规则确定；该图为结果展示，不是独立身份验证。')
# Whole-tissue UMAP on identical coordinates.
a=sc.read_h5ad(p.parent/'phase5_nonpdx_qc_20260916/results'/sample/'annotated.h5ad');ee=E.reindex(a.obs_names);xy=a.obsm['X_umap'];fig,axs=plt.subplots(1,3,figsize=(15,5.6));catplot(axs[0],xy,group(ee).to_numpy(),pal,'A  Reviewed identity | BC22');axs[0].legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,frameon=False,fontsize=8,markerscale=3)
sig=ee.infercnv_signal.to_numpy(float);valid=np.isfinite(sig);axs[1].scatter(xy[~valid,0],xy[~valid,1],s=3,c='#dedede',rasterized=True);q=axs[1].scatter(xy[valid,0],xy[valid,1],c=sig[valid],s=3,cmap='magma',vmin=0,vmax=max(.1,np.nanquantile(sig,.99)),rasterized=True);tidy(axs[1]);axs[1].set_title('B  inferCNV signal');fig.colorbar(q,ax=axs[1],shrink=.7,label='Mean absolute log2 signal')
cp=ee.copykat_call.fillna('Not evaluated').replace({'not_evaluable':'Not evaluated','diploid':'Diploid call','aneuploid':'Aneuploid call'});catplot(axs[2],xy,cp.to_numpy(),{'Not evaluated':'#d5d8dc','Diploid call':'#327DA8','Aneuploid call':'#D26D5B'},'C  CopyKAT call');axs[2].legend(loc='upper center',bbox_to_anchor=(.5,-.16),frameon=False,fontsize=9,markerscale=3);fig.subplots_adjust(bottom=.25,wspace=.35);save(fig,'Fig02_BC22_identity_CNV_UMAP','BC22全细胞UMAP，三面板使用相同坐标。A为复核身份；B为inferCNV信号，灰色为未评估；C为CopyKAT判定。信号色阶上限为0.1与本样本99%分位中的较大值，便于阅读。未定身份保留灰色，不将全部非CAF细胞视为恶性。')
# Descriptive distributions, no cell-level P values.
fig,axs=plt.subplots(2,3,figsize=(14,8));vgroup=['Immune ref.','Endo ref.','CAF cand.','Malig. cand.'];vcolors=colors[:4]
for ax,short in zip(axs.flat,['BC10','BC17','BC2','BC5','BC21','BC22']):
 e=E[E['sample']=='GSE152048__'+short];masks=[e.cnv_role.str.contains('immune',na=False),e.cnv_role.str.contains('endo',na=False),e.caf_strict,e.final_identity=='Malignant_OS_candidate']
 for j,mask in enumerate(masks):
  vals=e.loc[mask,'infercnv_signal'].dropna().to_numpy();n=len(vals)
  if n>=2 and np.ptp(vals)>0:
   obj=ax.violinplot(vals,positions=[j],widths=.75,showextrema=False)
   for bdy in obj['bodies']:bdy.set_facecolor(vcolors[j]);bdy.set_edgecolor('#45505b');bdy.set_alpha(.75)
  if n:ax.scatter(j,np.median(vals),c='black',s=14,zorder=3)
  ax.text(j,.97,f'n={n}',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=8)
 ax.set_xticks(range(4),vgroup,rotation=20,ha='right');ax.set_ylim(bottom=0);ax.set_title(short);ax.set_ylabel('Mean absolute log2 signal');ax.grid(axis='y',alpha=.2)
fig.suptitle('Inferred CNV signal by reference and candidate identity',fontsize=15);fig.text(.5,.015,'Cell distributions are descriptive; no cell-level hypothesis tests. Each panel is one tissue sample.',ha='center',fontsize=10);fig.tight_layout(rect=[0,.05,1,.94]);save(fig,'Fig03_CNV_signal_distributions','6个关键组织的inferCNV信号分布，点为中位数，n为该组可评估细胞数。各面板纵轴按样本数据显示，比较时注意刻度。按CNV参与定义的候选身份分组，组间差异不构成独立验证；不添加以细胞作重复的P值。')
# Method agreement tables.
fig,axs=plt.subplots(1,2,figsize=(11,4.9));rr=['low_signal','ambiguous_signal','elevated_signal'];cc=['diploid','aneuploid','not_evaluable']
for ax,short in zip(axs,['BC22','BC10']):
 e=E[E['sample']=='GSE152048__'+short];t=pd.crosstab(e.infercnv_call,e.copykat_call).reindex(index=rr,columns=cc,fill_value=0);im=ax.imshow(np.log1p(t.to_numpy()),cmap='Blues',vmin=0)
 for i in range(3):
  for j in range(3):ax.text(j,i,str(t.iloc[i,j]),ha='center',va='center',color='white' if np.log1p(t.iloc[i,j])>np.log1p(t.to_numpy().max())*.6 else '#253343',fontsize=13)
 ax.set_xticks(range(3),['Diploid','Aneuploid','Not evaluable']);ax.set_yticks(range(3),['Low','Ambiguous','Elevated']);ax.set(xlabel='CopyKAT call',ylabel='inferCNV call',title=f'{short} | method call counts');t.to_csv(root/'data'/f'{short}_method_counts.tsv',sep='\t')
fig.text(.5,.015,'Numbers are cell counts; color is log(1 + count). Unassessed inferCNV cells are outside this table.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.07,1,1]);save(fig,'Fig04_CNV_method_calls','BC22和BC10两种CNV方法的判定交叉表，保留模糊及CopyKAT不可评估类别；数字为细胞数，颜色为log(1+计数)。未获inferCNV判定者未列入。该图呈现一致与分歧，不能将两种表达推断方法视为独立DNA验证。')
# CAF state visualizations.
a=sc.read_h5ad(p/'integration_review/strict_CAF_harmony_review.h5ad');old=sc.read_h5ad(p/'analysis/GSE152048_strict_CAF_states.h5ad');patients=a.obs.patient_id.str.replace('GSE152048:','',regex=False).to_numpy();ps=sorted(set(patients));pp=dict(zip(ps,cols));xy=a.obsm['X_umap_harmony_theta2'];lab=a.obs['harmony_theta2_r0.5'].astype(str);z=np.where(lab=='0','H0 contractile candidate',np.where(lab=='2','H2 inflammatory candidate','Other candidate clusters'));sp={'Other candidate clusters':'#d7dbe0','H0 contractile candidate':'#DE8650','H2 inflammatory candidate':'#327DA8'}
corrected=pd.read_csv(p/'remediation_20260917/state_final_sensitivity/corrected_recluster_cells.tsv.gz',sep='\t',dtype={'corrected_state':str}).set_index('cell_id').reindex(a.obs_names);cl=corrected.corrected_state.to_numpy();csp={'0':'#DE8650','1':'#327DA8','2':'#6BA484','3':'#8F79A9','4':'#C6A255','5':'#5FA6A5'}
fig,axs=plt.subplots(2,2,figsize=(12.5,10));catplot(axs[0,0],old.obsm['X_umap'],patients,pp,'A  Before Harmony | patients',5);catplot(axs[0,1],xy,patients,pp,'B  After Harmony | patients',5);catplot(axs[1,0],xy,z,sp,'C  Frozen candidate states',5);catplot(axs[1,1],xy,cl,csp,'D  Corrected labels on fixed Harmony map',5)
for ax in axs.flat:ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=3,frameon=False,fontsize=8,markerscale=2)
fig.suptitle('Strict CAF candidate set | 2,047 cells from 6 patients',fontsize=15);fig.subplots_adjust(hspace=.55,bottom=.1,top=.9,wspace=.3);save(fig,'Fig05_CAF_Harmony_states_UMAP','2,047个严格CAF候选的UMAP。A为未整合坐标；B为Harmony强度2后的坐标，均按患者着色。C突出冻结H0/H2候选状态；D把去污染后重聚类标签投射到同一Harmony坐标，便于核对边界变化，D不是去污染后重新计算的UMAP。候选簇不等于已确认稳定亚型。')
# Continuous programs on fixed coordinates, corrected expression values.
prog=pd.read_csv(p/'remediation_20260917/program_stability/cell_program_expression.tsv.gz',sep='\t');prog=prog[prog.version=='corrected'].set_index('cell_id').reindex(a.obs_names)
fig,axs=plt.subplots(2,3,figsize=(14,8.4))
for ax,name in zip(axs.flat,['Contractile','Inflammatory','Matrix_remodeling','Antigen_presentation','Osteogenic','Mineral_calcium']):
 v=prog[name].to_numpy();order=np.argsort(v);im=ax.scatter(xy[order,0],xy[order,1],c=v[order],s=5,cmap='viridis',vmin=0,vmax=np.quantile(v,.99),rasterized=True,linewidths=0);tidy(ax);ax.set_title(name.replace('_',' '));fig.colorbar(im,ax=ax,shrink=.72,label='Mean log1p expression')
fig.suptitle('CAF functional programs after DecontX | fixed Harmony coordinates',fontsize=15);fig.tight_layout(rect=[0,0,1,.94]);save(fig,'Fig06_CAF_programs_UMAP','去污染后六个功能程序的连续表达图，使用固定Harmony坐标。程序值为覆盖基因的平均log1p标准化表达，色阶各截至99%分位。各程序基因构成及色阶不同，不比较跨程序颜色强弱；MHC-II或成骨程序表达本身不证明apCAF或成骨型CAF。')
# Patient composition before and after relabeling.
fig,axs=plt.subplots(1,2,figsize=(12,5.3))
for ax,arr,title in [(axs[0],z,'Frozen candidate-state labels'),(axs[1],np.where(cl=='0','H0 contractile candidate',np.where(cl=='1','H2 inflammatory candidate','Other candidate clusters')),'Matched clusters after DecontX')]:
 ct=pd.crosstab(patients,arr).reindex(index=ps,columns=list(sp),fill_value=0);fr=ct.div(ct.sum(1),axis=0);left=np.zeros(len(ct))
 for name,color in sp.items():ax.barh(ps,fr[name],left=left,color=color,label=name);left+=fr[name].to_numpy()
 for i,N in enumerate(ct.sum(1)):ax.text(1.015,i,f'n={N}',va='center',fontsize=9)
 ax.set(xlim=(0,1.17),xlabel='Fraction of strict CAF candidates',title=title);ax.set_xticks([0,.25,.5,.75,1],['0','25%','50%','75%','100%']);ax.invert_yaxis()
handles=[Line2D([0],[0],lw=7,color=c,label=k) for k,c in sp.items()];fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,fontsize=9);fig.tight_layout(rect=[0,.1,1,1]);save(fig,'Fig07_patient_state_composition','同一严格CAF集合中，冻结状态与去污染重聚类匹配簇的患者内比例。分母是每位患者严格CAF候选数，不是全肿瘤细胞数；变化来自标签重分配，不表示细胞在体内增殖。')
# Patient-level effect plots; descriptive/small-side samples distinguished.
d=pd.read_csv(p/'remediation_20260917/program_stability/within_patient_program_stability.tsv',sep='\t');fig,axs=plt.subplots(1,2,figsize=(11,5.2))
for ax,state,program,title,color in zip(axs,['H2','H0'],['Inflammatory','Contractile'],['Inflammatory candidate state','Contractile candidate state'],['#327DA8','#DE8650']):
 q=d[(d.version=='corrected')&(d.design=='reclustered_states')&(d.state==state)&(d.program==program)].copy();q['short']=q.patient.str.replace('GSE152048:','',regex=False);q=q.sort_values('short');ax.axvline(0,color='#999999',lw=1)
 for i,row in enumerate(q.itertuples()):
  x=row.mean_log_expression_difference;ax.plot([row._asdict()['resampling_2.5pct'] if False else q.iloc[i]['resampling_2.5pct'],q.iloc[i]['resampling_97.5pct']],[i,i],color=color,lw=2);ax.scatter(x,i,s=70,edgecolor=color,facecolor=color if row.eligible_primary_20_per_side else 'white',zorder=3)
 ax.set_yticks(range(len(q)),[f'{t}\n{int(n)} / {int(o)} cells' for t,n,o in zip(q.short,q.target_cells,q.other_cells)]);ax.invert_yaxis();ax.set(title=title,xlabel='Within-patient program expression difference');ax.grid(axis='x',alpha=.15)
fig.text(.5,.015,'Filled: both sides >=20 cells. Open: supplementary description. Bars: conditional cell-resampling range, not biological CI.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.1,1,1]);save(fig,'Fig08_patient_program_effects','校正后重新聚类标签下的逐患者炎症/收缩程序差异。实心为两侧均≥20细胞，空心为小细胞数侧补充描述；横线为500次条件细胞重采样的2.5%—97.5%分位范围，不是患者总体置信区间或显著性检验。纵轴同时列出目标/其他CAF数量。')
(root/'captions.json').write_text(json.dumps(captions,ensure_ascii=False,indent=2));(root/'COMPLETE').write_text('REPORT_FIGURES_RENDERED_PENDING_VISUAL_QA\n');print('ALL_FIGURES_RENDERED',len(captions))
