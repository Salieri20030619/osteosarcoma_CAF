from pathlib import Path
import pandas as pd,numpy as np,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase61_F3_internal_portraits_20261005'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def tab(d):return d.to_csv(sep='\t',index=False,float_format='%.6g',na_rep='NA')
assert (O/'QA_COMPLETE').exists();c=read(O/'tables/primary_245_cells.tsv.gz');allc=read(O/'tables/F3_cells.tsv.gz');x=read(O/'tables/all_GSEA_results.tsv.gz');marker=read(O/'tables/descriptive_marker_display.tsv');zf=np.load(O/'objects/primary_full_transcriptome.npz');z={k:zf[k] for k in zf.files};zf.close();gidx={g:i for i,g in enumerate(z['symbols'])};states=['C0','C1','C2'];state_colors={'C0':'#277DA1','C1':'#D2823A','C2':'#8A5A9E'};plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False});v=np.log1p(z['corrected']/z['corrected_library'][:,None]*1e4)
genes=['FBLN2','TFPI2','COL4A1','COL4A2','IL1R1','ISLR','CLU','COL14A1','CTHRC1','LOX','GJA1','MFAP2','ANXA2','LUM','CNN3','STAT1','MX1'];idx=[gidx[g] for g in genes];means=np.vstack([v[c.state.eq(st)][:,idx].mean(0) for st in states]);frac=np.vstack([(z['corrected'][c.state.eq(st)][:,idx]>=.5).mean(0) for st in states]);scaled=(means-means.mean(0))/np.maximum(means.std(0),1e-9)
fig,ax=plt.subplots(figsize=(12,3.6),layout='constrained')
for i,st in enumerate(states):im=ax.scatter(np.arange(len(genes)),np.repeat(i,len(genes)),s=10+frac[i]*210,c=scaled[i],cmap='RdBu_r',vmin=-1.5,vmax=1.5,edgecolors='#555',linewidth=.3)
ax.set_xticks(range(len(genes)),genes,rotation=50,ha='right');ax.set_yticks(range(3),[f'{s} (n={int((c.state==s).sum())} cells)' for s in states]);ax.invert_yaxis();ax.set_title('F3 internal marker distributions',weight='bold');fig.colorbar(im,ax=ax,label='Gene-scaled mean log1p(CP10k)');fig.supxlabel('Dot size: fraction detected at corrected count >=0.5. Descriptive markers, not patient-level significant DE genes.',fontsize=9);fig.savefig(O/'figures/01_internal_marker_distribution.png',dpi=180);plt.close(fig)
terms=['HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_HYPOXIA','HALLMARK_ANGIOGENESIS','HALLMARK_MYC_TARGETS_V1','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','REACTOME_TRANSLATION','REACTOME_COLLAGEN_FORMATION','NABA_CORE_MATRISOME','NABA_ECM_GLYCOPROTEINS','NABA_BASEMENT_MEMBRANES'];contexts=['corrected__'+s+'_vs_other_F3' for s in states]+['corrected__'+s+'_Fibro_vs_other_F3_Fibro' for s in states];d=x[x.context.isin(contexts)&x.Term.isin(terms)];d.to_csv(O/'tables/key_pathway_evidence.tsv',sep='\t',index=False);m=d.pivot(index='Term',columns='context',values='NES').reindex(index=terms,columns=contexts)
fig,ax=plt.subplots(figsize=(11.5,7),layout='constrained');im=ax.imshow(m,cmap='RdBu_r',vmin=-2.6,vmax=2.6,aspect='auto');ax.set_yticks(range(len(terms)),[t.replace('HALLMARK_','').replace('REACTOME_','').replace('NABA_','').replace('_',' ') for t in terms]);ax.set_xticks(range(6),['C0 all','C1 all','C2 all','C0 Fibro','C1 Fibro','C2 Fibro']);ax.axvline(2.5,color='white',lw=3)
for i,t in enumerate(terms):
 for j,ctx in enumerate(contexts):
  r=d[(d.Term==t)&(d.context==ctx)].iloc[0];star='**' if r.model_BH_FDR<.05 else '*' if r.GSEA_model_P<.05 else '';ax.text(j,i,f'{r.NES:+.2f}{star}',ha='center',va='center',fontsize=8,color='white' if abs(r.NES)>1.6 else 'black')
ax.set_title('Internal functional programs and lineage sensitivity',weight='bold');fig.colorbar(im,ax=ax,label='GSEA NES');fig.supxlabel('* Model P <0.05; ** model BH FDR <0.05. Gene-set permutation evidence, not biological-patient inference.',fontsize=9);fig.savefig(O/'figures/02_internal_programs.png',dpi=180);plt.close(fig)
pair=['CLU','ISLR','LOX','COL14A1','CTHRC1','MFAP2','ANXA2','STAT1','MX1'];pg=[];fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
for j,(ctx,label,color,offset) in enumerate([('corrected__C1_vs_C2_same_patient','All cells','#277DA1',-.12),('corrected__C1_vs_C2_same_patient_Fibro','Fibro only','#D2823A',.12)]):
 g=read(O/'tables/genes'/f'{ctx}.tsv.gz').set_index('symbol').loc[pair];g['context']=ctx;pg.append(g.reset_index());ax.scatter(g.mean_delta,np.arange(len(pair))+offset,label=label,c=color,s=45)
ax.set_yticks(range(len(pair)),pair);ax.invert_yaxis();ax.axvline(0,color='#999',lw=.8);ax.set_xlabel('C1 - C2 pseudobulk log2(CP10k + 0.1) difference');ax.set_title('Within one patient: C1 and C2 share matrix markers but differ in state',weight='bold');ax.legend(frameon=False);fig.supxlabel('17-439_00067: 40 vs84 cells; Fibro subset29 vs54. Descriptive comparison within one biological sample.',fontsize=9);fig.savefig(O/'figures/03_C1_C2_same_patient.png',dpi=180);plt.close(fig);pd.concat(pg).to_csv(O/'tables/C1_C2_marker_evidence.tsv',sep='\t',index=False)
fig,aa=plt.subplots(2,3,figsize=(11.5,7),layout='constrained')
for st in states:
 q=c[c.state==st];aa[0,0].scatter(q.UMAP1,q.UMAP2,c=state_colors[st],label=st,s=17)
aa[0,0].legend(frameon=False);aa[0,0].set_title('Frozen internal states')
for ax,gene in zip(aa.flat[1:],['FBLN2','COL4A1','CLU','LOX','MFAP2']):
 value=v[:,gidx[gene]];im=ax.scatter(c.UMAP1,c.UMAP2,c=value,cmap='viridis',s=17,vmin=0,vmax=max(float(np.quantile(value,.98)),.1));ax.set_title(gene);fig.colorbar(im,ax=ax,shrink=.7)
for ax in aa.flat:ax.set_xticks([]);ax.set_yticks([])
fig.suptitle('F3 internal expression distribution | SCPCP000017',weight='bold');fig.supxlabel('Corrected log1p(CP10k), gene-specific color scales; frozen UMAP. Points are cells, not independent statistical samples.',fontsize=9);fig.savefig(O/'figures/04_internal_feature_maps.png',dpi=180);plt.close(fig)
# Traceable working labels, with unchanged original state identifiers.
labels={'C0':'FBLN2_COL4_inflammatory_response_state','C1':'ISLR_CLU_LOX_matrix_related_state','C2':'MFAP2_ANXA2_MYC_translation_state'};allc['internal_working_description']=allc.state.map(labels);allc.to_csv(O/'tables/F3_internal_working_labels.tsv.gz',sep='\t',index=False)
summary=[]
for st in states:
 ctx='corrected__'+st+'_vs_other_F3';q=x[x.context.eq(ctx)];r=marker[marker.context.eq(ctx)];summary.append({'state':st,'all_cells':int(allc.state.eq(st).sum()),'primary_cells':int(c.state.eq(st).sum()),'description':labels[st],'display_genes':r.symbol.tolist(),'nominal_up_pathways':int(((q.GSEA_model_P<.05)&(q.NES>0)).sum()),'model_FDR_up_pathways':int(((q.model_BH_FDR<.05)&(q.NES>0)).sum())})
pd.DataFrame(summary).to_csv(O/'tables/state_portrait_summary.tsv',sep='\t',index=False)
text='''CAF项目 F3内部亚群表达与功能画像
2026-10-05 phase61

研究结论

F1–F5继续作为粗粒度工作标签。本轮直接刻画F3内部C0/C1/C2，得到三个可用的描述性亚群：C0为FBLN2/COL4相关、偏炎症应答的状态；C1为ISLR/CLU/LOX相关基质状态；C2为MFAP2/ANXA2相关、偏MYC和翻译的状态。C1/C2共享明显的基质表达背景，其内部差别更像状态偏向，不必强行对应两种标准CAF类型。

本轮以内部表达分布、全转录组排序和完整通路为主线；与其他F类的特异性和CAF正式命名不作为分析前置条件。保留原始C0/C1/C2编号，仅新增工作描述字段。

方法

沿用F3全部251个细胞和冻结的内部状态。C0/C1/C2全来源分别119/42/90个；主分析使用同一来源SCPCP000017的245个细胞，对应115/41/89个。每个内部群分别与其余F3细胞比较；另进行DecontX后的Fibro限定分析，以及17-439_00067同一患者内C1对C2的比较（40对84个，Fibro限定29对54个）。共12个全转录组对比。

基因表达用完整library归一化。输出汇总计数log2(CP10k+0.1)差、逐细胞log1p(CP10k)均值、检出率差与描述性AUC。三个主样本的canonical基因测量覆盖已逐项核对，缺测没有补成表达零。主要全群排序各含15036个至少在一侧表达的canonical基因；不同范围的排序规模保留在对比登记表。未按P筛选输入基因。

可读标志图按目标群至少10%检出、正表达差与正检出率差展示前12项；此规则用于展示，不是名义显著候选筛选。全部基因保留在完整表中。去污染计数检出用≥0.5，并同时输出>0检出率。

通路使用冻结版本的Hallmark、Reactome和Matrisome全部集合，执行36项GSEA，1000次基因集置换，固定种子20261005。按每个对比和通路库报告模型原始P、native FDR q-value及BH FDR；模型P<0.05可进入探索性通路解释。结果中的0*表示有限置换未观察到相应尾部事件，不能解释为真实概率为零。基因层面没有构造细胞级推断P，也没有将模型P冒充独立患者检验。

C0突出FBLN2和胶原IV相关标志，炎症应答程序更强

C0的代表表达标志为FBLN2、TFPI2、COL4A1、COL4A2、TIMP3、IL1R1。FBLN2在C0检出率74.8%，其余F3为4.6%；归一化汇总表达差+5.15，描述性AUC=0.866。TFPI2检出率51.3%对0.8%，COL4A1为77.4%对7.7%。这些标志可用于在当前对象中辨认C0。

完整Hallmark中IL6-JAK-STAT3、TNFα-NFκB和缺氧分别呈NES+1.97、+1.73、+1.70，均通过模型BH FDR；Fibro限定后方向及模型证据保留。C0因此适合描述为“FBLN2/COL4相关炎症应答状态”。这里的炎症应答不自动等同于分泌炎性因子最多。

虽然胶原IV标志突出，C0并不是全面基质生产最强的一群：核心Matrisome、胶原形成及翻译排序均较低。完整基底膜集合在整群只是正向（NES+1.43，模型P=0.0803），Fibro限定为模型名义信号（P=0.0295，BH FDR=0.0591）。因此名称保留“COL4相关”，不写成已经证实的完整基底膜功能增强。MTCO1P12、MTATP6P1等技术敏感转录本也进入自动表达排名，已保留全表，但不用于赋予亚群生物学功能名称。

C1保留ISLR CLU LOX基质状态，整体基质富集依赖组成

C1的代表标志包括ISLR、CLU、GJA1、COL14A1、CTHRC1、LOX、OLFML3。ISLR检出率68.3%对22.1%，CLU为65.9%对27.0%，LOX为58.5%对17.2%。核心Matrisome、ECM糖蛋白及胶原形成在C1整体比较中分别呈NES+1.60、+1.49、+1.43，均通过模型BH FDR。血管生成仅为模型名义信号（NES+1.63，P=0.00810，BH FDR=0.0578），不把它单独升级为血管生成型身份。

Fibro限定后，C1的ECM糖蛋白集合仍呈正向（NES+1.39），并保留模型BH FDR支持；核心Matrisome整体富集则未保持（NES−1.02，模型P=0.451），胶原形成也没有模型名义支持。该结果支持特定基质组分偏向，不能概括为全面基质活化。C1仍有可识别的ISLR/CLU/LOX表达组合；工作名称以该组合和完整组分结果为依据，不要求它匹配某个文献CAF分类。

C2突出MYC和翻译程序，干扰素应答是次级线索

C2的代表标志包括MFAP2、ANXA2、LUM、CNN3，并与C1共享MFAP4、CTHRC1、COL14A1等基质基因。MFAP2检出率80.9%对29.5%，ANXA2为77.5%对30.8%，LUM为92.1%对38.5%。

MYC V1与翻译的NES分别+2.11、+1.32；Fibro限定后分别+2.31、+1.74，均保持模型BH FDR支持。因此“MYC/翻译偏高”是C2本轮最清楚的内部状态特征。

整群干扰素α应答也呈正向富集（NES+1.71），但Fibro限定后没有模型名义支持（P=0.120）。其作为次级线索记录，不能把全部C2定名为干扰素型CAF。MYC/翻译富集也不等同于整群实际增殖速度更快。

同患者C1与C2显示共享基质背景上的状态差异

在17-439_00067中，C1对C2的CLU、ISLR、LOX表达差分别+1.13、+0.82、+0.44；Fibro限定后仍分别+1.40、+0.83、+0.26。相反，C2的MFAP2、ANXA2略高。CTHRC1在两群均高检出（95.0%与90.5%），COL14A1也均高（90.0%与82.1%），说明二者共享基质表达，而非截然不同的谱系。

同患者通路比较支持C2的MYC V1更高：以C1减C2排序时NES−1.53、模型P=0.00160、BH FDR=0.0267；Fibro限定后仍有模型名义方向（NES−1.33，P=0.0229，BH FDR=0.393）。整群干扰素α/γ也偏C2；Fibro限定仅γ保留模型名义证据。两群核心Matrisome在同患者整群对比中没有模型名义差异（P=0.797）。

因此，C1/C2作为两个内部工作状态可以继续保留：C1偏CLU/ISLR相关基质表达，C2偏MYC/翻译并伴部分应答程序。二者之间无需再强行套用不同CAF大类名称。

内部标志完整展示
'''
text+=tab(marker[marker.context.str.startswith('corrected')][['context','symbol','mean_delta','target_detection','reference_detection','descriptive_AUC']])
text+='\n关键通路模型证据\n'+tab(d[['context','Term','NES','GSEA_model_P','model_BH_FDR','FDR q-val']])
text+='\n同患者C1/C2标志表达\n'+tab(pd.concat(pg)[['context','symbol','mean_delta','target_detection','reference_detection','descriptive_AUC']])
text+='''
讨论与使用口径

本轮描述的是当前F3内部的表达结构与可检验功能偏向。C0主要来自一位患者，C1/C2主要来自另一位；同患者对比已将部分患者差异分离，但仍是一个患者内的描述。此信息集中用于界定外推范围，不作为继续分析内部亚群的门槛。

既有Fibro/Osteo表达复核和原始质量记录作为辅助列保留，不以CAF身份判定替代内部特征分析。上述工作名称不是固定细胞类型本体，不表示已完成蛋白、空间或功能验证。疾病命名的通路库条目仅表示原基因集来源，不能据“melanoma”“viral”等字样赋予本亚群其他疾病属性。

F3整体对其他F群缺少相对特异性，与F3内部存在上述差异是不同问题。本轮不再用F3对其他F类的P值决定是否描述内部亚群。

交付与验收

全部251个原始F3细胞和C0/C1/C2编号保留。tables/F3_internal_working_labels.tsv.gz新增工作描述，可直接用于后续分组。tables/genes包含12份全转录组描述表；tables/gsea及all_GSEA_results.tsv.gz包含36项完整通路结果，共21971行；原始/去污染、Fibro限定和同患者结果均保留。resources包含基因集、排序、比较登记、预设口径及输入指纹。

独立验收重建全部基因效应、直接核对240项描述性AUC，并验证三个主样本的基因测量覆盖。全部通路leading-edge成员核对通过；独立R重算36个模型BH FDR家族，最大差异约1.37e-14。没有生成不适用的独立患者差异P。

图1展示内部标志的表达量与检出比例；图2展示整体及Fibro限定通路；图3展示同患者C1/C2标志差异；图4展示冻结UMAP上的亚群和特征表达分布。
'''
(O/'delivery/F3内部亚群表达与功能画像_20261005.txt').write_text(text,encoding='utf-8-sig');print('REPORT_COMPLETE',flush=True)
