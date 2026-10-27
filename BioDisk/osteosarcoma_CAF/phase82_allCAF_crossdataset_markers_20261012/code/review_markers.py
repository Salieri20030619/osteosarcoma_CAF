from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase82_allCAF_crossdataset_markers_20261012'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
reg=read('MARKER_EVIDENCE_REGISTRY.tsv');ev=read('all_candidate_source_layer_design_statistics.tsv');pat=read('candidate_patient_expression.tsv');det=read('candidate_detection_and_Fibro_patient_expression.tsv');sens=read('candidate_patient_direction_sensitivity.tsv')
assert len(reg)==14 and (O/'logs/screen.exit_code').read_text().strip()=='0'
cols17='SCPCP000017__';cols23='SCPCP000023__'
reg['both_LOPO_direction_retained']=reg[cols17+'LOPO_direction_consistent']&reg[cols23+'LOPO_direction_consistent']
reg['both_doublet_nominal_retained']=(reg[cols17+'doublet_excluded_P']<.05)&(reg[cols23+'doublet_excluded_P']<.05)
reg['both_Fibro_direction_retained']=(reg[cols17+'Fibro_effect']*reg[cols17+'mean_lung_minus_bone']>0)&(reg[cols23+'Fibro_effect']*reg[cols23+'mean_lung_minus_bone']>0)
reg['paired_both_patients_same_direction']=np.where(reg.direction.eq('lung_higher'),reg.paired17_lung_higher_patients.eq(reg.paired17_N),reg.paired17_bone_higher_patients.eq(reg.paired17_N))
reg['all_Fibro_and_pair_directions_retained']=reg.both_Fibro_direction_retained&reg.paired_both_patients_same_direction
reg['interpretation']='Cross-dataset nominal expression marker in original expanded candidate pool; neither source FDR-significant'
save(reg,'MARKER_EVIDENCE_REGISTRY.tsv')
checks=[]
for r in ev[(ev.sensitivity=='all_cells')&(ev.design=='independent')&ev.P.lt(.05)].itertuples():
    z=pat[(pat.layer==r.layer)&(pat.dataset==r.dataset)&(pat.design=='independent')&(pat.symbol==r.symbol)]
    x=z[z.site.eq('Bone')].pseudobulk_log2_CP10k_plus01.dropna().to_numpy();y=z[z.site.eq('Lung')].pseudobulk_log2_CP10k_plus01.dropna().to_numpy();a=np.r_[x,y];obs=abs(y.mean()-x.mean());hits=0;total=0
    for ix in itertools.combinations(range(len(a)),len(y)):
        m=np.zeros(len(a),bool);m[list(ix)]=True;hits+=abs(a[m].mean()-a[~m].mean())>=obs-1e-12;total+=1
    pp=hits/total;checks.append(dict(layer=r.layer,dataset=r.dataset,symbol=r.symbol,P_stored=r.P,P_recomputed=pp,error=abs(pp-r.P)))
assert max(r['error'] for r in checks)<1e-9
pd.DataFrame(checks).to_csv(O/'QA/candidate_exact_P_verification.tsv',sep='\t',index=False)
read('candidate_existing_membrane_annotation.tsv').to_csv(O/'tables/existing_annotation_not_a_complete_membrane_inventory.tsv',sep='\t',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
reg=reg.sort_values(['direction','symbol'],ascending=[False,True]).reset_index(drop=True)
y=np.arange(len(reg));fig,ax=plt.subplots(figsize=(12,8.5))
for shift,ds,color,label in [(-.12,cols17,'#326b9b','SCPCP000017 (bone 12, lung 3)'),(.12,cols23,'#ce6840','SCPCP000023 (bone 4, lung 3)')]:
    ax.hlines(y+shift,reg[ds+'LOPO_effect_min'],reg[ds+'LOPO_effect_max'],color=color,alpha=.65,lw=2)
    ax.scatter(reg[ds+'mean_lung_minus_bone'],y+shift,s=38,color=color,label=label,zorder=3)
ax.axvline(0,color='#555',lw=.8);ax.set_yticks(y,reg.symbol);ax.invert_yaxis();ax.set_xlabel('Mean lung minus bone: log2(CP10k + 0.1)');ax.set_title('14 markers nominally replicated in both datasets',fontsize=16,pad=15);ax.legend(frameon=False,loc='lower right')
fig.text(.015,.035,'All markers: raw P < 0.05 in both datasets, same direction. Source-specific FDR: 0.6883 and 1.0000.\nSegments: leave-one-patient-out effect ranges, not confidence intervals. RGS2-AS1 = original AL390957.1 (same ENSG ID).',fontsize=9)
fig.tight_layout(rect=[0,.085,1,.97]);fig.savefig(O/'figures/01_crossdataset_marker_effects.png',dpi=180);plt.close(fig)

# Display independent evidence types without converting them into an arbitrary score.
mat=[]
for r in reg.itertuples():
    sign=1 if r.direction=='lung_higher' else -1
    def direction(v):return 1 if v*sign>0 else (-1 if v*sign<0 else 0)
    paired=[r.paired17_lung_higher_patients,r.paired17_bone_higher_patients]
    good=paired[0] if sign>0 else paired[1];bad=paired[1] if sign>0 else paired[0]
    paircode=1 if good==r.paired17_N else (-1 if bad==r.paired17_N else .5 if good and bad else 0)
    mat.append([1,1,1 if r.supported17_nominal_same_direction else 0,1 if r.both_doublet_nominal_retained else 0,direction(getattr(r,cols17+'Fibro_effect')),direction(getattr(r,cols23+'Fibro_effect')),paircode])
fig,ax=plt.subplots(figsize=(12,8));cmap=ListedColormap(['#b85050','#ededed','#e5b765','#3f8c86']);norm=BoundaryNorm([-1.5,-.5,.25,.75,1.5],4)
ax.imshow(mat,cmap=cmap,norm=norm,aspect='auto');ax.set_yticks(y,reg.symbol);ax.set_xticks(range(7),['Dataset17\nnominal','Dataset23\nnominal','Supported17\nnominal same','Doublet\nretained','Fibro17\ndirection','Fibro23\ndirection','Paired17\nboth donors'],rotation=30,ha='right');ax.set_title('Evidence agreement and limitations for every candidate',fontsize=15,pad=12)
ax.legend(handles=[Patch(color='#3f8c86',label='Agreement'),Patch(color='#e5b765',label='Mixed paired directions'),Patch(color='#ededed',label='No added support / zero'),Patch(color='#b85050',label='Opposite')],loc='upper center',bbox_to_anchor=(.5,-.2),ncol=2,frameon=False)
fig.text(.015,.015,'Fibro restriction: descriptive only (dataset17 bone 11 / lung 3; dataset23 bone 2 / lung 2). Paired17 N=2.\nColumns reflect different evidence types; supported/expanded layers and sensitivity subsets are not independent validation datasets.',fontsize=9)
fig.tight_layout(rect=[0,.14,1,.97]);fig.savefig(O/'figures/02_candidate_evidence_agreement.png',dpi=180);plt.close(fig)

selected=['RTKN','BEST1','RGS2-AS1','EBF1'];fig,axes=plt.subplots(4,2,figsize=(11,12),sharex=True)
for i,symbol in enumerate(selected):
    for j,ds in enumerate(['SCPCP000017','SCPCP000023']):
        ax=axes[i,j];z=pat[(pat.layer=='expanded')&(pat.dataset==ds)&(pat.design=='independent')&(pat.symbol==symbol)]
        for x,side,color in [(0,'Bone','#326b9b'),(1,'Lung','#ce6840')]:
            v=z[z.site.eq(side)].pseudobulk_log2_CP10k_plus01.to_numpy();xx=np.full(len(v),x)+np.linspace(-.1,.1,len(v));ax.scatter(xx,v,color=color,s=35);ax.hlines(v.mean(),x-.19,x+.19,color='black',lw=1.5)
        r=reg[reg.symbol.eq(symbol)].iloc[0];prefix=ds+'__';ax.set_title(f'{symbol} | {ds}\nP={r[prefix+"P"]:.4f}; FDR={r[prefix+"FDR"]:.3f}',fontsize=10);ax.set_xticks([0,1],['Bone','Lung']);ax.set_ylabel('Patient log2(CP10k + 0.1)');ax.set_xlim(-.4,1.4)
fig.suptitle('Patient-level expression in the two source datasets',fontsize=16)
fig.text(.015,.015,'Each point is one biological patient; black line is group mean. Display selection: three markers with additional supported-layer\nnominal evidence plus the only bone-higher marker. All14 markers and all patient values remain in the accompanying tables.',fontsize=9)
fig.tight_layout(rect=[0,.065,1,.955]);fig.savefig(O/'figures/03_selected_marker_patient_expression.png',dpi=180);plt.close(fig)

summary=dict(primary_markers=len(reg),lung_higher=int(reg.direction.eq('lung_higher').sum()),bone_higher=int(reg.direction.eq('bone_higher').sum()),extra_supported17_nominal=reg.loc[reg.supported17_nominal_same_direction,'symbol'].tolist(),both_Fibro_direction_retained=reg.loc[reg.both_Fibro_direction_retained,'symbol'].tolist(),paired_both_same=reg.loc[reg.paired_both_patients_same_direction,'symbol'].tolist(),all_Fibro_and_pair_directions_retained=reg.loc[reg.all_Fibro_and_pair_directions_retained,'symbol'].tolist(),all_LOPO_directions_retained=bool(reg.both_LOPO_direction_retained.all()),all_doublet_nominal_retained=bool(reg.both_doublet_nominal_retained.all()),independent_candidate_P_checks=len(checks),max_P_error=max(r['error'] for r in checks),source_P_FDR_unchanged=True,new_inferential_tests=0,visual_review_pending=True)
(O/'QA/NUMERIC_ACCEPTANCE.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
report='''整体CAF候选池骨肺部位差异的跨数据集marker筛选
2026年10月12日

两个数据集复现14个同向名义显著marker
在SCPCP000017与SCPCP000023的扩展候选层中，14个marker均达到原始P<0.05并同向变化，13个肺侧升高、1个骨侧升高。全部保留为跨数据集名义复现候选；两个来源的FDR分别为0.688345和1，未达到校正后显著。这里的marker表示整体候选池的部位相关RNA表达标志，不是已经验证的CAF特异身份或流式分选标志。

数据与方法
沿用Phase81全部CAF候选池骨肺比较，支持层与扩展层分别分析，不跨来源合并患者，不把两层当独立数据集。主比较为SCPCP000017骨12/肺3、SCPCP000023骨4/肺3名独立患者。筛选规则为同层、不同数据集、有效P<0.05且效应同向；保留上下调，不以FDR、膜蛋白注释或抗体可用性作入口。
使用冻结的来源canonical特征映射，不根据最小P挑选转录本；原始P/FDR/家族不改动，未生成合并P或在筛选后的14个基因上重算FDR。效应为肺减骨的患者平均log2(CP10k+0.1)差。逐患者剔除、成纤维谱系限定和检测率用于描述证据；细胞不作为独立生物学重复。

全部14个主候选
'''
compact=pd.DataFrame(dict(marker=reg.symbol,direction=reg.direction,effect17=reg[cols17+'mean_lung_minus_bone'],P17=reg[cols17+'P'],FDR17=reg[cols17+'FDR'],effect23=reg[cols23+'mean_lung_minus_bone'],P23=reg[cols23+'P'],FDR23=reg[cols23+'FDR']))
save(compact,'MAIN_14_MARKERS_FOR_REVIEW.tsv');report+=compact.to_string(index=False)+'\n'
report+='\n额外证据与优先复核角度\n'
report+='RTKN、BEST1、RGS2-AS1在来源17的支持层也有同向P<0.05，为同来源候选定义变化的补充支持，不是第三个独立验证来源。三者可优先用于追踪整体候选池表达，但RTKN与BEST1在来源23的Fibro限定集合中均无表达差，RGS2-AS1的两名配对患者方向一正一负、均值与独立组相反。不能据此直接确认为稳定CAF特异marker。\n'
report+='所有14个marker在两个来源逐患者剔除后均保留原效应方向，去双细胞后P<0.05和方向均保留。该证据说明结果不是仅由单个患者或已标记双细胞独自决定，但不是独立验证。\n'
report+='两个来源Fibro限定后均保留同向效应的候选为：'+ '、'.join(summary['both_Fibro_direction_retained'])+'。此处仅为描述性效应，来源17有效骨11/肺3、来源23仅骨2/肺2；后者细胞覆盖稀疏，不能把这一层升级为独立显著验证。其余8个候选在来源23的Fibro限定记录中效应为0；该事实不能抹去原池的名义差异，但限制其CAF谱系归属。\n'
report+='CALHM5在来源17两名配对患者均为肺侧较高，但来源23的Fibro限定效应为0。EBF1在两个独立组均骨侧较高，却在两名配对患者均肺侧较高，属于关键方向矛盾。当前14个中没有同时得到两来源Fibro方向和两名配对患者一致方向支持的候选；该检查用于解释证据，不作为删除14个P<0.05候选的新门槛。\n'
secondary=read('secondary_cross_layer_replication.tsv');conflicts=read('crossdataset_nominal_direction_conflicts.tsv')
report+='\n跨层补充与矛盾\n支持层来源17对扩展层来源23还得到9条同向名义复现记录，其中6个是主14名单之外的补充候选：HGF、VLDLR、MSRB2、LETM1、SFXN4、CCDC137。因候选层不同，单列而不混入同层主名单。\n'+secondary.to_string(index=False)+'\n'
report+='TCF7L2在上述跨层比较中虽两边P<0.05但方向相反，不作为同向复现marker。\n'+conflicts.to_string(index=False)+'\n'
report+='\n名称核对\n项目规范映射将ENSG00000285280记为RGS2-AS1，Phase81原表符号为AL390957.1；二者是同一特征，P和效应未改变。完整特征及原名称对照见stable_feature_and_symbol_crosswalk.tsv。\n'
report+='\n讨论\n扩展候选池保留身份未解决及成骨样细胞；整个候选池的名义marker尚不能等同于成纤维细胞特异marker。跨数据集原始P重复、患者效应方向、配对矛盾及谱系敏感性分别解释。不同项目的样本编号和记录患者编号未见重叠，但仅凭编号不宣称独立招募已完全验证。本轮未作新的公共数据下载、独立实验验证或抗体选择。既有受体/膜注释表只是有限覆盖，缺失注释不等于非膜蛋白；不据此排除表达候选。\n'
report+='\n交付文件\nMARKER_EVIDENCE_REGISTRY.tsv汇总14个主候选的原P/FDR/家族及各类支持；all_candidate_source_layer_design_statistics.tsv保留主、补充与反向候选在全部来源/层/配对/去双细胞条件下的统计；患者表达、检测率、Fibro限定效应、名称对应及全部阴性证据另表保存。三图覆盖全部14个效应、逐项证据一致性，以及3个额外支持层候选与唯一骨侧高候选的逐患者表达。\n'
(O/'delivery/整体CAF_跨数据集marker筛选_20261012.txt').write_text(report,encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print('FIBRO_COVERAGE\n'+det[(det.layer=='expanded')&(det.dataset=='SCPCP000023')&(det.subset=='Fibro_only')][['site','donor_record','cells','library_UMIs']].drop_duplicates().to_string(index=False))
