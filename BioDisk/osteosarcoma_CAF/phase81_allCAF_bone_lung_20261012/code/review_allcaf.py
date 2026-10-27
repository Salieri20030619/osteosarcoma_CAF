from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,itertools,math,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase81_allCAF_bone_lung_20261012'
def read(name):return pd.read_csv(O/'tables'/name,sep='\t',low_memory=False)
def save(d,name):d.to_csv(O/'tables'/name,sep='\t',index=False)
def bh(p):
    p=np.array(p,float);ii=np.argsort(p);s=np.minimum.accumulate((p[ii]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[ii]=np.minimum(s,1);return q
def exact(x,y,paired):
    if paired:
        d=(y-x)[np.isfinite(x)&np.isfinite(y)]
        if len(d)<2 or not np.any(d):return np.nan
        obs=abs(d.mean());null=[abs(np.mean(d*s)) for s in itertools.product([-1,1],repeat=len(d))]
    else:
        x=x[np.isfinite(x)];y=y[np.isfinite(y)]
        if min(len(x),len(y))<2:return np.nan
        a=np.r_[x,y]
        if np.all(a==a[0]):return np.nan
        obs=abs(y.mean()-x.mean());null=[]
        for ix in itertools.combinations(range(len(a)),len(y)):
            mask=np.zeros(len(a),bool);mask[list(ix)]=True;null.append(abs(a[mask].mean()-a[~mask].mean()))
    return np.mean(np.array(null)>=obs-1e-12)
assert (O/'logs/analysis.exit_code').read_text().strip()=='0'
receipt=json.loads((O/'COMPUTATION_COMPLETE.json').read_text())
assert receipt['original_inputs_unchanged']
genes=read('all_genes_statistics.tsv.gz');programs=read('all_programs_statistics.tsv.gz');proportions=read('all_proportions_statistics.tsv.gz');roster=read('patient_expression_roster.tsv');summ=read('comparison_results_summary.tsv');reg=read('comparison_registry.tsv')
bh_checks=[]
for df,kind in [(genes,'genes'),(programs,'programs'),(proportions,'proportions')]:
    for fam,g in df.groupby('family'):
        z=g[g.P.notna()]
        if len(z):bh_checks.append(dict(kind=kind,family=fam,valid_tests=len(z),error=float(np.max(np.abs(bh(z.P)-z.FDR.to_numpy())))))
pd.DataFrame(bh_checks).to_csv(O/'QA/independent_BH_checks.tsv',sep='\t',index=False)
probes=[];alignment=[]
for row in summ[summ.sensitivity.eq('all_cells')].itertuples():
    with np.load(O/'objects'/f'{row.context_id}.npz',allow_pickle=True) as z:
        x=z['bone_expression'];y=z['lung_expression'];features=z['features'].astype(str);b=z['bone_patients'].astype(str);l=z['lung_patients'].astype(str)
        xprog=z['bone_programs'];yprog=z['lung_programs'];terms=z['programs'].astype(str)
    alignment.append(dict(context_id=row.context_id,disjoint_or_paired=(set(b)==set(l) if row.design=='paired' else not set(b)&set(l)),bone_N_matches=len(b)==row.bone_N,lung_N_matches=len(l)==row.lung_N))
    combos=2**len(b) if row.design=='paired' else math.comb(len(b)+len(l),len(l))
    if combos>65536:continue
    for df,kind,col,features_,bx,ly in [(genes,'gene','feature',features,x,y),(programs,'program','program',terms,xprog,yprog)]:
        g=df[df.context_id.eq(row.context_id)].set_index(col)
        chosen=list(g[g.P.notna()].sort_values('P').head(3).index)
        if kind=='gene':chosen+=list(features[::max(1,len(features)//5)])
        else:chosen+=list(g.index[g.index.str.startswith('customer|')])
        for f in dict.fromkeys(chosen):
            j=np.flatnonzero(features_==f)[0];p=exact(bx[:,j],ly[:,j],row.design=='paired');old=g.loc[f,'P']
            probes.append(dict(context_id=row.context_id,kind=kind,feature=f,stored_P=old,independent_P=p,error=abs(p-old) if np.isfinite(p) and pd.notna(old) else np.nan,NA_matches=np.isnan(p)==pd.isna(old)))
pd.DataFrame(probes).to_csv(O/'QA/independent_P_checks.tsv',sep='\t',index=False);pd.DataFrame(alignment).to_csv(O/'QA/patient_alignment.tsv',sep='\t',index=False)
assert max(r['error'] for r in bh_checks)<1e-9
assert pd.DataFrame(probes).error.max()<1e-9 and pd.DataFrame(probes).NA_matches.all()
assert all(r['disjoint_or_paired'] and r['bone_N_matches'] and r['lung_N_matches'] for r in alignment)
primary=summ[summ.sensitivity.eq('all_cells')].copy();primary=primary.sort_values(['layer','dataset','design'])
testable=primary[primary.status.eq('testable')].copy()
cov=[]
for r in testable.itertuples():
    for site,ids in [('Bone',r.bone_patients),('Lung',r.lung_patients)]:
        for donor in ids.split(';'):
            z=roster[(roster.sensitivity=='all_cells')&(roster.layer==r.layer)&(roster.dataset==r.dataset)&(roster.donor_record==donor)&(roster.site==site)]
            assert len(z)==1
            cov.append(dict(context_id=r.context_id,layer=r.layer,dataset=r.dataset,design=r.design,site=site,donor_record=donor,cells=int(z.cells.iloc[0]),library_UMIs=float(z.library_UMIs.iloc[0])))
cov=pd.DataFrame(cov);save(cov,'main_patient_coverage.tsv')
customer=programs[(programs.sensitivity=='all_cells')&programs.program.str.startswith('customer|')].copy();save(customer,'all_customer_program_results.tsv')
main_g=genes[(genes.sensitivity=='all_cells')&genes.nominal];main_p=programs[(programs.sensitivity=='all_cells')&programs.nominal]
top=main_g.sort_values(['context_id','P','feature']).groupby('context_id').head(12);save(top,'display_top12_genes_per_context.tsv')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
labels=[f"{'S' if r.layer=='supported' else 'E'} | {r.dataset} | {r.design}" for r in testable.itertuples()];yy=np.arange(len(testable))
fig,ax=plt.subplots(1,2,figsize=(14,max(5,.6*len(testable)+2)))
ax[0].barh(yy-.17,testable.bone_N,.32,color='#326b9b',label='Bone');ax[0].barh(yy+.17,testable.lung_N,.32,color='#ce6840',label='Lung');ax[0].set_yticks(yy,labels);ax[0].invert_yaxis();ax[0].set_xlabel('Biological patients');ax[0].legend(frameon=False)
ax[1].barh(yy,testable.nominal_up,color='#ce6840',label='Higher in lung');ax[1].barh(yy,-testable.nominal_down,color='#326b9b',label='Lower in lung');ax[1].set_yticks(yy,[]);ax[1].invert_yaxis();ax[1].axvline(0,color='gray',lw=.8);ax[1].set_xlabel('Nominal gene records (raw P < 0.05)');ax[1].legend(frameon=False)
for i,r in enumerate(testable.itertuples()):ax[1].text(r.nominal_up+max(10,testable.nominal_up.max()*.02),i,f'{r.nominal_up}/{r.nominal_down}; FDR hits: {r.FDR_genes}',va='center',fontsize=8)
ax[1].set_xlim(-max(50,testable.nominal_down.max()*1.15),max(100,testable.nominal_up.max()*1.5))
fig.suptitle('All CAF candidates | anatomical bone versus lung',fontsize=16)
fig.text(.015,.015,'S: supported; E: expanded. All F1-F5 pooled within each layer; layers overlap. Paired N represents the same donors at both sites.',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(O/'figures/01_allCAF_comparisons.png',dpi=180);plt.close(fig)
cp=customer[customer.context_id.isin(testable.context_id)]
order=testable.context_id.tolist();terms=sorted(cp.program.unique());matrix=cp.pivot(index='context_id',columns='program',values='mean_lung_minus_bone').reindex(index=order,columns=terms)
fig,ax=plt.subplots(figsize=(max(9,len(terms)*1.4),max(5,len(order)*.55+2)))
v=np.nanmax(np.abs(matrix.to_numpy()));v=max(v,.01);im=ax.imshow(matrix.to_numpy(),cmap='RdBu_r',vmin=-v,vmax=v,aspect='auto')
for i,ctx in enumerate(order):
    for j,term in enumerate(terms):
        r=cp[(cp.context_id==ctx)&(cp.program==term)].iloc[0]
        mark='**' if r.FDR<.05 else ('*' if r.P<.05 else '')
        ax.text(j,i,f'{r.mean_lung_minus_bone:+.2f}{mark}',ha='center',va='center',fontsize=9,color='white' if abs(r.mean_lung_minus_bone)>v*.55 else '#222')
ax.set_yticks(range(len(order)),labels);ax.set_xticks(range(len(terms)),[t.split('|')[1].replace('_customer','').replace('_',' ') for t in terms],rotation=35,ha='right');ax.set_title('Patient-level customer programs | lung minus bone',pad=15);fig.colorbar(im,ax=ax,label='Mean log1p(CP10k) score difference',shrink=.75)
fig.text(.015,.015,'* raw P < 0.05; ** BH FDR < 0.05. No mark means no nominal evidence; all exact P/FDR and effective N are retained in tables.',fontsize=9)
fig.tight_layout(rect=[0,.055,1,.98]);fig.savefig(O/'figures/02_allCAF_customer_programs.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(12,max(5,.55*len(order)+1)))
for i,ctx in enumerate(order):
    for j,(side,color) in enumerate([('Bone','#326b9b'),('Lung','#ce6840')]):
        z=cov[(cov.context_id==ctx)&(cov.site==side)];offset=np.linspace(-.07,.07,len(z)) if len(z)>1 else [0]
        ax.scatter(z.cells,np.full(len(z),i)+(j-.5)*.23+offset,c=color,s=28,alpha=.8,label=side if i==0 else None)
ax.set_xscale('log');ax.set_yticks(range(len(order)),labels);ax.invert_yaxis();ax.set_xlabel('CAF candidate cells per patient at each anatomical site');ax.set_title('Patient coverage after pooling all CAF classes');ax.legend(frameon=False);fig.tight_layout();fig.savefig(O/'figures/03_allCAF_patient_coverage.png',dpi=180);plt.close(fig)

summarycols=['layer','dataset','design','bone_N','lung_N','nominal_up','nominal_down','FDR_genes','nominal_programs','FDR_programs']
signif_customer=customer[customer.P.lt(.05)]
report='''全部CAF候选细胞的骨病灶与肺病灶比较
2026年10月12日

本轮范围
按先前确认的方案一，以病灶实际解剖取材部位划分骨组和肺组，不以是否转移作为入组条件。本轮将每层F1-F5全部合并后重新计算，不再先分大类或内部群。支持层1224和扩展层8174分别分析，细胞成员及标签不变。扩展层包含身份未解决或成骨样候选，不能等同于已确认纯CAF。

方法
按来源分别，以患者为生物学重复。同患者同部位汇总原始计数，完整文库归一化后计算log2(CP10k+0.1)；程序使用患者内细胞的平均log1p(CP10k)及冻结面板。两部位均有候选细胞者作配对，其余作独立组比较。每侧不足2名独立患者或配对少于2名时不计算推断P。效应均为肺减骨。
采用已审计的精确患者置换/符号翻转，超过预登记枚举规模时使用固定种子的Monte Carlo。按来源、候选层、设计、敏感性和假设类型分别对有效P作BH；无效检验保留NA且不进入BH，旧NA按1计数的口径另存为诊断列。本规则在新分析前冻结。名义候选仍按原始P<0.05选择，不要求FDR显著。比例分母为完整来源捕获细胞。去双细胞结果、F类组成及全部不可评估比较独立保存。

整体CAF比较结果
'''
report+=primary[summarycols].to_string(index=False)+'\n'
report+=f'\n全细胞口径共{len(testable)}项表达比较具备检验条件。名义基因记录{len(main_g)}条，其中FDR<0.05为{int(main_g.FDR_significant.sum())}条；名义患者级程序{len(main_p)}条，其中FDR<0.05为{int(main_p.FDR_significant.sum())}条。数字按来源、层和设计计记录，不能作为独立发现数。\n'
report+='\n客户预设程序的名义信号\n'
if len(signif_customer):report+=signif_customer[['layer','dataset','design','bone_N','lung_N','program','mean_lung_minus_bone','P','FDR']].to_string(index=False)+'\n'
else:report+='当前数据没有原始P<0.05的客户预设程序；完整方向和效应仍保留。\n'
prop_primary=proportions[proportions.sensitivity.eq('all_cells')]
report+='\n整体候选比例\n'+prop_primary[prop_primary.P.notna()][['layer','dataset','design','bone_N','lung_N','bone_mean_fraction','lung_mean_fraction','mean_lung_minus_bone','P','FDR']].to_string(index=False)+'\n'
report+='\n基因与程序的完整证据\n所有基因特征保留完整P/FDR、效应和有效患者数；名义结果关联去双细胞效应与对应配对结果。top12表仅用于阅读展示，不额外筛除候选。ORA为基因列表富集模型，其P不是患者功能效应P。\n'
report+='\n解释边界\n合并所有F类增加每名患者的细胞覆盖，也使表达差异同时包含候选群体组成和群内表达变化；不能把整体差异直接归给某个CAF亚型。来源分别分析、两层互相重叠，支持层和扩展层不能算独立验证。转移和治疗史留作注释，解剖部位关联不等于转移因果作用。部分来源仍缺有效患者重复。\n'
report+='\n与上轮的关系\n上轮为F1/F3/F4及内部群的比较，本轮为整个CAF候选池，两者估计对象不同，不能只按显著条目数量评价哪个结果更正确。FDR有效检验家族的口径变化已明确记录，配套legacy_NA_as_one_FDR_diagnostic列可以辨别其影响。既有Phase80结果完整保留。\n'
report+='\n附件\nall_genes_statistics、all_programs_statistics、all_proportions_statistics为全量；main_nominal_*_with_sensitivity为名义结果及配对/去双细胞证据；patient_Fclass_composition为组成；comparison_registry为全部可评估及不可评估登记。三张图展示患者数、名义基因、客户程序及患者细胞覆盖。\n'
(O/'delivery/全部CAF候选_骨肺病灶比较_20261012.txt').write_text(report,encoding='utf-8')
qa=dict(status='numeric_QA_passed_visual_review_pending',independent_gene_program_P_probes=len(probes),max_P_error=float(pd.DataFrame(probes).error.max()),BH_families_checked=len(bh_checks),max_BH_error=max(r['error'] for r in bh_checks),patient_alignment_verified=True,source_inputs_unchanged=receipt['original_inputs_unchanged'],figures=3,visual_review_pending=True)
(O/'QA/NUMERIC_ACCEPTANCE.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps(qa))
print(primary[summarycols].to_string(index=False))
print('CUSTOMER_SIGNALS\n'+signif_customer[['layer','dataset','design','bone_N','lung_N','program','mean_lung_minus_bone','P','FDR']].to_string(index=False))
