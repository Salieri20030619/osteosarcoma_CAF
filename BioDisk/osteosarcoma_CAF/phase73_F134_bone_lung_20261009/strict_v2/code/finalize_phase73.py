from pathlib import Path
import os
os.environ['MPLBACKEND']='Agg'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,datetime,itertools,math,shutil,zipfile,traceback
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');ROOT=P/'phase73_F134_bone_lung_20261009';O=ROOT/'strict_v2'
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
def dump(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8388608),b''):h.update(b)
 return h.hexdigest()
def bhcheck(d):
 maxerr=0.;families=0
 for _,z in d.groupby('family'):
  p=z.P.fillna(1).to_numpy();ix=np.argsort(p);v=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[ix]=np.minimum(v,1.)
  valid=z.P.notna().to_numpy()
  if valid.any():maxerr=max(maxerr,float(np.max(np.abs(q[valid]-z.FDR.to_numpy()[valid]))))
  assert z.loc[~z.P.notna(),'FDR'].isna().all();families+=1
 assert maxerr<1e-12
 return families,maxerr
def main():
 assert (O/'core.exit_code').read_text().strip()=='0' and (O/'surface.exit_code').read_text().strip()=='0'
 qa=[];rows=[];primary_summary=[]
 files=sorted((O/'tables').glob('*__all_gene_statistics.tsv.gz'))
 for f in files:
  d=read(f);nf,err=bhcheck(d);qa.append(dict(file=f.name,families=nf,max_FDR_error=err))
  for key,g in d.groupby('context_id',sort=False):
   valid=g[g.P.notna()]
   if valid.empty:continue
   z=np.load(O/'objects'/(key+'.npz'));features=z['features'].astype(str);idx={v:i for i,v in enumerate(features)}
   choices=list(dict.fromkeys([valid.sort_values('P').feature.iloc[0],valid.iloc[len(valid)//2].feature,valid.iloc[-1].feature]+(['ENSG00000147257'] if 'ENSG00000147257' in set(valid.feature) else [])))
   for fe in choices:
    j=idx[fe];x=z['bone_expression'][:,j];y=z['lung_expression'][:,j];paired=g.design.iloc[0]=='paired'
    if paired:
     mask=np.isfinite(x)&np.isfinite(y);delta=y[mask]-x[mask];effect=float(delta.mean());s=np.array(list(itertools.product([-1,1],repeat=len(delta))));p=float(np.mean(np.abs(s@delta)>=abs(delta.sum())-len(delta)*1e-12))
    else:
     x=x[np.isfinite(x)];y=y[np.isfinite(y)];effect=float(y.mean()-x.mean());n=len(x);k=len(y);allv=np.r_[x,y]
     assert math.comb(n+k,k)<=65536
     total=float(allv.sum());observed=abs(effect);hits=0;den=0
     for comb in itertools.combinations(range(n+k),k):
      sy=float(allv[list(comb)].sum());null=abs(sy/k-(total-sy)/n);hits+=null>=observed-1e-12;den+=1
     p=hits/den
    rec=g[g.feature.eq(fe)].iloc[0];assert abs(effect-rec.mean_lung_minus_bone)<1e-10 and abs(p-rec.P)<1e-12,(key,fe,p,rec.P)
    rows.append(dict(context_id=key,feature=fe,P_recorded=rec.P,P_recomputed=p,effect_error=abs(effect-rec.mean_lung_minus_bone)))
   if g.sensitivity.iloc[0]=='all_cells' and g.level.iloc[0]=='main':primary_summary.append(dict(layer=g.layer.iloc[0],dataset=g.dataset.iloc[0],target=g.target.iloc[0],design=g.design.iloc[0],bone_N=int(g.bone_N.max()),lung_N=int(g.lung_N.max()),valid_genes=len(valid),nominal_up=int((valid.P.lt(.05)&valid.mean_lung_minus_bone.gt(0)).sum()),nominal_down=int((valid.P.lt(.05)&valid.mean_lung_minus_bone.lt(0)).sum()),FDR_significant=int(valid.FDR.lt(.05).sum())))
 for f in ['all_program_statistics.tsv.gz','all_proportion_statistics.tsv','all_ORA_statistics.tsv.gz','surface_background_statistics.tsv']:
  d=read(O/'tables'/f);nf,err=bhcheck(d);qa.append(dict(file=f,families=nf,max_FDR_error=err))
 save(pd.DataFrame(rows),'independent_P_effect_QA.tsv');save(pd.DataFrame(qa),'independent_FDR_QA.tsv');save(pd.DataFrame(primary_summary),'main_class_results_summary.tsv')
 # Direct reproduction of the reference PDF's supported F3 GPC3 result.
 old=read(P/'phase33_client_report_20260926/tables/strict_lung_site_all_classes_surface.tsv')
 old=old[(old.layer=='supported')&(old.target=='F3')&(old.symbol=='GPC3')&(~old.paired)].iloc[0]
 new=read(O/'tables/supported__SCPCP000017__all_gene_statistics.tsv.gz')
 new=new[(new.sensitivity=='all_cells')&(new.target=='F3')&(new.symbol=='GPC3')&(new.design=='independent')].iloc[0]
 assert new.bone_N==old.n_bone==6 and new.lung_N==old.n_lung==2
 assert abs(new.P-old.P)<1e-12 and abs(new.median_lung_minus_bone-old.median_effect)<1e-10
 reference_check=dict(reference='20260928 PDF page9, supported F3 GPC3',bone_N=6,lung_N=2,P=new.P,median_effect=new.median_lung_minus_bone,original_P_reproduced=True,new_FDR_family_documented=True)
 dump(O/'QA/REFERENCE_REPRODUCTION.json',reference_check)
 # Leave-one-patient-out effect stability is descriptive and creates no new P values.
 markers=read(O/'tables/surface_marker_evidence_review.tsv');lopo=[]
 for r in markers.to_dict('records'):
  z=np.load(O/'objects'/(r['context_id']+'.npz'));j=list(z['features'].astype(str)).index(r['feature']);x=z['bone_expression'][:,j];y=z['lung_expression'][:,j]
  vals=[]
  for side,vec in [('Bone',x),('Lung',y)]:
   for i in range(len(vec)):
    xx=np.delete(x,i) if side=='Bone' else x;yy=np.delete(y,i) if side=='Lung' else y
    if np.isfinite(xx).any() and np.isfinite(yy).any():vals.append(float(np.nanmean(yy)-np.nanmean(xx)))
  r.update(LOPO_effect_min=min(vals),LOPO_effect_max=max(vals),LOPO_direction_consistent=all(np.sign(v)==np.sign(r['mean_lung_minus_bone']) for v in vals))
  lopo.append(r)
 markers=pd.DataFrame(lopo);save(markers,'surface_marker_evidence_review_with_LOPO.tsv')
 # Human-flow product evidence; no claim that these products identify CAF specifically.
 products=[
  dict(gene='GPC1',product='MAB45192',clone='1019718',evidence='Human A549 membrane-protein flow demonstration; unconjugated antibody with secondary detection',url='https://www.rndsystems.com/products/human-glypican-1-antibody-1019718_mab45192',role='F1 exploratory observation'),
  dict(gene='GPC3',product='FAB2119A',clone='307801',evidence='Human HepG2 membrane-associated flow demonstration; APC conjugate',url='https://www.rndsystems.com/products/human-glypican-3-apc-conjugated-antibody-307801_fab2119a',role='Supported-F3 exploratory observation; retains original report candidate'),
  dict(gene='SDC3',product='FAB3539A',clone='polyclonal',evidence='Human A172 flow demonstration; APC conjugate',url='https://resources.rndsystems.com/pdfs/datasheets/fab3539a.pdf?v=20240325',role='F3/F4 shared candidate; cannot distinguish these classes alone'),
  dict(gene='SUSD2',product='749522',clone='W5C5',evidence='Human flow qualified in manufacturer sheet',url='https://www.bdbiosciences.com/content/dam/bdb/product_assets/product_pdf/singcolorpureantibody/bb700/pdf_1/749522.pdf',role='Supported-F3 auxiliary candidate; low RNA detection'),
  dict(gene='FLRT2',product='MAB2877/FAB2877 series',clone='367205',evidence='Reviewed manufacturer catalog lists Western blot; human flow suitability not established by these pages',url='https://www.rndsystems.com/target/flrt2/antibodies',role='F4 RNA candidate; do not assume conjugation implies flow validation'),
 ]
 save(pd.DataFrame(products),'flow_reagent_evidence_review.tsv')
 # Figures use every eligible patient; displays do not create selection thresholds.
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':170})
 prop=read(O/'tables/patient_site_proportions.tsv.gz');stats=read(O/'tables/all_proportion_statistics.tsv')
 fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
 for row,layer in enumerate(['supported','expanded']):
  for col,target in enumerate(['F1','F3','F4']):
   ax=axes[row,col];z=prop[(prop.sensitivity=='all_cells')&(prop.layer==layer)&(prop.dataset=='SCPCP000017')&(prop.target==target)&(prop.denominator=='layer')]
   for k,site in enumerate(['Bone','Lung']):
    v=z[z.site==site].fraction.dropna().to_numpy()*100;xx=k+np.linspace(-.07,.07,len(v)) if len(v) else []
    ax.scatter(xx,v,c=['#3974a8','#bd6547'][k],s=28,alpha=.8);ax.plot([k-.15,k+.15],[np.mean(v),np.mean(v)],c='black',lw=1) if len(v) else None
   ax.set_xticks([0,1],['Bone','Lung']);ax.set_title(layer+' '+target);ax.set_ylabel('% within frozen layer')
   if layer=='supported' and target=='F4':ax.text(.5,.6,'No F4 class in reference layer',ha='center',transform=ax.transAxes,fontsize=9)
 fig.suptitle('F1/F3/F4 patient-level captured proportions | SCPCP000017\nEach point is one patient-site; paired and independent tests are reported separately')
 fig.savefig(O/'figures/01_patient_proportions.png');plt.close(fig)
 chosen=[('expanded','F1','GPC1'),('supported','F3','GPC3'),('expanded','F4','SDC3'),('expanded','F1','CD74')]
 fig,axes=plt.subplots(1,4,figsize=(14,4.8),constrained_layout=True)
 for ax,(layer,target,gene) in zip(axes,chosen):
  r=markers[(markers.layer==layer)&(markers.target==target)&(markers.symbol==gene)].iloc[0];z=np.load(O/'objects'/(r.context_id+'.npz'));j=list(z['features'].astype(str)).index(r.feature)
  for k,site in enumerate(['bone','lung']):
   v=z[site+'_expression'][:,j];ax.scatter(k+np.linspace(-.08,.08,len(v)),v,s=35,c=['#3974a8','#bd6547'][k]);ax.plot([k-.15,k+.15],[np.nanmean(v)]*2,c='black')
  ax.set_xticks([0,1],['Bone','Lung']);ax.set_title(f'{layer} {target}\n{gene} | P={r.P:.4g}; FDR={r.FDR:.2g}');ax.set_ylabel('Patient pseudobulk log2(CP10k + 0.1)')
 fig.suptitle('Illustrative site-associated surface candidates | independent patient groups')
 fig.savefig(O/'figures/02_surface_patient_expression.png');plt.close(fig)
 bg=read(O/'tables/surface_background_statistics.tsv');bg=bg[bg.site.eq('Lung')]
 fig,ax=plt.subplots(figsize=(12,6.5),constrained_layout=True);labels=[];matrix=[]
 backgrounds=sorted(bg.background.unique())
 for layer,target,gene in chosen:
  r=bg[(bg.layer==layer)&(bg.target==target)&(bg.symbol==gene)].set_index('background');matrix.append(r.mean_target_minus_background.reindex(backgrounds).to_numpy());labels.append(f'{layer} {target} | {gene}')
 arr=np.array(matrix);lim=max(1,float(np.nanmax(np.abs(arr))));im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim)
 ax.set_yticks(range(len(labels)),labels);ax.set_xticks(range(len(backgrounds)),backgrounds,rotation=50,ha='right');fig.colorbar(im,ax=ax,label='Target - background\nlog2(CP10k + 0.1)')
 ax.set_title('Lung-site RNA background screen | descriptive effect, not measured flow purity')
 fig.savefig(O/'figures/03_marker_background_effects.png');plt.close(fig)
 programs=read(O/'tables/nominal_programs.tsv');cust=programs[(programs.sensitivity=='all_cells')&programs.program.str.startswith('customer|')]
 genes=pd.DataFrame(primary_summary)
 prose='''F1/F3/F4：保留原分类的骨—肺两层比较（2026-10-09）

本轮结论
原报告的支持层1224、扩展层8174及F1–F5标签完整保留。本轮比较F1/F3/F4及已有固定内部状态的骨—肺差异，分别分析患者级比例和群内表达/功能。配对与非配对患者分开，来源分开，不限于旧B1病史标签。

一、比例比较未发现名义显著变化
主类占全部捕获细胞、占对应候选层，以及内部状态占父类三种比例均已计算。当前预设比较没有原始P<0.05的条目。零分子保留，零分母为不可评估。这里是捕获细胞比例，不是组织中的绝对丰度。

二、表达差异主要位于扩展F1、支持F3、扩展F4
主研究SCPCP000017中，独立患者比较：扩展F1骨8/肺3；支持F3骨6/肺2；扩展F4骨11/肺2。完整名义上下调数量见main_class_results_summary.tsv。当前基因结果均未达到FDR<0.05；按原始P<0.05保留探索性候选，不把FDR阴性写成无信号。扩展F3仍只有骨2/肺1，无有效推断性检验；支持F3与扩展F3的成员不同，不能互相替代。
支持F3的GPC3结果复现原报告：中位/平均差+3.1544，骨6/肺2，P=0.035714、FDR=1。扩展F1的apCAF程序及扩展F4的iCAF程序呈名义升高，不能仅据程序命名确认新的CAF类型或免疫抑制功能。F4_C9内部状态的上调基因有基质相关富集线索；基因集富集模型的P/FDR与患者级程序检验分开保存，不能互相替代。

三、表面候选与检测可行性
全部名义表面候选共85条亚群/基因记录、69个唯一基因，全部回到同来源原始全细胞背景核查。完整结果见surface_marker_evidence_review_with_LOPO.tsv；未用自拟纯度或表达比例阈值删减候选。
F1的GPC1、支持F3的GPC3，以及F3/F4的SDC3可作为优先观察线索。这些基因有厂家人源流式应用资料，但仍没有本项目蛋白验证。GPC3延续原报告观察方向；SDC3在F3和F4均有信号，不能单独区分二者。F1的CD74/CSF1R/IFNGR2虽有名义差异，但存在免疫或其他谱系背景表达；不把它们直接作为CAF特异分选标志。FLRT2保留RNA线索，当前核对的厂家产品资料没有确立人体表面流式适用性。产品来源见flow_reagent_evidence_review.tsv。尚未形成经过验证的F1/F3/F4分选组合，未完成的蛋白/抗体性能不冒称已验证。

四、关键解释边界
扩展F4独立比较的肺侧仅2名患者，其中一名只有1个F4细胞；扩展F1独立比较肺侧3名患者，其中两名分别仅2个、1个F1细胞。支持F3肺侧两名患者分别有6个、114个细胞。患者数和细胞覆盖均应随候选解释。去双细胞敏感性和逐患者剔除的描述性效应范围已保存；这些不增加独立样本量。骨—肺取材关联不等于肺转移因果效应，也不直接证明复发、基质物理包绕或免疫排斥。

五、原报告及数据完整性
原报告第4页分类数量全部复现；第9页支持F3 GPC3的患者数、效应和P值复现。支持层没有独立F4，未人为补出。主类标签未改变；内部状态复用已有冻结映射，未重新调聚类参数。
对照原报告时发现旧派生表把Paraspinal归为Bone。正式strict_v2改用已验收的clinical_site_review口径，将椎旁样本留作Unknown_or_mixed；此前A/B1/B2标签未变。首轮计算完整保留为superseded，仅strict_v2为本轮正式输出。表面筛选使用完整原始膜蛋白注释及原报告候选清单，未局限于旧108个主类候选基因。

方法与文件
单基因：患者×部位原始计数汇总，以完整文库计算log2(CP10k+0.1)，效应为肺减骨。功能程序：冻结基因集内的平均log1p(CP10k)，按患者汇总。比例：同患者同部位合并分子分母。检验沿用报告中的精确患者置换/符号翻转及预先登记的Monte Carlo规则；缺重复或无变化不产生P值。FDR按基因、功能程序、比例、富集、背景比较分别在预定义家族内计算。所有可评估来源及不可评估原因均保留。
主要文件：frozen_target_registry.tsv；expression_comparison_registry.tsv；all_proportion_statistics.tsv；各来源all_gene_statistics.tsv.gz；all_program_statistics.tsv.gz；all_ORA_statistics.tsv.gz；surface_marker_evidence_review_with_LOPO.tsv；patient_site_proportions.tsv.gz；expression_patient_roster.tsv。
本轮完成条件：退出码、原始计数核查、独立P/FDR复核、报告结果复现、marker全背景筛查与图件检查。完整数据和中间证据随压缩包保存。
'''
 (O/'delivery/RESULTS_ZH.txt').write_text(prose,encoding='utf-8')
 (O/'delivery/RESULTS_ZH.txt').write_text(prose+'\n主类结果摘要\n'+genes.to_string(index=False)+'\n\n名义客户程序\n'+cust[['layer','target','bone_N','lung_N','program','mean_lung_minus_bone','P','FDR']].to_string(index=False)+'\n',encoding='utf-8')
 dump(O/'QA/FINAL_NUMERIC_ACCEPTANCE.json',dict(status='PASS',P_and_effect_checks=len(rows),FDR_files=len(qa),max_FDR_error=max(q['max_FDR_error'] for q in qa),reference_GPC3_reproduced=True,figures_created=3,visual_review_pending=True,full_marker_background_completed=True))
 print('FINAL_NUMERIC_ACCEPTANCE',len(rows),'P_EFFECT_CHECKS',len(qa),'FDR_FILES',flush=True)
if __name__=='__main__':
 try:main()
 except BaseException:
  traceback.print_exc();(O/'finalize.exit_code').write_text('1\n');raise
 else:(O/'finalize.exit_code').write_text('0\n')
