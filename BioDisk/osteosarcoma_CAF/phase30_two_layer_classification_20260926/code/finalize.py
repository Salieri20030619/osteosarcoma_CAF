from pathlib import Path
import json,datetime,hashlib,shutil
import numpy as np,pandas as pd,anndata as ad
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase30_two_layer_classification_20260926'
now=datetime.datetime.now().astimezone().isoformat()
assert (O/'classification.exit_code').read_text().strip()=='0' and (O/'CLASSIFICATION_COMPLETE').exists()
classes=pd.read_csv(O/'tables/class_profiles.tsv',sep='\t')
anno=pd.read_csv(O/'tables/class_annotation_frozen.tsv',sep='\t')
source=pd.read_csv(O/'tables/layer_membership_all11379.tsv.gz',sep='\t',low_memory=False)
assert len(source)==11379 and source.cell_id.is_unique
s=pd.read_csv(O/'tables/supported_classified_cells.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
e=pd.read_csv(O/'tables/expanded_classified_cells.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
assert len(s)==1224 and len(e)==8174 and set(s.index)<=set(e.index)
assert s.dual_reference.sum()==e.dual_reference.sum()==114
assert s.class_id.notna().all() and e.class_id.notna().all()
assert source.in_extended.sum()==8174 and source.in_supported.sum()==1224
st=pd.read_csv(O/'tables/class_program_statistics.tsv',sep='\t')
sf=pd.read_csv(O/'tables/class_surface_gene_statistics.tsv',sep='\t')
for file in [O/'tables/class_program_statistics.tsv']+list((O/'tables').glob('*_all_gene_statistics.tsv.gz')):
 d=pd.read_csv(file,sep='\t')
 assert d.loc[d.n_independent_donors<2,'P'].isna().all()
 test=d.P.notna();assert d.loc[test,'P'].between(0,1).all()
 assert (d.loc[test,'P']>=np.exp2(1-d.loc[test,'n_independent_donors'])-1e-12).all()
 for _,g in d.groupby('family'):
  calc=multipletests(g.P.fillna(1),method='fdr_bh')[1];valid=g.P.notna()
  assert np.allclose(g.loc[valid,'FDR'],calc[valid])
 assert np.array_equal(d.nominal_exploratory.to_numpy(),(test&(d.P<.05)&(d.median_delta>0)).to_numpy())
summary=pd.read_csv(O/'tables/statistics_summary.tsv',sep='\t')
nominal_surface=sf[sf.nominal_exploratory].copy()
nominal_surface['priority_basis']=np.where(nominal_surface.symbol.isin(['PDGFRB','MCAM'])&nominal_surface.class_id.eq('F2'),'same-direction nominal support in both nested layers','one-layer nominal lead')
nominal_surface.to_csv(O/'delivery/流式候选名义显著证据.tsv',sep='\t',index=False)
pivot=classes.pivot(index=['class_id','class_name'],columns='layer',values='cells').fillna(0).astype(int)
pivot.to_csv(O/'delivery/新版双层分类数量.tsv',sep='\t')
cross=pd.read_csv(O/'tables/class_correspondence_counts.tsv',sep='\t',index_col=0)
agreement=sum(cross.loc[k,k] for k in cross.index if k in cross.columns)/cross.to_numpy().sum()
my_names={
'supported_myCAF':{0:'COL4/NOTCH3基底膜-收缩状态',1:'MYH11/ACTA2收缩状态',2:'POSTN/FN1基质活化状态',3:'MGP/翻译应激状态',4:'极小混合群（不单独建立亚型）'},
'expanded_myCAF':{0:'COL1/COL11基质沉积状态',1:'ACTA2/MYH11收缩状态',2:'AQP1/IGFBP5重塑状态',3:'POSTN/VCAN基质活化状态',4:'翻译/增殖相关状态',5:'MGP/热休克状态'}}
myrows=[]
for layer,mapping in my_names.items():
 pr=pd.read_csv(O/'tables'/f'{layer}_cluster_profiles.tsv',sep='\t')
 for r in pr.to_dict('records'):myrows.append(dict(layer=layer,cluster=r['cluster'],cells=r['cells'],donor_records=r['donor_records'],expression_state=mapping[int(r['cluster'])],validated_subtype=False,no_main_class_P_transferred=True))
pd.DataFrame(myrows).to_csv(O/'tables/myCAF_internal_state_review.tsv',sep='\t',index=False)
colors={'F1':'#3274a1','F2':'#e1812c','F3':'#3a923a','F4':'#c03d3e','F5':'#777777'}
english={'F1':'Inflammatory / adventitial-like','F2':'Contractile / active matrix','F3':'Matrix remodeling','F4':'Osteogenic / matrix overlap','F5':'Mixed / unresolved'}
fig,axes=plt.subplots(1,2,figsize=(16,6),layout='constrained')
for ax,layer,obs in zip(axes,['supported','expanded'],[s,e]):
 a=ad.read_h5ad(O/'objects'/f'{layer}_integrated.h5ad',backed='r')
 for cl in colors:
  ix=np.flatnonzero(obs.loc[a.obs_names,'class_id'].eq(cl).to_numpy())
  if len(ix):ax.scatter(a.obsm['X_umap'][ix,0],a.obsm['X_umap'][ix,1],s=5,label=f'{cl}: {english[cl]} ({len(ix)})',color=colors[cl],rasterized=True)
 ax.set_title(f'{layer.capitalize()} layer: {len(obs):,} candidates');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(fontsize=7,frameon=True,loc='best');a.file.close()
fig.savefig(O/'figures/two_layer_classification.png',dpi=170);plt.close(fig)
# Preserve the pre-precision QA and reconcile its summary with the final exact statistics.
qf=O/'delivery/CLASSIFICATION_QA.json'
shutil.copy2(qf,O/'precision_audit_before_exact/CLASSIFICATION_QA.json')
q=json.loads(qf.read_text());q['statistics_pools']=summary.where(pd.notna(summary),None).to_dict('records');q['precision']='All existing tests exact, maximum15 donors';q['verified_at']=now
qf.write_text(json.dumps(q,ensure_ascii=False,indent=2,default=str))
qa={'at':now,'supported_cells':1224,'expanded_cells':8174,'dual_reference_cells':114,'layers_nested':True,'classified_all_research_cells':True,'no_additional_exclusion_of_expanded_6950':True,'class_mapping_sha256':hashlib.sha256((O/'tables/class_annotation_frozen.tsv').read_bytes()).hexdigest(),'all_patient_P_bounds_and_FDR_families_verified':True,'all_existing_tests_exact':True,'shared_cell_class_agreement':float(agreement),'shared_cell_cluster_ARI':0.6493967154973032,'main_UMAP_and_both_myCAF_UMAP_visual_QA':'PASS','classification_summary_plot':'rendered; visual inspection pending','old_analysis_preserved':True,'114_cells_separately_clustered':False}
(O/'delivery/FINAL_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
# Read exact final values directly into the report.
def value(layer,cls,feature,field,kind='program'):
 d=st if kind=='program' else sf
 col='feature' if kind=='program' else 'symbol'
 return d[(d.layer==layer)&(d.pool=='SCPCP000017')&(d.class_id==cls)&(d[col]==feature)].iloc[0][field]
mdtable='| 新主类 | 支持候选层 | 扩展探索层 |\\n|---|---:|---:|\\n'.replace('\\n',chr(10))
for (cid,name),r in pivot.iterrows():mdtable+=f"| {cid} {name} | {r.get('supported',0):,} | {r.get('expanded',0):,} |"+chr(10)
text=f"""# CAF双层研究与新版候选分类（2026-09-26）

## 已恢复支持层与扩展层并行研究

本轮按用户最新范围执行：1,224个支持候选和8,174个扩展探索候选为两个主要研究层；114个双方法支持细胞只作证据对照，未单独分群。扩展层包括此前保留的6,950个身份未决候选，未再因成骨、增殖或应激标签统一剔除。3,205个竞争谱系细胞及全部107份样本的原始证据继续保留。

两层分别覆盖46和59份样本、41和52个供者记录，均涉及9个来源。供者记录不直接作为跨来源合并的独立N。两层分别得到12和24个表达群，再根据完整marker与程序合并为四种表达表型及一个未定类别；命名与映射在查看新P值前冻结。

{mdtable}

F4在支持层没有单独形成该主类，不能解释为支持层没有成骨基因表达。F4的5,304个扩展候选继续作为研究对象，身份未决与成骨程序分别记录。F5是保留的混合/未定类别，不是新的细胞谱系。增殖、线粒体及应激状态为独立标签。

## 收缩-活化基质类提供跨层一致的流式候选线索

SCPCP000017中，F2相对同患者、同取材部位的其他候选：
- 支持层：10位独立患者，myCAF程序9/10同向升高，中位差+1.3653，P=0.00390625，FDR=0.125，属于名义显著。
- 扩展层：12位独立患者，myCAF程序11/12同向升高，中位差+0.6253，P=0.00146484375，FDR={value('expanded','F2','program_myCAF_customer','FDR'):.6g}；基质程序12/12升高，中位差+1.3517，P=0.00048828125，FDR={value('expanded','F2','program_mCAF_matrix','FDR'):.6g}。

程序差的单位为平均log1p CP10k分数。两层共享细胞，上述为分层一致性，不是独立队列验证。

F2的表面候选：
- PDGFRB：支持层10/10同向，P=0.001953125、FDR=1；扩展层11/12同向，P=0.0029296875、FDR={value('expanded','F2','PDGFRB','FDR','gene'):.6g}。
- MCAM：支持层9/10同向，P=0.00390625、FDR=1；扩展层10/12同向，P=0.0048828125、FDR={value('expanded','F2','MCAM','FDR','gene'):.6g}。
- THY1在扩展层、JAG1在支持层提供单层名义显著线索；完整效应与P/FDR见“流式候选名义显著证据.tsv”。

因此，PDGFRB和MCAM可优先进入小规模分析型流式的表达/共表达及身份验证，THY1、JAG1可作为次级比较候选。F2仍有血管周/周细胞程序重叠，不将表面阳性直接命名为已确认myCAF。RNA门控覆盖及背景混入已分别给出，不将其当作实测流式纯度。

## 基质重塑与成骨重叠类保留独立研究价值

支持层F3的基质程序在9位患者中8位同向，中位差+0.5656，P=0.0234375，FDR=0.2125，为探索性信号。

扩展层F4的成骨程序在15位患者中9位同向，中位差+0.1157，P={value('expanded','F4','program_Osteogenic_customer','P'):.8g}，FDR={value('expanded','F4','program_Osteogenic_customer','FDR'):.6g}，为名义显著；CD36提供候选表面线索。保留该方向，不据此确认这些细胞全部是CAF或其成骨程序具有转移因果作用。

F1的DCN、COL15A1、CFD、CXCL12、SFRP2及炎症亚结构支持描述性分类；未将所有主类都包装为统计显著类别。名义显著用于候选优先级，不用于反复调整群体直至获得显著。

## 两层对应和myCAF内部探索已完成

共同1,224个细胞的群体ARI=0.6494，主类一致率={agreement:.1%}。F1有466/478个支持层细胞在扩展层仍属F1；F2为67/90，F3为242/343。F3有86个细胞在扩展层落入F2，提示收缩与活化基质边界会随候选范围改变。

按新分类F2分别进行内部重新建模：支持层90个细胞得到5个计算群，扩展层225个得到6个计算群。它们来自当前层内定义，不以旧230个myCAF-like细胞作数量上限。支持层一个计算群只有3个细胞，保留记录但不单独建立生物学亚型；其余按MYH11/ACTA2收缩、POSTN/FN1基质活化、MGP应激等描述。主类的P值未转移给内部小群。

## 方法与完整证据

两层使用相同冻结参数独立拟合：全原始文库分母logCP10k、共同实测15,391基因、HVG2000、PCA50、Harmony研究×技术、15近邻，Leiden1主模型/0.5边界敏感性。表达分类采用marker和群体程序；未添加细胞UMI、基因数、线粒体比例、marker数或参考细胞数入口门槛。

推断在同研究、同患者、同部位内比较某类与其他候选，再平均同患者的多个匹配部位。基因使用log2(CP10k+0.1)伪bulk差，程序使用完整来源的表达分数。样本不足以形成生物学重复的行仅留效应，不计算P。全部实际检验的独立患者数不超过15，已用全部符号组合精确计算双侧P；此前Monte Carlo估计留档，未改变类别、患者、统计量或假设家族。基因与程序分别在预定义同类家族内BH校正，完整保留P/FDR。

分类、marker发现及检验来自同一发现数据；名义显著不等于独立验证、稳定新谱系或转移关联已成立。当前新分类可作为实验候选框架，临床A/B1/B2与配对证据仍按真实资料范围解释。

主要文件：
- tables/layer_membership_all11379.tsv.gz：嵌套层成员与背景归属。
- tables/class_annotation_frozen.tsv：查看新P值前冻结的群到主类映射。
- tables/supported_classified_cells.tsv.gz、expanded_classified_cells.tsv.gz：全部逐细胞新分类。
- tables/class_program_statistics.tsv、*_all_gene_statistics.tsv.gz：完整P/FDR。
- tables/independent_patient_comparison_roster.tsv：每项比较的患者与细胞量。
- tables/surface_full_background_donor_screen.tsv：16个表面marker及全部136个单/双组合。
- tables/myCAF_internal_state_review.tsv：两层内部状态复核。
- delivery/FINAL_QA.json：完整性与统计验收。

S525解除内存上限恢复任务仍独立于这两个冻结研究层；本轮未等待或擅自中止该任务。历史分析均保留。
"""
(O/'delivery/双层候选分类与流式线索_20260926.md').write_text(text,encoding='utf-8')
# Keep active project state aligned with the user's two-layer scope.
statefile=P/'PROJECT_STATE.json';backup=P/'PROJECT_STATE.before_phase30_20260926.json';assert not backup.exists();shutil.copy2(statefile,backup)
state=json.loads(statefile.read_text())
state.update(current_stage=O.name,current_task='TWO_LAYER_CAF_CLASSIFICATION_WITH_NOMINAL_EXPLORATORY_SUPPORT',status='PHASE30_TWO_LAYER_CLASSIFICATION_COMPLETE',current_status='PHASE30_TWO_LAYER_CLASSIFICATION_COMPLETE',updated_at=now,dynamic_status_file=str(O/'STATUS.json'),latest_report=str(O/'delivery/双层候选分类与流式线索_20260926.md'),latest_handoff=str(O/'delivery/双层候选分类与流式线索_20260926.md'),next_action='Use supported and expanded layers as joint primary exploratory analyses. Prioritize PDGFRB/MCAM F2 expression and identity validation; nominal P<0.05 may support candidates. Keep114 as reference only.',remaining_acceptance='No unfinished phase30 computation; candidate lineages, metastasis association and flow protein specificity are not independently established.')
state[O.name]={'status':'COMPLETE','path':str(O),'supported':1224,'expanded':8174,'reference_only':114,'classes':['F1','F2','F3','F4','F5'],'nominal_significance_allowed':True,'FDR_retained_not_required_for_candidate_nomination':True,'primary_scope_includes_identity_unresolved_expanded_cells':True}
tmp=statefile.with_name(statefile.name+'.tmp');tmp.write_text(json.dumps(state,ensure_ascii=False,indent=2));tmp.replace(statefile)
(O/'delivery/CODE_SHA256.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (O/'code').glob('*.py')},indent=2))
(O/'DELIVERY_READY').write_text(now+chr(10));print(json.dumps(qa,ensure_ascii=False));print(pivot.to_string())
