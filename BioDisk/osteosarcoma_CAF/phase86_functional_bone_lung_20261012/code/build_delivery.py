from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase86_functional_bone_lung_20261012'
def read(n):return pd.read_csv(O/'tables'/n,sep='\t',low_memory=False)
def dump(x,p):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
qa=json.loads((O/'QA/INDEPENDENT_QA.json').read_text());assert qa['status']=='PASS'
f=read('all_functional_bone_lung_statistics.tsv');p=read('patient_site_programs.tsv');a=read('all_marker_functional_associations.tsv.gz');counts=read('scope_patient_cell_counts.tsv')
main=f[(f.matrix=='raw')&(f.minimum_cells==1)&(f.design=='independent')&(f.endpoint=='mean_score')].copy()
main.to_csv(O/'tables/PRIMARY_FUNCTIONAL_RESULTS.tsv',sep='\t',index=False)
reference=a[(a.scope=='Fibro_reviewed')&(a.dataset=='SCPCP000017')&a.gene.isin(['GPC1','CD248','THY1','PDGFRA','PDGFRB','FAP','PDPN','CD36','ANTXR1'])].copy()
reference.to_csv(O/'tables/REFERENCE_SURFACE_MARKERS.tsv',sep='\t',index=False)
reg=read('ALL14_retention_registry.tsv')
reg=reg.rename(columns={'surface_annotation':'included_in_existing_surface_annotation'})
reg['annotation_note']='False means absent from this annotation list, not proof of absence from plasma membrane; no antibody suitability inferred'
reg.to_csv(O/'tables/ALL14_retention_registry.tsv',sep='\t',index=False)
refs=[dict(gene='GPC1',source='R&D Systems',product='MAB45192, clone 1019718',url='https://www.rndsystems.com/products/human-glypican-1-antibody-1019718_mab45192',evidence='Official human flow cytometry application; A549 cell-line example; not osteosarcoma CAF validation'),dict(gene='CD248',source='BioLegend',product='Purified anti-human CD248 Endosialin',url='https://www.biolegend.com/en-ie/products/purified-anti-human-cd248-endosialin-antibody-22563',evidence='Official datasheet FC quality tested; not osteosarcoma CAF specificity or panel validation')]
dump(refs,O/'resources/FLOW_REAGENT_REFERENCE.json')
programs=list(json.loads((O/'resources/PROGRAMS.json').read_text()))
contexts=[('Fibro_reviewed','SCPCP000017'),('Fibro_reviewed','SCPCP000023'),('expanded_candidates','SCPCP000017'),('expanded_candidates','SCPCP000023')]
values=np.array([[main[(main.program==pr)&(main.scope==sc)&(main.dataset==ds)].iloc[0].mean_lung_minus_bone for sc,ds in contexts] for pr in programs])
fig,ax=plt.subplots(figsize=(10,8));lim=max(.1,float(np.nanmax(np.abs(values))));im=ax.imshow(values,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
ax.set_yticks(range(len(programs)),[x.replace('_',' ') for x in programs],fontsize=10)
labels=[]
for sc,ds in contexts:
 row=main[(main.scope==sc)&(main.dataset==ds)].iloc[0];labels.append(f'{ds[-3:]} | '+('Fibro' if sc=='Fibro_reviewed' else 'All candidates')+f'\nBone n={int(row.bone_N)}, Lung n={int(row.lung_N)}')
ax.set_xticks(range(4),labels,fontsize=9)
for i,pr in enumerate(programs):
 for j,(sc,ds) in enumerate(contexts):
  r=main[(main.program==pr)&(main.scope==sc)&(main.dataset==ds)].iloc[0]
  ax.text(j,i,f'{r.mean_lung_minus_bone:+.3f}\nP={r.P:.3f}',ha='center',va='center',fontsize=8,color='black')
ax.set_title('Bone vs Lung: patient-level functional scores\nLung minus Bone; no nominal P < 0.05',fontsize=13,pad=16)
fig.colorbar(im,ax=ax,label='Difference in mean log1p(CP10k) score',fraction=.035,pad=.04)
fig.tight_layout();fig.savefig(O/'figures/01_functional_bone_lung.png',dpi=170);fig.savefig(O/'figures/01_functional_bone_lung.pdf');plt.close(fig)
sel=['collagen_assembly','collagen_crosslinking','inflammation','IL6_JAK_STAT3']
fig,axes=plt.subplots(2,2,figsize=(10,7))
for ax,pr in zip(axes.ravel(),sel):
 r=main[(main.scope=='Fibro_reviewed')&(main.dataset=='SCPCP000017')&(main.program==pr)].iloc[0]
 g=p[(p.scope=='Fibro_reviewed')&(p.dataset=='SCPCP000017')&(p.matrix=='raw')&(p.program==pr)]
 for j,(site,ids,color) in enumerate([('Bone',r.bone_patients.split(';'),'#4575b4'),('Lung',r.lung_patients.split(';'),'#d73027')]):
  vals=g[(g.site==site)&g.donor_record.isin(ids)].sort_values('donor_record').mean_score.to_numpy()
  ax.scatter(j+np.linspace(-.09,.09,len(vals)),vals,color=color,s=40,zorder=3);ax.plot([j-.15,j+.15],[vals.mean()]*2,color=color,lw=2)
 ax.set_xticks([0,1],['Bone','Lung']);ax.set_title(pr.replace('_',' ')+f'\nP={r.P:.4f}; FDR={r.FDR:.4f}');ax.set_ylabel('Patient mean score');ax.grid(axis='y',alpha=.2)
fig.suptitle('Reviewed Fibro cells, SCPCP000017: each dot is one patient',fontsize=13);fig.tight_layout();fig.savefig(O/'figures/02_patient_scores.png',dpi=170);fig.savefig(O/'figures/02_patient_scores.pdf');plt.close(fig)
focus=reference[(reference.program=='collagen_assembly')&reference.gene.isin(['GPC1','CD248'])]
report=['骨—肺临床分组下的 CAF 功能状态与表面标志分析（2026-10-12）','',
'结论：本轮已按骨、肺临床分组完成新的患者级计算。预设功能程序及高评分细胞比例尚未显示名义显著的骨—肺差异。发现了可用于研究基质功能状态的表面基因线索，但当前不能将其命名为肺转移相关 CAF 的流式标志。',
'', '一、比较前提与纳入来源',
'“队列/组别”指骨、肺临床分组；SCPCP000017 和 SCPCP000023 是分别检验的两个数据来源。沿用方案一，不按有无肺转移史筛选骨病灶；A/B1/B2 标签保留，未与本轮部位分组混用。SCPCP000018、GSE152048、GSE318840 仅作既往记录，不参与本轮决策。14 个候选全部保留。',
'', '二、骨—肺功能差异尚未建立',
f'本轮固定 14 项程序，包括胶原合成、装配、交联、降解、ECM 组织/降解、收缩、炎症、TGF-beta、IL6-JAK-STAT3 等。主分析为 Fibro 审核细胞，完整 CAF 候选池为身份敏感性分析。合计 {len(f)} 条检验登记（包含重叠的原始/去污染、样本覆盖、配对设计和比例敏感性，并非独立研究），其中 {int(f.P.notna().sum())} 条可计算 P 值，{int(f.P.lt(.05).sum())} 条原始 P<0.05。',
'患者平均程序评分和高评分细胞比例均未出现名义显著条目。这不是仅因 FDR 未通过而否定信号：原始 P 本身未达到筛选标准。未检出不等于两类组织生物学等同。',
'主要 Fibro 原始评分（效应为肺减骨；评分单位为平均 log1p(CP10k)）：',
main[main.scope=='Fibro_reviewed'][['dataset','program','bone_N','lung_N','mean_lung_minus_bone','P','FDR']].to_string(index=False),
'', '三、存在基质状态标志，但与骨—肺关联是不同证据',
'开放检索现有膜蛋白注释及预设候选，并纳入全部14个保留基因。比较同一患者、同一部位内功能评分较高与较低的细胞；每侧至少5个细胞。先在患者内部计算效应，再以独立患者作检验。测试基因如果参与对应程序定义，则不将它的关联作为独立标志证据。',
'S17 的 Fibro 限定分析中，以下两项均由12名患者支持检验，与胶原装配高评分状态呈正向名义关联：']
for r in focus.itertuples():report.append(f'{r.gene}：高评分减低评分的平均表达效应 {r.mean_high_minus_low:+.3f}（log2(CP10k+0.1) 差），P={r.P:.6g}，FDR={r.FDR:.6g}；正向患者 {int(r.positive_patients)}/{int(r.N)}。')
report+=['这说明可以找到描述基质相关状态的标志线索；它们尚不满足本次核心的骨—肺差异关联。程序高低端对比得到的表达差异不等于骨—肺差异，也不是独立外部验证。全部其他名义标志和完整统计均在附件表中，未仅保留这两项。',
'', '四、流式检测处于候选可行性阶段',
'GPC1 和 CD248 均可查到厂家的人源流式应用资料，链接见下。厂家应用不代表在骨肉瘤 CAF 中已经验证，也不证明 CAF 特异性。当前没有经过本项目验证的分选组合；普通流式计数也不能直接测量胶原包绕和组织增厚。',
refs[0]['url'],refs[1]['url'],
'旧膜蛋白注释覆盖不完整。14个候选表中的“未在现有注释中收录”不等于“不是膜蛋白”；例如膜蛋白身份、胞外表位及活细胞抗体适用性须逐项核对。文献 LRRC15 等候选若未收录于该检索表，不能把未被检索解读为 RNA 阴性。',
'', '五、主要限制与解释',
'S23 在当前 Fibro 审核池中只有骨7个、肺19个细胞，主比较为2对2名患者。它不能提供可靠的细分状态标志重复验证；按每侧至少5个细胞作高低评分对照时，也没有足够的患者内部细胞对照。完整候选池结果单列，不冒充纯CAF验证。',
'独立主比较与配对比较按端点可用的患者记录拆分，不将同一患者同时放入骨、肺两个独立组。复核时发现初版采用了更严格的“临床清单双部位患者全部排除”规则，现已与既有Phase81端点可用配对规则统一；初版完整保留为敏感性记录，两种口径均无名义功能差异。',
'功能评分仍可能受到 RNA 检出复杂度、技术和细胞状态组成影响，相关诊断表单独交付；高评分阈值为事先固定的患者等权第75百分位，不表示自然存在的离散亚群。',
'骨—肺部位差异不能代替复发或耐药终点。已保留来源中的病程及治疗字段，未将复发标签混入本轮分组，也未从炎症或 ECM 转录直接推断免疫抑制、耐药或物理包绕。',
'', '六、记录与验收',
'原 F1–F5、原内部簇、细胞身份记录和原始计数均保留。11项继承程序在原始/去污染矩阵中重建误差最大约1e-15；新临床检验及标志检验由独立实现复核，完整FDR家族重现。14个候选仍全部纳入后续研究。',
'本轮完成：分组前提更新、功能状态的骨—肺比较、患者内部表面标志关联搜索、统计复核及结果归档。尚未完成：新增临床样本验证、蛋白流式验证、功能干预实验。']
(O/'delivery/CAF骨肺分组_功能状态与流式候选结果_20261012.txt').write_text('\n'.join(report)+'\n',encoding='utf-8')
summary=dict(phase=O.name,status='analysis_and_QA_complete_visual_review_pending',completed_at=datetime.now(timezone.utc).isoformat(),clinical_groups=['Bone','Lung'],primary_sources=['SCPCP000017','SCPCP000023'],record_only_sources=['SCPCP000018','GSE152048','GSE318840'],programs=14,registered_functional_comparisons=len(f),valid_functional_tests=int(f.P.notna().sum()),nominal_functional_tests=int(f.P.lt(.05).sum()),all14_retained=True,new_computation=True,validated_flow_panel=False,QA=qa,report=str(O/'delivery/CAF骨肺分组_功能状态与流式候选结果_20261012.txt'),output_dir=str(O),primary_pairing='Phase81_endpoint_available',literature_reagents=refs)
dump(summary,O/'delivery/SUMMARY.json')
files=[q for folder in ['tables','resources','QA','figures','delivery','code'] for q in (O/folder).rglob('*') if q.is_file() and q.suffix not in ['.zip','.pyc']]
mf={str(q.relative_to(O)):dict(bytes=q.stat().st_size,sha256=sha(q)) for q in files};dump(mf,O/'delivery/FILE_MANIFEST.json')
archive=O/'delivery/CAF_functional_bone_lung_20261012.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for name in mf:z.write(O/name,name)
 z.write(O/'delivery/FILE_MANIFEST.json','delivery/FILE_MANIFEST.json')
summary.update(archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size)
dump(summary,O/'DELIVERY_READY.json');print(json.dumps(summary,ensure_ascii=False))
