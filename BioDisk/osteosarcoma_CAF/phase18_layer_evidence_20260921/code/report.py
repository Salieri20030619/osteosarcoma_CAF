from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
O=Path('C:/codex_test/CAF_evidence_20260921')
def read(n):return pd.read_csv(O/'tables'/n,sep='\t')
L=read('program_evidence_ledger.tsv');C=read('cluster_evidence_cards.tsv');G=read('gene_candidate_evidence.tsv');R=read('cluster_surface_ranking.tsv');CE=read('clinical_program_effects.tsv');PA=read('paired_site_cluster_abundance.tsv');X=read('cross_layer_consistency.tsv');RF=read('separately_refit_program_consistency.tsv')
assert len(L)==260 and len(C)==52 and L[['layer','cluster','program']].duplicated().sum()==0
assert not X.loc[~X.population_changed,'cross_layer_positive_pass'].any()
for r in L[L.evidence_level=='L3'].to_dict('records'):
 fixed=X[(X.layer==r['layer'])&(X.cluster==r['cluster'])&(X.program==r['program'])]
 refit=RF[(RF.layer==r['layer'])&(RF.cluster==r['cluster'])&(RF.program==r['program'])]
 assert fixed.cross_layer_positive_pass.any() or refit.refit_positive_pass.any()
 assert r['evaluable_known_donor_records']>=3 and r['datasets']>=2 and r['positive_fraction']>=.75 and r['median_within_donor_delta']>=.1
unknown=read('unknown_patient_correction.tsv');assert not unknown.counts_as_biological_replicate.any()
effects=read('within_donor_program_effects.tsv');assert not effects.donor_record.str.contains('UNCONFIRMED|UNKNOWN_PATIENT|unknown',case=False).any()
for i,r in C.iterrows():
 genes=G[(G.layer==r.layer)&(G.cluster==r.cluster)&(G.evidence_level=='G2_exploratory_recurrent')].sort_values('median_effect',ascending=False)
 C.loc[i,'recurrent_gene_candidates']=';'.join(genes.symbol.head(8));C.loc[i,'recurrent_gene_count']=len(genes)
 surf=R[(R.layer==r.layer)&(R.cluster==r.cluster)].sort_values('ranking_score',ascending=False)
 if len(surf):
  z=surf.iloc[0];C.loc[i,'top_RNA_surface_combination']=z.markers;C.loc[i,'surface_evidence_level']='S1_RNA_only';C.loc[i,'surface_median_target_coverage']=z.median_coverage;C.loc[i,'surface_median_known_target_fraction']=z.median_known_target_fraction;C.loc[i,'surface_median_CNV_high_fraction']=z.median_CNV_high_fraction
 else:C.loc[i,'surface_evidence_level']='S0_no_evaluable_target'
 C.loc[i,'identity_basis']={'strict':'frozen strict acceptance; not independent fibroblast truth','supported':'at least one low CNV, no high conflict, previous reference checks pass','expanded':'broad molecular candidates; CNV conflicts and unresolved cells retained with labels'}[r.layer]
 C.loc[i,'next_step']='review identity and competing lineage first' if pd.notna(r.identity_risk_tags) and str(r.identity_risk_tags) else ('prioritize source/donor gene review and clinical follow-up' if r.highest_program_level in ['L2','L3'] else 'retain as descriptive candidate; seek stronger replication')
 C.loc[i,'source_tables']='program_evidence_ledger.tsv; gene_candidate_evidence.tsv; clinical_program_effects.tsv; cluster_surface_ranking.tsv'
C.to_csv(O/'tables/cluster_evidence_cards.tsv',sep='\t',index=False)
programs=['myCAF_customer','iCAF_customer','apCAF_customer','Osteogenic_customer','Mineral_customer']
fig,axs=plt.subplots(1,3,figsize=(13,10),gridspec_kw={'width_ratios':[1,1,1]})
cmap=ListedColormap(['#eeeeee','#c6dbef','#6baed6','#08519c']);norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4)
for ax,layer in zip(axs,['strict','supported','expanded']):
 d=L[L.layer==layer].copy();d['level_number']=d.evidence_level.str[1:].astype(int);p=d.pivot(index='cluster',columns='program',values='level_number').reindex(columns=programs).sort_index()
 ax.imshow(p.to_numpy(),cmap=cmap,norm=norm,aspect='auto');ax.set_title(layer+' / program evidence');ax.set_xticks(range(5),['myCAF','iCAF','apCAF','Osteo','Mineral'],rotation=45,ha='right');ax.set_yticks(range(len(p)),['C'+str(z) for z in p.index])
 for y in range(len(p)):
  for x in range(5):v=int(p.iloc[y,x]);ax.text(x,y,'L'+str(v),ha='center',va='center',fontsize=8,color='white' if v==3 else 'black')
fig.suptitle('Internal exploratory evidence only; identity and clinical support assessed separately',fontsize=12)
fig.text(.5,.012,'L0: not evaluable  |  L1: descriptive  |  L2: donor-record recurrence  |  L3: cross-layer sensitivity support',ha='center',fontsize=10)
fig.tight_layout(rect=[0,.04,1,.97]);fig.savefig(O/'figures/program_evidence_levels.png',dpi=180);plt.close(fig)
counts=L.groupby(['layer','evidence_level']).size().unstack(fill_value=0).reindex(['strict','supported','expanded']).reindex(columns=['L0','L1','L2','L3','L4'],fill_value=0).fillna(0).astype(int)
counts.to_csv(O/'tables/evidence_level_counts.tsv',sep='\t')
lines=['# CAF 三层探索与证据支持记录','',
'## 支持层群7和群9形成可继续推进的程序线索，临床与分选证据仍为描述层','',
'已探索严格层10群、支持层18群、广义拓展层24群，覆盖全部52条样本记录的既有联合分析范围。完成260项客户程序证据记录、776条描述性基因候选记录、52群×136种RNA单/双表面组合，以及已有同患者骨—肺和任何转移标签比较。所有群保留，分级用于记录证据，不作为删除探索结果的门槛。','',
'支持层群7的myCAF程序在19条可评估已知供者记录中18条方向为正，中位群内减群外效应+0.548 logCP10k；与广义拓展层对应群的重新拟合对照支持L3。其与严格子集的对照仍只有3条记录、2条同方向，不能将支持层—广义层一致写成严格层完全验证。FN1、COL6A3、TIMP1、CALD1、POSTN等基因也有重复方向支持。',
'支持层群9的矿化和iCAF面板达到L3；myCAF与成骨面板为L2。DCN、SFRP2、LUM、COL3A1、SERPINF1、MGP等表达线索支持继续细化该群。不同程序的证据等级单独记录，不把一个程序的支持转给全部程序。',
'支持层群3的myCAF程序达到L3，但71.7%的细胞带血管壁程序，属于需优先澄清身份的收缩相关群。广义拓展层若干成骨程序也达到L2/L3，同时CNV冲突比例很高；程序重复性不能替代CAF身份确认。','',
'## 分级显示32项程序具备跨供者或跨层探索性支持','',
'|分析层|L0 不可评估|L1 描述|L2 跨供者重复|L3 跨层支持|L4 外部验证|','|---|---:|---:|---:|---:|---:|']
for layer,r in counts.iterrows():lines.append('|'+layer+'|'+'|'.join(str(int(r[k])) for k in counts.columns)+'|')
lines+=['','L1包括单样本/表达描述、效应不一致或尚未达到预设重复规则的记录；负向效应另有明确状态，并不等于“没有表达”。L2要求至少3条已知供者记录、至少2个来源、至少75%群内减群外方向为正，中位差≥0.1 logCP10k。L3在L2基础上，要求真正改变细胞集合的跨层对照，或共享细胞Jaccard≥0.3的另层重新拟合对应群，至少3条共同供者记录中≥75%两侧均正且两侧中位效应≥0.1。完全相同细胞的重复比较不升级。',
'这些是本轮流程明确记录的探索性记账规则，0.1和75%不代表已校准的生物学或临床界值。严格层仍具有原身份筛选证据；本轮严格层程序未达到L2不否定这些细胞已经满足既定身份规则。',
'独立患者跨来源去重未完成，表中使用“已知供者记录”，不能直接作为合并后的独立患者N。程序分级是同来源内部证据，所有L4均为0。','',
'## 临床关联保持在C1描述层，主要方向和未支持结果一并记录','',
'ScPCA两侧各至少10个myCAF的同患者集合：支持层1组、广义拓展层2组、严格层0组。该门槛用于本轮描述的表达稳定性，不是生物学重复数；4组元数据配对全部保留在表中，包括不满足条件者。',
'在支持层唯一满足集合门槛的CCG1106_006配对中，群7两侧各3个myCAF，尚不足以进行该群程序比较；群9两侧均无该群myCAF。不能因整个集合可配对就宣称群7或群9获得了骨—肺重复验证。',
'GSE162454任何转移标签比较中，支持层群7为Yes 2例、No 1例，myCAF程序中位差+0.368；群9为Yes 2例、No 3例，myCAF +0.252、成骨 +0.183、矿化仅+0.004，iCAF −0.176、apCAF −0.255。完整方向保留。这些均为描述性差异，没有把“任何转移”改写为确认的肺转移A/B1分组，也没有计算推断性P值。',
'临床支持代码：C0不可评估；C1可描述的同患者部位或来源内标签比较。本轮未建立独立验证的转移相关性，不将程序L3当作临床L3。','',
'## 表面组合仅达到S1 RNA筛查支持','',
'支持层群7的FAP+THY1组合在16条有≥10个目标细胞的供者/部位记录中，覆盖中位数31.2%，门内已知目标比例3.6%；群9在18条记录中分别为46.5%和11.6%。这些结果提供后续组合改进的起点，当前数据不支持将其作为可直接使用的高纯度分选panel。',
'RNA门控分母为51个有拓展候选样本的全部QC细胞，同一已知供者同一部位的多个标本合并。所有136个组合及完整门内CNV高信号/血管壁比例均保留，未只输出排名靠前者。零目标样本也保留在该分母范围中。未知患者不进入供者汇总。原始目标门控计数与完整来源矩阵逐样本、逐组合一致；无目标支持者记S0，其余仍为S1 RNA-only，未做表面蛋白或分选验证。','',
'## 完整证据总账与复核入口','',
'- `tables/cluster_evidence_cards.tsv`：全部52群的身份基础、风险比例、程序最高等级、临床级别、候选基因和RNA组合。',
'- `tables/program_evidence_ledger.tsv`：260项程序主张，逐项列等级、方向、效应、供者记录数、对应层和相反证据。',
'- `tables/within_donor_program_effects.tsv`、`cross_layer_consistency.tsv`、`separately_refit_program_consistency.tsv`：所有可评估及不可评估的程序效应和分级依据。',
'- `tables/clinical_program_effects.tsv`、`paired_site_cluster_abundance.tsv`：配对与任何转移标签探索，保留零值和不可评估项。',
'- `tables/gene_candidate_evidence.tsv`：描述性基因候选的供者方向复核，G1发现描述/G2重复探索；无独立验证。',
'- `tables/cluster_surface_ranking.tsv`、`cluster_surface_donor_screen.tsv.gz`：每群所有RNA组合、完整分母及混入记录。',
'- `figures/program_evidence_levels.png`：三层程序支持级别图。','',
'## 数据与统计处理记录','',
'本轮修正了上一轮将未确认患者样本共享UNCONFIRMED供者字段的问题：它们保留为逐样本显示记录，但不计入供者重复支持；三层计数和对应样本见unknown_patient_correction.tsv。新的供者证据以本阶段为准，原阶段文件不覆盖。',
'程序分数沿用完整来源基因空间并检查≥80%覆盖；基因候选来自共同实测基因矩阵，排名为描述性群均值差。未增加生物学重复、未做细胞层面的疾病检验、未用聚类清晰度代替临床证据。本轮未做新的聚类参数搜索、未把探索分级作为细胞纳入新阈值。',
'解释边界集中记录：同来源重新拟合及分层并非外部验证；表达程序与聚类基因可有重叠；某些群的血管壁、肿瘤或其他竞争谱系身份仍需复核；配对存在取样时点及治疗差异；无表面蛋白实测。当前可推进的是分层候选线索及其可追溯证据。']
(O/'delivery/三层探索与证据支持记录.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
# Keep a readable complete register rather than selecting only positive programs.
full=['# 全部程序证据登记','', '|层|群|程序|等级|可评估已知供者记录|正向|负向|中位效应|状态|','|---|---|---|---|---:|---:|---:|---:|---|']
for r in L.sort_values(['layer','cluster','program']).to_dict('records'):
 eff='NA' if pd.isna(r['median_within_donor_delta']) else f'{r["median_within_donor_delta"]:.3f}'
 full.append(f'|{r["layer"]}|{r["cluster"]}|{r["program"]}|{r["evidence_level"]}|{r["evaluable_known_donor_records"]}|{r["positive_records"]}|{r["negative_records"]}|{eff}|{r["evidence_status"]}|')
(O/'delivery/全部程序证据登记.md').write_text('\n'.join(full)+'\n',encoding='utf-8')
qa=json.loads((O/'delivery/QA_local.json').read_text());qa.update(remote=json.loads((O/'delivery/QA_remote.json').read_text()),L3_basis_validated=True,unknown_patients_excluded=True,identity_risk_tags_verified=True,all_claims_retained=True,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
assert (O/'markers_surface.exit_code').read_text().strip()=='0'
(O/'delivery/QA.json').write_text(json.dumps(qa,indent=2))
manifest={str(p.relative_to(O)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['tables','delivery','figures','code'] for p in (O/folder).glob('*') if p.is_file() and p.name!='SHA256.json'}
(O/'delivery/SHA256.json').write_text(json.dumps(manifest,indent=2));(O/'ANALYSIS_COMPLETE').write_text('LAYER_EVIDENCE_EXPLORATION_COMPLETE\n')
print(json.dumps(qa),flush=True)
