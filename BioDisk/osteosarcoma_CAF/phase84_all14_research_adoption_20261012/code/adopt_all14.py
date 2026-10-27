from pathlib import Path
from datetime import datetime, timezone
from zipfile import ZipFile, ZIP_DEFLATED
import csv, json, hashlib, shutil

ROOT=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
PHASE='phase84_all14_research_adoption_20261012'
OUT=ROOT/PHASE
for sub in ['tables','resources','delivery','state_backups']:
    (OUT/sub).mkdir(parents=True,exist_ok=True)
now=datetime.now(timezone.utc).isoformat()
stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
genes=['NEBL','RTKN','NDUFS7','PNKD','TOR1AIP1','GRID2','TAF10','ZNF526','BEST1','CALHM5','PEF1-AS1','RGS2-AS1','AL096854.1','EBF1']
main_sources=['SCPCP000017','SCPCP000023']
record_sources=['SCPCP000018','GSE152048','GSE318840']
def read_tsv(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))
def write_json(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def write_tsv(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
baseline=ROOT/'phase82_allCAF_crossdataset_markers_20261012/tables/MARKER_EVIDENCE_REGISTRY.tsv'
literature=ROOT/'phase83_upstream_and_direction_20261012/tables/ALL14_LITERATURE_DIRECTION_INTEGRATION.tsv'
stats={r['symbol']:r for r in read_tsv(baseline)}
lit={r['gene']:r for r in read_tsv(literature)}
assert set(genes)==set(stats)==set(lit), (list(stats),list(lit))
questions={
'NEBL':('骨架与黏附','核对 NEBL/LASP2 亚型，研究其表达与 CAF 黏附、迁移表型的关系。'),
'RTKN':('骨架与黏附','围绕 Rho 相关骨架与迁移功能设计研究；基因身份限定为 RTKN，不借用 RTKN2 证据。'),
'NDUFS7':('线粒体与代谢','研究呼吸链复合体 I 相关代谢状态；区分线粒体通路背景证据与 NDUFS7 单基因功能证据。'),
'PNKD':('线粒体与代谢','以既往其他癌种的代谢研究提出乳酸及相关信号假设，检验其是否适用于 CAF。'),
'TOR1AIP1':('核结构与转录调控','区分 LAP1 亚型，研究核形态、迁移与组织定位；记录不同肿瘤模型的作用差异。'),
'GRID2':('膜蛋白与信号探索','先确认 CAF 中的表达和细胞定位，再探索受体相关功能；不把预后关联当作自噬机制证据。'),
'TAF10':('核结构与转录调控','研究其与转录状态及增殖表型的关系，评估既往肿瘤模型机制在 CAF 中的适用性。'),
'ZNF526':('核结构与转录调控','研究表达与 CAF 状态及存活表型的关系，区分遗传变异关联和表达效应。'),
'BEST1':('膜蛋白与信号探索','研究细胞定位及离子、体积调节与增殖相关表型；不预设骨基质功能。'),
'CALHM5':('膜蛋白与信号探索','先确认表达、定位及检测可行性，再探索功能；不借用其他 CALHM 家族成员的机制。'),
'PEF1-AS1':('lncRNA 探索','按稳定基因标识确认转录本、表达和定位，随后研究 CAF 状态关联；不借用 PEF1 功能。'),
'RGS2-AS1':('lncRNA 探索','统一 AL390957.1/LINC01032/RGS2-AS1 标识，确认转录本及定位，再开展功能研究。'),
'AL096854.1':('lncRNA 探索','结合既往骨肉瘤关联记录，研究其细胞来源、表达及 CAF 状态关系；不预设预后或转移功能。'),
'EBF1':('核结构与转录调控','结合骨肉瘤细胞及骨间质研究，研究其与骨病灶 CAF 状态的关系，区分不同细胞背景。')}
rows=[]
for gene in genes:
    s,l=stats[gene],lit[gene]
    row=dict(gene=gene,feature_ids=s['feature_ids'],aliases=l['aliases'],retain=True,next_stage_status='纳入下一步研究',decision_sources=';'.join(main_sources),direction=s['direction'],evidence_label='两个数据集同向名义显著；探索性候选',research_axis=questions[gene][0],next_research_question=questions[gene][1])
    for source in main_sources:
        for col in ['bone_N','lung_N','mean_lung_minus_bone','P','FDR','family']:
            key=source+'__'+col;row[key]=s[key]
    for key in ['evidence_type','literature_context','literature_finding','transfer_boundary','reference_titles','reference_identifiers','reference_urls','searched_at']:
        row[key]=l[key]
    row['incomplete_sources_used_for_retention_or_ranking']=False
    row['literature_role']='确定研究问题和实验依据；文献不足不剔除候选'
    rows.append(row)
write_tsv(OUT/'tables/ALL14_NEXT_STAGE_RESEARCH_REGISTRY.tsv',rows)
write_json(OUT/'tables/ALL14_NEXT_STAGE_RESEARCH_REGISTRY.json',rows)
roles=[dict(source=x,role='主要统计依据' if x in main_sources else '仅记录',use_for_retention=x in main_sources,use_for_ranking=x in main_sources,note='保留原始 P、FDR、效应及样本量' if x in main_sources else '不以其同向、反向、零差异或缺测升降级候选') for x in main_sources+record_sources]
write_tsv(OUT/'tables/DATA_SOURCE_ROLES.tsv',roles)
policy=dict(version='2026-10-12_all14_user_adoption',requested_action='全部 14 个候选保留并纳入下一步研究',adopted_at=now,comparison='所有 CAF 候选细胞合并；按骨骼病灶与肺病灶分组，不按患者是否发生肺转移再筛选',included_genes=genes,excluded_genes=[],primary_statistical_sources=main_sources,record_only_sources=record_sources,record_only_sources_affect_retention=False,record_only_sources_affect_ranking=False,literature_role='优先用于研究方向及机制假设；按骨肉瘤/CAF/其他癌种、功能实验/关联/注释区分证据',literature_absence_is_exclusion=False,superseded_recommendation='phase83 中依据三个不完整来源方向一致或相反所作的候选升降级建议不再生效；历史结果保留',statistical_rule='两个主要来源各自 P < 0.05 为名义显著；保留 FDR，不将名义显著写成校正后显著',status='研究名单及方案已登记；尚未执行新增功能实验',baseline_cell_pool_changed=False)
write_json(OUT/'resources/RESEARCH_SELECTION_POLICY.json',policy)
for p in [baseline,literature]:shutil.copy2(p,OUT/'resources'/('HISTORICAL_'+p.name))
phase83=ROOT/'phase83_upstream_and_direction_20261012'
for rel in ['tables/MARKER_DIRECTION_SUMMARY.tsv','tables/SCPCP000018_late_recovery_direction_sensitivity.tsv','resources/LITERATURE_CURATED_RECORDS.json','resources/LITERATURE_REFERENCES.tsv','resources/LITERATURE_SEARCH_TERMS.tsv']:
    p=phase83/rel
    if p.exists():shutil.copy2(p,OUT/'resources'/('ARCHIVAL_'+p.name))
(OUT/'resources/ARCHIVE_README.txt').write_text('历史方向结果、上游核查及原始联合解释仅用于追溯。三个不完整来源不参与当前候选保留或优先级判断。当前决策以 RESEARCH_SELECTION_POLICY.json 和 ALL14_NEXT_STAGE_RESEARCH_REGISTRY 为准。主要来源的配对及细胞身份敏感性结果仍是科学记录，不据此取消本次全部纳入决定。\n',encoding='utf-8')
parts=['14 个 CAF 候选全部进入下一阶段研究（2026-10-12）','',
'本次决定：全部 14 个候选保留。三个人群信息或重复不足的数据集仅作记录，不参与候选取舍和优先级判断。文献用于确定研究方向，不以缺少既往文献剔除候选。',
'纳入名单：'+ '、'.join(genes)+'。',
'主要统计依据：SCPCP000017 与 SCPCP000023；两者均按解剖病灶部位比较骨骼与肺病灶。原始 P、FDR、效应量及独立样本数完整保留。14 个候选均为两来源同向名义显著，当前并非 FDR 显著。',
'仅作记录：SCPCP000018、GSE152048、GSE318840。既往根据这三个来源方向结果作出的升降级建议撤回，包括对 RTKN、TAF10、TOR1AIP1 的方向性降级。历史数值和上游核查不删除。',
'文献使用：优先参考与研究问题匹配的原始研究；骨肉瘤或骨间质研究、CAF 直接研究、其他癌种功能研究、表达/预后关联和基因注释分别注明。文献为研究提供依据，不替代本项目细胞中的功能验证。',
'', '下一阶段研究安排',
'1. 全部 14 个基因进入共同的表达及细胞来源核实：保留患者层面的两来源效应，确认 CAF 身份、组织定位和转录本/亚型。现有细胞身份及配对敏感性记录继续保留，不作为本次剔除条件。',
'2. 按下列研究问题组织功能研究准备。分组表示研究方向，不表示候选等级；新增实验尚未执行。']
for r in rows:parts.append(f"   {r['gene']}（{r['research_axis']}）：{r['next_research_question']}")
parts+=['3. 后续结果对全部 14 个候选分别记录，区分表达支持、定位支持和功能支持，不以三个不完整来源的方向结果重新设置门槛。','', '逐基因文献证据与统计信息']
for r in rows:
    parts += ['',r['gene'],r['literature_finding'],'证据适用范围：'+r['transfer_boundary'],'参考：'+r['reference_urls']]
parts+=['','本次交付范围：更新候选研究名单、来源角色、逐基因研究问题和项目恢复记录；未重算差异表达，未修改原始数据或冻结的 8174 细胞集合，未声称功能验证完成。']
report=OUT/'delivery/14个CAF候选_全部纳入下一阶段研究_20261012.txt'
report.write_text('\n'.join(parts)+'\n',encoding='utf-8')
# Verify adoption and preservation of existing statistical evidence.
assert len(rows)==14 and all(r['retain'] for r in rows)
assert all(not r['incomplete_sources_used_for_retention_or_ranking'] for r in rows)
assert all(not r['use_for_retention'] and not r['use_for_ranking'] for r in roles if r['source'] in record_sources)
for r in rows:
    for source in main_sources:
        for col in ['bone_N','lung_N','mean_lung_minus_bone','P','FDR','family']:
            key=source+'__'+col;assert r[key]==stats[r['gene']][key]
        assert float(r[source+'__P'])<.05
assert sha(baseline)==sha(OUT/'resources'/('HISTORICAL_'+baseline.name))
qa=dict(status='passed',all14_retained=True,record_only_sources_excluded_from_decision=True,primary_statistics_unchanged=True,original_files_unchanged=True,checked_at=now)
write_json(OUT/'delivery/QA.json',qa)
record=dict(phase=PHASE,status='completed_verified_research_roster_adopted',updated_at=now,output_dir=str(OUT),retained_genes=genes,retained_count=14,excluded_count=0,selection_policy=policy,report=str(report),research_registry=str(OUT/'tables/ALL14_NEXT_STAGE_RESEARCH_REGISTRY.tsv'),qa=qa,new_functional_experiments_executed=False,local_delivery_verified=False)
states=['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_CLASSIFICATION_DIFFERENCE_RESUME_CURRENT.json','CAF_ALLCAF_CROSSDATASET_MARKERS_CURRENT.json','CAF_UPSTREAM_DIRECTION_LITERATURE_CURRENT.json']
write_json(OUT/'delivery/COMPLETED.json',record)
mf={str(p.relative_to(OUT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for folder in ['tables','resources','delivery'] for p in (OUT/folder).rglob('*') if p.is_file()}
write_json(OUT/'delivery/FILE_MANIFEST.json',mf)
archive=OUT/'delivery/CAF_all14_research_adoption_20261012.zip'
with ZipFile(archive,'w',ZIP_DEFLATED) as z:
    for name in list(mf)+['delivery/FILE_MANIFEST.json']:z.write(OUT/name,name)
record.update(archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size)
for name in states:
    p=ROOT/name
    if not p.exists():continue
    shutil.copy2(p,OUT/'state_backups'/f'{stamp}_{name}')
    obj=json.loads(p.read_text(encoding='utf-8-sig'))
    obj['all14_research_adoption']=record
    obj['active_research_selection_policy']=policy
    obj['previous_incomplete_source_direction_rankings_status']='superseded_record_only'
    if name in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json']:
        for key in ['current_phase','current_stage','phase']:obj[key]=PHASE
        for key in ['latest_task','last_completed_task','last_completed_analysis','latest_research_results']:obj[key]=record
        obj['active_task']=dict(status='completed',phase=PHASE,description='14 个候选全部纳入；研究名单及逐基因研究问题已登记')
    obj['updated_at']=now
    write_json(p,obj)
write_json(ROOT/'CAF_ALL14_RESEARCH_CURRENT.json',record)
print(json.dumps(record,ensure_ascii=False))
