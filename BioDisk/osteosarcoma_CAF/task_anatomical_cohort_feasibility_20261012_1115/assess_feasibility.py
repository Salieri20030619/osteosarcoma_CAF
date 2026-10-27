import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(name, rows):
    if not rows:
        raise ValueError(f'No records: {name}')
    with (ROOT / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

manifest = read(ROOT / 'clinical_sample_manifest.tsv')
sm = {r['sample']: r for r in manifest}
assert len(sm) == len(manifest)
reg = read(ROOT / 'expression_comparison_registry.tsv')
reg = [r for r in reg if r['sensitivity'] == 'all_cells']
summary = {r['context_id']: r for r in read(ROOT / 'cohort_population_evidence_summary.tsv')}
targets = read(ROOT / 'frozen_target_registry.tsv')
layers = {layer: read(ROOT / f'{layer}_fixed_cells_with_sites.tsv.gz') for layer in ['supported', 'expanded']}
assert len(layers['supported']) == 1224
assert len(layers['expanded']) == 8174

clinical = []
for ds in sorted({r['dataset'] for r in manifest}):
    rows = [r for r in manifest if r['dataset'] == ds]
    for group, pred in [
        ('Bone_all', lambda r: r['strict_site_group'] == 'Bone'),
        ('Lung_all', lambda r: r['strict_site_group'] == 'Lung'),
        ('Bone_documented_lung_metastasis_history', lambda r: r['strict_site_group'] == 'Bone' and r['AB_group'] == 'B1'),
        ('Bone_explicit_A', lambda r: r['strict_site_group'] == 'Bone' and r['AB_group'] == 'A'),
        ('Bone_history_not_established', lambda r: r['strict_site_group'] == 'Bone' and r['AB_group'] not in ['A', 'B1']),
    ]:
        z = [r for r in rows if pred(r)]
        clinical.append(dict(dataset=ds, group=group, samples=len(z), source_scoped_patients=len({r['donor_record'] for r in z}),
                             sample_ids=';'.join(r['sample'] for r in z)))
write('clinical_group_availability.tsv', clinical)

feasibility = []
for r in reg:
    s = summary.get(r['context_id'], {})
    x = dict(r)
    for k in ['bone_cells_by_patient', 'lung_cells_by_patient', 'bone_total_cells', 'lung_total_cells', 'valid_P_features',
              'nominal_up', 'nominal_down', 'FDR_significant', 'nominal_programs']:
        x[k] = s.get(k, '')
    x['interpretation'] = {'testable': 'Exploratory if independent biological N<10; keep effect/P/FDR and cell coverage',
                           'direction_only': 'No valid inferential test with current replication; effect/direction only',
                           'one_or_both_groups_absent': 'Group absent; no expression contrast'}.get(r['status'], r['status'])
    feasibility.append(x)
write('all_class_state_bone_lung_feasibility.tsv', feasibility)
write('testable_class_state_bone_lung.tsv', [r for r in feasibility if r['status'] == 'testable'])

# Audit the proposed third cohort using the existing documented B1 label only.
# Do not infer metastasis-free status from a missing B1 label, or assume timing.
third = []
for t in targets:
    cells = layers[t['layer']]
    key = 'analysis_class' if t['level'] == 'main' else 'analysis_internal'
    z = [c for c in cells if c[key] == t['target']]
    for ds in sorted({r['dataset'] for r in manifest}):
        groups = defaultdict(Counter)
        for c in z:
            m = sm[c['sample']]
            if m['dataset'] != ds:
                continue
            if m['strict_site_group'] == 'Bone':
                groups['Bone_all'][m['donor_record']] += 1
                g = 'Bone_B1' if m['AB_group'] == 'B1' else ('Bone_A' if m['AB_group'] == 'A' else 'Bone_unknown')
                groups[g][m['donor_record']] += 1
            elif m['strict_site_group'] == 'Lung':
                groups['Lung_all'][m['donor_record']] += 1
        b, l = groups['Bone_B1'], groups['Lung_all']
        paired = set(b) & set(l)
        bx, lx = set(b) - paired, set(l) - paired
        if not any(groups.values()):
            continue
        third.append(dict(layer=t['layer'], dataset=ds, target=t['target'], parent=t['parent'], level=t['level'],
                          Bone_all_N=len(groups['Bone_all']), Lung_all_N=len(l), Bone_B1_N=len(b),
                          Bone_explicit_A_N=len(groups['Bone_A']), Bone_unknown_history_N=len(groups['Bone_unknown']),
                          B1_Lung_paired_N=len(paired), B1_Lung_unpaired_B1_N=len(bx), B1_Lung_unpaired_Lung_N=len(lx),
                          B1_cells=sum(b.values()), B1_cells_per_patient=';'.join(f'{d}:{b[d]}' for d in sorted(b)),
                          B1_Lung_paired_testable=len(paired) >= 2,
                          B1_Lung_unpaired_testable=len(bx) >= 2 and len(lx) >= 2,
                          overlapping_Bone_all_vs_B1_invalid=True,
                          timing_note='Documented history; pulmonary metastasis already present at bone sampling not necessarily established'))
write('third_cohort_class_state_availability.tsv', third)

accepted = [r for r in feasibility if r['status'] == 'testable']
assert all(int(r['bone_N']) >= 2 and int(r['lung_N']) >= 2 for r in accepted)
counts = Counter(r['status'] for r in feasibility)
out = dict(checked_at=datetime.now(timezone.utc).isoformat(), status='completed',
           scope='F1/F3/F4 and all existing internal states, original supported/expanded layers',
           selected_primary='Bone_all vs Lung_all, independent of metastatic history',
           third_cohort='Documented B1 bone subgroup; nested in Bone_all, not an independent third arm',
           missing_history_is_not_metastasis_negative=True, patient_units_source_scoped=True,
           existing_P_and_FDR_unchanged=True, new_inferential_tests=0,
           expression_comparisons=len(feasibility), status_counts=dict(counts),
           main_testable_contexts=sum(r['level']=='main' for r in accepted),
           internal_testable_contexts=sum(r['level']=='internal' for r in accepted),
           next_analysis_basis='Reuse completed Phase73 strict_v2 and Phase74 effects, P/FDR and patient/source consistency; no blanket exclusion for FDR failure',
           limitation='Eligibility is per endpoint. Expression feasibility must not be substituted for proportion feasibility.',
           source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'clinical_sample_manifest.tsv', ROOT/'expression_comparison_registry.tsv', ROOT/'frozen_target_registry.tsv', ROOT/'supported_fixed_cells_with_sites.tsv.gz', ROOT/'expanded_fixed_cells_with_sites.tsv.gz', ROOT/'cohort_population_evidence_summary.tsv']})
proportions = [r for r in read(ROOT/'cohort_all_proportion_statistics.tsv') if r['sensitivity'] == 'all_cells']
prop_valid = [r for r in proportions if r['P'] not in ['', 'NA', 'nan']]
out['proportion_records_with_valid_P'] = len(prop_valid)
out['proportion_nominal_records'] = sum(float(r['P']) < 0.05 for r in prop_valid)
assert out['proportion_nominal_records'] == 0
assert len({r['context_id'] for r in reg}) == len(reg)
write('proportion_comparisons_with_valid_P.tsv', prop_valid)
report = [
    'CAF 解剖病灶部位比较可行性与接续结果',
    '核查日期 2026-10-12',
    '',
    '结论：方案1作为主比较。全部明确骨病灶与全部明确肺病灶按解剖取材部位分组，不因转移史缺失或是否转移而排除。F1/F3/F4及原有内部状态中，所有可评估比较均保留。',
    '方案2不是三个独立队列：明确肺转移史的骨病灶是全部骨病灶的子集。保留为补充子集比较，不把同一骨样本重复放入独立三组检验。',
    '若需要互斥三组，应为明确无肺转移的骨病灶、明确有肺转移史的骨病灶、肺病灶；未知病史不能归无转移。现有明确A仅见GSE162454/GSE250015，B1仅见SCPCP000017/SCPCP000023，缺少同来源两类骨病灶的有效对照；不能将来源差异解释为转移史效应。骨取材时与肺转移发生时的先后仍按原始证据保留。',
    '',
    '已完成结果及本次工作',
    '复用已验收Phase73 strict_v2及Phase74，完整保留原始P、FDR、检验家族、方向及细胞覆盖。本次新增逐群可比性与第三队列重叠核查，没有重复运行差异计算或生成新的推断P。',
    'F1/F3/F4主类和固定内部状态共702个来源/层/群体/设计登记组合：10个具备表达检验条件，29个仅作方向描述，663个一侧或两侧群体缺失。702不是独立实验次数或独立患者数。',
    '可检验不代表名义显著。两名配对患者的精确双侧检验最小P为0.5；2对2独立组亦难形成P<0.05，仍保留效应。',
    '',
    '全部10个具备表达检验条件的比较',
    '层\t来源\t类别或状态\t设计\t骨患者N\t肺患者N\t名义上调基因\t名义下调基因',
]
for r in accepted:
    report.append('\t'.join(str(r[k]) for k in ['layer','dataset','target','design','bone_N','lung_N','nominal_up','nominal_down']))
report += [
    '',
    '结果解释',
    '主要来源SCPCP000017：扩展F1为8对3，162条名义基因；原内部C4为8对3，162条名义基因。父子列表共享141条，不能计作独立验证。抗原呈递程序的唯一配对患者方向与独立组结果相反。',
    '支持F3为6对2，758条名义上调基因；功能结果以Reactome线索为主，客户预设面板、Hallmark和Matrisome无名义差异。扩展F3为2对1，缺有效重复；支持层与扩展层不能互相替代。',
    '扩展F4独立组为11对2，794条名义基因。炎症样程序的2名配对患者方向均与独立组相反；其他来源SCPCP000023可比较3对3，未出现名义基因。主要来源肺侧两名患者分别仅1与642个F4细胞。',
    '扩展F4_C9为5对2，137条名义基因并有基质程序线索；肺侧两名患者仅1与3个细胞，应作为稀疏覆盖的探索性结果。F4_C1为2对配对，SCPCP000023的F4_C3为2对2，均保留已完成的阴性结果。',
    '上述骨肺基因及患者级程序均未通过相应FDR阈值。原始P<0.05条目保留为名义显著/探索性结果；缺有效重复不补造P。',
    f'比例端点另行核查：{len(prop_valid)}条原始全细胞口径比例统计具有有效P，P<0.05为0条；比例可比性不能由表达细胞存在与否替代。',
    '',
    '明确肺转移史骨病灶子集',
    '主来源扩展F1只有1名B1患者（3细胞），支持F1为同一患者（2细胞），仅1对方向描述。扩展F3没有B1细胞；支持F3有1名B1患者（1细胞），不具备两组推断条件。扩展F4有2名B1患者（91和11细胞），均有肺取材，可保留2对探索性配对结果。',
    '在当前F1/F3/F4及内部状态里，B1与全部肺病灶的非配对补充比较没有达到每侧有效重复条件；其配对可评估部分已被骨肺配对结果覆盖。',
    '',
    '执行口径',
    '按来源分别分析，以患者为统计单位；同患者两部位作配对并从对应端点非配对分析中移除。保留转移、治疗、细胞覆盖、谱系/CNV和双细胞证据；不新增细胞数或FDR硬门槛。不改变F类及内部标签。',
    '本次只接续已确定的F1/F3/F4队列差异任务，未扩展到F2/F5重算或流式标志筛选。',
    '',
    '文件',
    'all_class_state_bone_lung_feasibility.tsv 完整表达可比性登记及证据；testable_class_state_bone_lung.tsv 全部10项；third_cohort_class_state_availability.tsv 第三队列逐群可比性；clinical_group_availability.tsv 临床台账；proportion_comparisons_with_valid_P.tsv 比例有效检验；FEASIBILITY_RECEIPT.json 本次核查与输入哈希。',
]
assert not any(r['B1_Lung_unpaired_testable'] for r in third)
(ROOT/'COHORT_FEASIBILITY_AND_RESULTS_ZH.txt').write_text('\n'.join(report)+'\n', encoding='utf-8')
(ROOT / 'FEASIBILITY_RECEIPT.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k != 'source_sha256'}, ensure_ascii=False, indent=2))
print('\nTESTABLE EXPRESSION CONTEXTS')
for r in accepted:
    print('\t'.join(str(r[k]) for k in ['layer','dataset','target','design','bone_N','lung_N','nominal_up','nominal_down']))
print('\nDOCUMENTED B1 MAIN-CLASS COUNTS')
for r in third:
    if r['level'] == 'main' and r['Bone_B1_N']:
        print(json.dumps(r, ensure_ascii=False))
print('\nCLINICAL BONE HISTORY')
for r in clinical:
    if r['group'] in ['Bone_documented_lung_metastasis_history','Bone_explicit_A'] and r['samples']:
        print(json.dumps(r, ensure_ascii=False))
