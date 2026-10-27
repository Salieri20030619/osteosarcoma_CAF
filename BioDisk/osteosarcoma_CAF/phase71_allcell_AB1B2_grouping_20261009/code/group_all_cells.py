from pathlib import Path
import csv, gzip, json, hashlib, collections, datetime, contextlib, shutil, traceback, zipfile
import h5py

P = Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O = P / 'phase71_allcell_AB1B2_grouping_20261009'
GROUPS = ['A', 'B1', 'B2', 'Unassigned']

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(f):
    h = hashlib.sha256()
    with open(f, 'rb') as s:
        for block in iter(lambda: s.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def read(f):
    op = gzip.open if str(f).endswith('.gz') else open
    with op(f, 'rt', encoding='utf-8-sig') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def dump(f, obj):
    f.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def table(f, rows, fields=None):
    rows = list(rows)
    fields = fields or list(rows[0])
    with f.open('w', encoding='utf-8', newline='') as h:
        w = csv.DictWriter(h, fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def group(label):
    if label == 'A_reported_nonmetastasis':
        return 'A', 'reported_nonmetastatic_followup_not_established'
    if label == 'A_candidate_reported_no_metastasis_timing_unknown':
        return 'A', 'reported_nonmetastatic_candidate_timing_unknown'
    if label == 'B1_history_supported_timing_unresolved':
        return 'B1', 'same_patient_lung_metastasis_history_sampling_order_unresolved'
    if label == 'B2':
        return 'B2', 'documented_lung_metastasis_specimen'
    return 'Unassigned', 'insufficient_or_outside_AB1B2_evidence:' + label

def run():
    assert not (O / 'COMPLETE.json').exists(), 'Completed output exists; do not overwrite'
    for d in ['tables', 'groups', 'resources', 'QA', 'delivery', 'logs', 'code']:
        (O / d).mkdir(parents=True, exist_ok=True)
    sources = {
        'clinical': P / 'phase33_client_report_20260926/tables/final_clinical_sample_manifest.tsv',
        'strict_site': P / 'phase44_patient_mechanism_20261001/clinical_reporting_extension/tables/clinical_manifest_with_strict_site_group.tsv',
        'annotation': P / 'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz',
        'source_manifest': P / 'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',
        'pairs': P / 'phase31_customer_completion_20260926/resources/bone_lung_pairs_frozen.tsv',
        'phase66': P / 'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',
    }
    hashes = {k: {'path': str(f), 'bytes': f.stat().st_size, 'sha256': sha(f)} for k, f in sources.items()}
    dump(O / 'resources/INPUT_MANIFEST.json', hashes)
    for k in ['clinical', 'strict_site', 'source_manifest', 'pairs']:
        shutil.copy2(sources[k], O / 'resources' / sources[k].name)
    protocol = {
        'frozen_at': now(), 'scope': 'All 876648 source-called cells in all 107 discovery samples',
        'operation': 'Sample-level clinical mapping inherited by every cell; no expression-based clinical grouping',
        'A': ['A_reported_nonmetastasis', 'A_candidate_reported_no_metastasis_timing_unknown'],
        'B1': ['B1_history_supported_timing_unresolved'], 'B2': ['B2'],
        'unknown_policy': 'All other original labels remain Unassigned with original evidence',
        'timing_policy': 'A is reported status, not proven lifetime absence; B1 is history-supported, not confirmed synchronous metastasis',
        'no_cell_filter': ['lineage', 'CNV', 'doublet', 'UMI', 'genes', 'mitochondrial fraction', 'old8174 membership'],
        'source_annotations': 'Phase55 annotations retained as reference; Phase66 membership and its repaired doublet evidence added separately',
        'patient_key': 'dataset plus verified original donor_record; blinded patient_id is not the clinical identifier',
        'pairing': 'Same-source donor with B1 and B2; multiple same-side specimens form one patient-side; contemporaneity not established',
        'cross_study_patient_independence': 'Not established; report source-scoped donor records',
        'matrix_policy': 'Group-specific sample manifests point to existing full_called_counts.h5ad; exact cell IDs verified; counts not copied or changed',
        'no_new_clustering_or_inferential_statistics': True,
    }
    dump(O / 'resources/GROUPING_PROTOCOL_FROZEN.json', protocol)
    clinical = read(sources['clinical'])
    cm = {r['sample']: r for r in clinical}
    sm = {r['sample']: r for r in read(sources['source_manifest'])}
    strict = {r['sample']: r for r in read(sources['strict_site'])}
    assert len(cm) == len(clinical) == len(sm) == 107 and set(cm) == set(sm) == set(strict)
    for s in cm:
        assert cm[s]['clinical_group'] == strict[s]['clinical_group']
        if cm[s]['clinical_group'] == 'B2':
            assert strict[s]['strict_site_group'] == 'Lung'
    overlay_rows = read(sources['phase66'])
    overlay = {r['cell_id']: r for r in overlay_rows}
    assert len(overlay) == len(overlay_rows) == 8174
    counts = collections.Counter()
    lineage_counts = collections.Counter()
    sample_ids = collections.defaultdict(set)
    seen_original, seen_blind = set(), set()
    overlay_found = 0
    extra = ['AB_group', 'AB_evidence_status', 'clinical_group_original', 'clinical_patient_id', 'clinical_donor_record', 'clinical_site', 'clinical_strict_site', 'clinical_treatment', 'clinical_disease_timing', 'clinical_primary_or_metastasis', 'clinical_group_evidence', 'phase66_candidate_member', 'phase66_doublet_assessed', 'phase66_doublet_positive']
    index_fields = ['original_cell_id', 'cell_id', 'sample', 'dataset', 'clinical_patient_id', 'clinical_donor_record', 'AB_group', 'AB_evidence_status', 'clinical_group_original']
    with contextlib.ExitStack() as stack:
        src = stack.enter_context(gzip.open(sources['annotation'], 'rt', encoding='utf-8'))
        rd = csv.DictReader(src, delimiter='\t')
        assert not set(extra) & set(rd.fieldnames)
        out = stack.enter_context(gzip.open(O / 'tables/fullcell_group_annotations.tsv.gz', 'wt', encoding='utf-8', newline='', compresslevel=3))
        master = csv.DictWriter(out, rd.fieldnames + extra, delimiter='\t')
        master.writeheader()
        writers = {}
        for g in GROUPS:
            h = stack.enter_context(gzip.open(O / 'groups' / (g + '_cells.tsv.gz'), 'wt', encoding='utf-8', newline='', compresslevel=3))
            writers[g] = csv.DictWriter(h, index_fields, delimiter='\t', extrasaction='ignore')
            writers[g].writeheader()
        for i, r in enumerate(rd, 1):
            s = r['sample']; c = cm[s]; g, ev = group(c['clinical_group'])
            cid = r['original_cell_id']; bid = r['cell_id']
            assert cid not in seen_original and bid not in seen_blind
            assert r['dataset'] == c['dataset'] and r['original_donor_record'] == c['donor_record']
            seen_original.add(cid); seen_blind.add(bid); sample_ids[s].add(cid)
            o = overlay.get(cid)
            if o:
                assert o['sample'] == s
                overlay_found += 1
            r.update(dict(zip(extra, [g, ev, c['clinical_group'], c['patient_id'], c['donor_record'], c['site'], strict[s]['strict_site_group'], c['treatment'], c['provider_disease_timing'], c['primary_or_metastasis'], c['group_evidence'], str(o is not None), o['doublet_assessed_any'] if o else '', o['doublet_positive_any'] if o else ''])))
            master.writerow(r); writers[g].writerow(r)
            counts[s] += 1
            lineage_counts[(g, c['dataset'], r['review_lineage_code'])] += 1
            if i % 200000 == 0:
                print('Mapped', i, 'cells', flush=True)
    assert sum(counts.values()) == len(seen_original) == len(seen_blind) == 876648
    assert set(counts) == set(cm) and overlay_found == 8174
    print('Verifying exact cell IDs against 107 original sample matrices', flush=True)
    sample_rows, matrix_qa = [], []
    for s, c in cm.items():
        assert counts[s] == int(c['full_called_cells'])
        path = P / 'phase26_client_reanalysis_20260922/inputs' / s / 'full_called_counts.h5ad'
        with h5py.File(path, 'r') as h:
            key = h['obs'].attrs['_index']
            if isinstance(key, bytes): key = key.decode()
            ids = h['obs'][key].asstr()[:]
            assert len(ids) == len(set(ids)) == counts[s] and set(ids) == sample_ids[s]
            assert 'X' in h
        g, ev = group(c['clinical_group'])
        row = dict(c)
        row.update(AB_group=g, AB_evidence_status=ev, strict_site_group=strict[s]['strict_site_group'], grouped_cells=counts[s], full_called_counts_h5ad=str(path), original_source_counts_path=sm[s]['counts_path'])
        sample_rows.append(row)
        matrix_qa.append({'sample': s, 'n_cells': counts[s], 'exact_ID_match': True, 'path': str(path), 'bytes': path.stat().st_size})
    table(O / 'tables/sample_group_manifest.tsv', sample_rows)
    table(O / 'QA/sample_matrix_ID_verification.tsv', matrix_qa)
    table(O / 'tables/unassigned_sample_manifest.tsv', [r for r in sample_rows if r['AB_group'] == 'Unassigned'])
    patient = collections.defaultdict(list)
    for r in sample_rows: patient[(r['dataset'], r['donor_record'], r['AB_group'])].append(r)
    patient_rows = []
    for (ds, donor, g), rr in sorted(patient.items()):
        patient_rows.append({'dataset': ds, 'donor_record': donor, 'AB_group': g, 'samples': ';'.join(r['sample'] for r in rr), 'n_samples': len(rr), 'cells': sum(r['grouped_cells'] for r in rr), 'sites': ';'.join(sorted(set(r['site'] for r in rr)))})
    table(O / 'tables/patient_group_manifest.tsv', patient_rows)
    summary = []
    study_rows = []
    for g in GROUPS:
        z = [r for r in sample_rows if r['AB_group'] == g]
        summary.append({'AB_group': g, 'samples': len(z), 'source_scoped_donor_records': len(set(r['donor_record'] for r in z)), 'cells': sum(r['grouped_cells'] for r in z), 'studies': ';'.join(sorted(set(r['dataset'] for r in z)))})
        table(O / 'groups' / (g + '_sample_matrix_manifest.tsv'), z)
        for ds in sorted(set(r['dataset'] for r in z)):
            q = [r for r in z if r['dataset'] == ds]
            study_rows.append({'AB_group': g, 'dataset': ds, 'samples': len(q), 'source_scoped_donors': len(set(r['donor_record'] for r in q)), 'cells': sum(r['grouped_cells'] for r in q)})
    table(O / 'tables/group_summary.tsv', summary)
    table(O / 'tables/study_group_summary.tsv', study_rows)
    table(O / 'tables/group_lineage_reference_counts.tsv', [dict(AB_group=g, dataset=ds, reference_lineage=lin, cells=n) for (g, ds, lin), n in sorted(lineage_counts.items())])
    pairs = []
    for r in read(sources['pairs']):
        a = r['bone_samples'].split(';'); b = r['lung_samples'].split(';')
        assert all(cm[x]['donor_record'] == r['donor_record'] for x in a + b)
        assert all(group(cm[x]['clinical_group'])[0] == 'B1' for x in a)
        assert all(group(cm[x]['clinical_group'])[0] == 'B2' for x in b)
        pairs.append(dict(r, B1_full_cells=sum(counts[x] for x in a), B2_full_cells=sum(counts[x] for x in b), statistical_unit='one_source_scoped_patient_not_one_sample_or_cell'))
    derived_pairs = {d for _, d, g in patient if g == 'B1'} & {d for _, d, g in patient if g == 'B2'}
    assert derived_pairs == {r['donor_record'] for r in pairs} and len(pairs) == 4
    table(O / 'tables/paired_B1_B2_patient_manifest.tsv', pairs)
    assert {r['AB_group']: r['cells'] for r in summary} == {'A':45706, 'B1':37610, 'B2':303543, 'Unassigned':489789}
    qa = {'checked_at': now(), 'status': 'PASS', 'cells': 876648, 'samples': 107, 'all_input_cells_retained': True, 'original_IDs_unique': True, 'blinded_IDs_unique': True, 'all_107_sample_matrix_ID_sets_match': True, 'group_assignment_independent_of_cell_identity': True, 'phase66_metadata_matched_cells': overlay_found, 'no_counts_modified': True, 'no_P_values_generated': True, 'source_scoped_pairs': 4, 'group_summary': summary}
    dump(O / 'QA/ACCEPTANCE.json', qa)
    (O / 'delivery/README.txt').write_text('全细胞 A/B1/B2 初步分组\n\n全部107份样本、876648个细胞均保留。临床分组由样本元数据映射到细胞，不使用CAF、CNV、双细胞或表达阈值筛选。\nA包括报告无转移的1例及报告无转移但时点不明的4例，保留原始标签。B1为同患者有肺转移样本支持的骨样本，未确认同步转移。B2为明确肺转移灶。其余样本保留为Unassigned，不强行填补。\n\n表格入口：tables/sample_group_manifest.tsv、group_summary.tsv、patient_group_manifest.tsv、paired_B1_B2_patient_manifest.tsv。全细胞注释：tables/fullcell_group_annotations.tsv.gz。各组精确细胞索引与完整表达矩阵路径在groups/。分组不需要重复复制原始表达矩阵：各样本整体属于一个组，107份原始矩阵的细胞ID已逐一核对。原始计数未修改。\n\nPhase55旧表达/CNV注释仅作为参考并保留；额外标记Phase66候选集合及其修复后的双细胞证据，不以这些字段排除细胞。临床患者ID使用clinical_patient_id/clinical_donor_record；原patient_id为盲审编号。\n\n患者数按来源内记录报告，跨来源独立性尚未完全核实。4条B1-B2配对记录分属两个研究各2条；多份肺样本不增加患者数。\n本轮只分组，尚未重新识别CAF、聚类、差异表达、富集或流式marker筛选。\n', encoding='utf-8')
    files = [f for sub in ['tables', 'groups', 'resources', 'QA', 'code', 'delivery'] for f in (O/sub).rglob('*') if f.is_file() and f.suffix != '.zip']
    manifest = [{'path':str(f.relative_to(O)), 'bytes':f.stat().st_size, 'sha256':sha(f)} for f in files]
    dump(O / 'delivery/FILE_MANIFEST.json', manifest)
    archive = O / 'delivery/CAF_allcell_AB1B2_grouping_20261009.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=1) as z:
        for f in files + [O / 'delivery/FILE_MANIFEST.json']:
            z.write(f, str(f.relative_to(O)))
    with zipfile.ZipFile(archive) as z: assert z.testzip() is None
    completed = {'completed_at': now(), 'status': 'completed_verified', 'scope': 'All-cell clinical grouping only', 'acceptance': str(O / 'QA/ACCEPTANCE.json'), 'archive': str(archive), 'archive_sha256':sha(archive), 'group_summary': summary}
    dump(O / 'COMPLETE.json', completed)
    print(json.dumps(completed, ensure_ascii=False), flush=True)

if __name__ == '__main__':
    try:
        run()
    except BaseException:
        if O.exists(): (O / 'exit_code').write_text('1\n')
        traceback.print_exc()
        raise
    else:
        (O / 'exit_code').write_text('0\n')
