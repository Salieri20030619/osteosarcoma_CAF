from pathlib import Path
import csv, gzip, json, hashlib, collections, datetime
import math

ROOT = Path(__file__).resolve().parent
P = ROOT / 'input_phase89'
OUT = ROOT / 'audit_tables'
OUT.mkdir(exist_ok=True)

def read(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def yes(x): return str(x).lower() == 'true'
def count(rows, key): return dict(collections.Counter(r[key] for r in rows))
def write(name, rows):
    if not rows: return
    with (OUT / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)
def quantile(vals, q):
    vals = sorted(vals); t = (len(vals)-1)*q; i = int(t)
    return vals[i] + (vals[min(i+1,len(vals)-1)]-vals[i])*(t-i)

cell = read(P/'tables/cell_state_mapping.tsv.gz')
groups = read(P/'tables/sample_cluster_state_mapping.tsv')
cv = read(P/'tables/leave_source_out_predictions.tsv')
meta = read(ROOT/'reference_metadata.tsv.gz')
cal = json.loads((P/'resources/CALIBRATION.json').read_text(encoding='utf-8'))
complete = json.loads((P/'MAPPING_COMPLETE.json').read_text(encoding='utf-8'))
checks = {}
checks['unique_query_ids_1268'] = len(cell)==len({r['cell_id'] for r in cell})==1268
checks['completion_counts_match'] = count(cell,'mapping_evidence')==complete['cell_mapping_counts']
checks['unique_groups_51'] = len(groups)==len({(r['sample'],r['cluster']) for r in groups})==51
checks['no_formal_admission'] = all(not yes(r['formal_baseline_admission']) for r in cell)
checks['reference_8174_33states'] = len(meta)==8174 and len({r['subgroup_id'] for r in meta})==33
checks['reference_doublets_26'] = sum(yes(r['doublet_positive_any']) for r in meta)==26
oldquery=read(ROOT.parent/'CAF_expanded_review_20261015/final_cell_review_decisions.tsv.gz')
oldquery={r['cell_id']:r for r in oldquery if yes(r['expanded_review_retained'])}
checks['exact_reviewed_query_id_set']=set(oldquery)=={r['cell_id'] for r in cell}
risk_changes=[]
for r in cell:
    for k,v in oldquery[r['cell_id']].items():
        if k not in r or r[k]!=v:
            # Pandas CSV round-trips can change only the spelling of numeric values.
            try:
                a,b=float(r[k]),float(v)
                if math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12): continue
            except (ValueError,KeyError): pass
            risk_changes.append((r['cell_id'],k,v,r.get(k)))
checks['all_upstream_query_columns_preserved']=not risk_changes
gm = {(r['sample'],r['cluster']):r for r in groups}
bad=[]
for r in cell:
    g=gm[r['sample'],r['cluster']]
    st=g['mapping_evidence']=='state_compatible' and r['template_state']==r['neighbor_state']==g['template_state'] and not yes(r['neighbor_state_tie'])
    cl=g['mapping_evidence'] in ('state_compatible','class_only_compatible') and r['template_class']==r['neighbor_class']==g['template_class'] and not yes(r['neighbor_class_tie'])
    expected='state_compatible' if st else 'class_only_compatible' if cl else 'unresolved'
    if r['mapping_evidence']!=expected: bad.append(r['cell_id'])
checks['cell_decisions_match_recorded_rules']=not bad
calstates=set(cal['calibratable_states'])
badgroups=[]
for g in groups:
    st=g['template_state']==g['neighbor_state'] and not yes(g['neighbor_state_tie']) and float(g['similarity'])>=cal['state_similarity_lower05'] and float(g['margin'])>=cal['state_margin_lower05'] and g['template_state'] in calstates
    cl=g['template_class']==g['neighbor_class'] and not yes(g['neighbor_class_tie']) and float(g['similarity'])>=cal['class_similarity_lower05']
    expected='state_compatible' if st else 'class_only_compatible' if cl else 'unresolved'
    if g['mapping_evidence']!=expected: badgroups.append((g['sample'],g['cluster']))
checks['group_decisions_match_recorded_rules']=not badgroups
quantile_errors={}
for key,values in [('state_similarity_lower05',[float(r['similarity']) for r in cv if yes(r['correct_state'])]),('state_margin_lower05',[float(r['margin']) for r in cv if yes(r['correct_state'])]),('class_similarity_lower05',[float(r['similarity']) for r in cv if yes(r['correct_class'])])]:
    quantile_errors[key]=abs(quantile(values,.05)-cal[key])
    # Source predictions were float32 and the TSV serializes fewer digits.
    checks[key+'_reproduced']=quantile_errors[key]<1e-7
for r in cv:
    st=r['template_state']==r['neighbor_state'] and not yes(r['neighbor_state_tie']) and float(r['similarity'])>=cal['state_similarity_lower05'] and float(r['margin'])>=cal['state_margin_lower05'] and r['template_state'] in calstates
    cl=r['template_class']==r['neighbor_class'] and not yes(r['neighbor_class_tie']) and float(r['similarity'])>=cal['class_similarity_lower05']
    r['audit_rule']='state_compatible' if st else 'class_only_compatible' if cl else 'unresolved'
write('heldout_predictions_with_mapping_rule.tsv',cv)
diagnostics=[]
for state in sorted({r['subgroup_id'] for r in meta}):
    true=[r for r in cv if r['true_state']==state]
    pred=[r for r in cv if r['template_state']==state]
    accepted=[r for r in pred if r['audit_rule']=='state_compatible']
    qc=[r for r in cell if r['mapped_state']==state]
    diagnostics.append(dict(state=state,true_profiles=len(true),true_dual_correct=sum(yes(r['dual_correct_state']) for r in true),predicted_profiles=len(pred),accepted_profiles=len(accepted),accepted_correct=sum(r['true_state']==state for r in accepted),query_state_cells=len(qc),query_reference_fail=sum(not yes(r['reference_QC_pass']) for r in qc),query_single_cell_below_similarity=sum(float(r['similarity'])<cal['state_similarity_lower05'] for r in qc),query_single_cell_below_margin=sum(float(r['margin'])<cal['state_margin_lower05'] for r in qc)))
write('state_diagnostics.tsv',diagnostics)
patient=[]
for ds,pat in sorted({(r['dataset'],r['patient_id']) for r in cell}):
    rr=[r for r in cell if (r['dataset'],r['patient_id'])==(ds,pat)]
    patient.append(dict(dataset=ds,patient_label=pat,libraries='|'.join(sorted({r['sample'] for r in rr})),cells=len(rr),reference_pass=sum(yes(r['reference_QC_pass']) for r in rr),state_compatible=sum(r['mapping_evidence']=='state_compatible' for r in rr),class_only=sum(r['mapping_evidence']=='class_only_compatible' for r in rr),unresolved=sum(r['mapping_evidence']=='unresolved' for r in rr),site='|'.join(sorted({r['site'] for r in rr}))))
write('source_patient_inventory.tsv',patient)
donors=collections.defaultdict(set)
for r in meta: donors[r['donor_record']].add(r['dataset'])
leaks=[]; repeated=[]
for r in cv:
    parts=r['neighbor_donors'].split('|')
    assert len(parts)%2==0
    nn=['|'.join(parts[i:i+2]) for i in range(0,len(parts),2)]
    assert all(d in donors for d in nn)
    if len(nn)!=len(set(nn)): repeated.append(r)
    for d in nn:
        if r['heldout_source'] in donors[d]: leaks.append((r['heldout_source'],d))
checks['no_heldout_source_in_neighbors']=not leaks
checks['unique_neighbor_donor_records']=not repeated
checks['donor_record_single_source']=all(len(v)==1 for v in donors.values())
checks['reference_hash_matches_phase89']=hashlib.sha256((ROOT/'reference_metadata.tsv.gz').read_bytes()).hexdigest()==next(v for k,v in json.loads((P/'resources/INPUT_HASHES.json').read_text(encoding='utf-8'))['before'].items() if 'cells_repaired' in k)
accepted=[r for r in cv if r['audit_rule']=='state_compatible']
classaccepted=[r for r in cv if r['audit_rule']=='class_only_compatible']
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,all_checks_pass=all(checks.values()),quantile_absolute_errors=quantile_errors,quantile_tolerance=1e-7,cv_profiles=len(cv),cv_template_state_correct=sum(yes(r['correct_state']) for r in cv),cv_template_class_correct=sum(yes(r['correct_class']) for r in cv),cv_accepted_state=len(accepted),cv_accepted_state_correct=sum(yes(r['correct_state']) for r in accepted),cv_accepted_class_only=len(classaccepted),cv_accepted_class_only_correct=sum(yes(r['correct_class']) for r in classaccepted),calibration_evaluation_same_profiles=True,query_counts=count(cell,'mapping_evidence'),query_reference_pass=sum(yes(r['reference_QC_pass']) for r in cell),mapped_states=[r for r in diagnostics if r['query_state_cells']],patient_inventory=patient)
result['upstream_column_change_examples']=risk_changes[:5]
(ROOT/'AUDIT_COMPUTATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
