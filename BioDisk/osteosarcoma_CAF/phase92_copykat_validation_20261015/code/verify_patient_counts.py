from pathlib import Path
import datetime, hashlib, json
import pandas as pd

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase92_copykat_validation_20261015'
out=O/'patient_count_verification';out.mkdir(exist_ok=True)
run=json.loads((O/'COMPUTATION_RECONCILED.json').read_text(encoding='utf-8'))
assert run['computation_exit_code']==0 and not run['missing_samples']
assert len(run['sample_accounting'])==10
old=pd.read_csv(P/'audit_expanded_admission_20261015/tables/final_cell_review_decisions.tsv.gz',sep='\t',encoding='utf-8').set_index('cell_id')
lung=pd.read_csv(P/'phase91_lung_patient_validation_20261015/tables/all_cell_validation.tsv.gz',sep='\t',encoding='utf-8').set_index('cell_id')
assert old.index.is_unique and lung.index.is_unique
rows=[];cell_rows=[];hashes={}
for sm,receipt in sorted(run['sample_accounting'].items()):
    assert receipt['exit_code']==0 and receipt['completion_marker']
    f=O/'tables'/(sm+'_copykat_full_evidence.tsv.gz')
    c=pd.read_csv(f,sep='\t',encoding='utf-8').set_index('cell_id')
    assert c.index.is_unique and c.patient_id.nunique()==1 and set(c['sample'])=={sm}
    assert len(c)==receipt['full_cells']
    v=json.loads((O/'resources'/(sm+'_VALIDITY_RECHECK.json')).read_text(encoding='utf-8'))
    assert v['sample']==sm and v['full_cells']==len(c) and v['patient_id']==str(c.patient_id.iloc[0])
    valid=bool(v['CopyKAT_valid_for_binary_identity'])
    is_lung=sm.startswith('GSE')
    support=c.full_scDblFinder_class.eq('singlet') & c.native_copykat_call.eq('diploid') & c.provisional_lineage.eq('Fibroblast')
    if is_lung:
        a=lung.reindex(c.index);assert a['sample'].notna().all()
        reviewed=a.cluster.astype(str).isin(['5','7','11'])
        ambient=a.decontX_lineage_winner.eq('Fibroblast') & a.decontX_detected_Fib_genes_ge1.ge(2)
        strict=pd.Series(False,index=c.index)
        source='GSE234187';site='Lung'
    else:
        a=old.reindex(c.index)
        reviewed=a.expanded_review_retained.eq(True)
        ambient=a.corrected_lineage_supported.eq(True) & a.corrected_det_Fib.ge(2)
        strict=a.supported_CAF_candidate.eq(True)
        support &= ~a.audit_doublet_positive.eq(True)
        source='NCC_Zheng2024';site='Unknown'
    selected=support & reviewed & valid
    stronger=selected & ambient
    selected_strict=selected & strict
    d=c[['sample','patient_id','native_copykat_call','provisional_lineage','full_scDblFinder_class']].copy()
    d['source']=source;d['site']=site;d['method_reference_pass']=valid
    d['reviewed_identity_pool']=reviewed;d['native_diploid_fib_singlet_before_reference_gate']=support & reviewed
    d['reference_validated_reviewed_candidate']=selected
    d['ambient_supported_reference_validated_candidate']=stronger
    d['legacy_strict_supported_reference_validated_candidate']=selected_strict
    d['formal_admission']=False
    cell_rows.append(d)
    rows.append({'sample':sm,'patient_id':str(c.patient_id.iloc[0]),'source':source,'site':site,'full_cells':len(c),'reference_pass':valid,'reviewed_diploid_fib_singlets_before_reference_gate':int((support&reviewed).sum()),'reference_validated_reviewed_candidates':int(selected.sum()),'ambient_supported_candidates':int(stronger.sum()),'legacy_strict_supported_candidates':int(selected_strict.sum()),'formal_admission':False})
    hashes[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
s=pd.DataFrame(rows)
assert s[s.source.eq('NCC_Zheng2024')].patient_id.nunique()==8
assert s.set_index('sample').loc['NCC__BS4','patient_id']==s.set_index('sample').loc['NCC__BS14','patient_id']
assert s[s.source.eq('GSE234187')].patient_id.nunique()==1
patients=s.groupby(['source','patient_id'],as_index=False).agg(libraries=('sample',lambda x:';'.join(sorted(x))),library_count=('sample','size'),reference_pass_libraries=('reference_pass','sum'),full_cells=('full_cells','sum'),reviewed_candidates=('reference_validated_reviewed_candidates','sum'),ambient_supported_candidates=('ambient_supported_candidates','sum'),legacy_strict_supported_candidates=('legacy_strict_supported_candidates','sum'))
cells=pd.concat(cell_rows);assert cells.index.is_unique
cells.to_csv(out/'cell_count_evidence.tsv.gz',sep='\t',encoding='utf-8',index_label='cell_id')
s.to_csv(out/'library_counts.tsv',sep='\t',encoding='utf-8',index=False)
patients.to_csv(out/'patient_counts.tsv',sep='\t',encoding='utf-8',index=False)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'patient_count_verification_complete_admission_pending','libraries':len(s),'full_cells':int(s.full_cells.sum()),'source_counts':{},'formal_new_patients':0,'scientific_admission_complete':False,'cross_study_independent_patient_total_established':False,'warnings':['Patient4 libraries BS4 and BS14 count once.','Reference pass is method QC, not patient admission.','Expanded reviewed candidates include uncertain identities; ambient support is additional evidence, not automatic formal CAF confirmation.','Lung reference insufficiency and BS9 reference classification failure remain distinct.','NCC anatomy is unverified; do not add these patients to Bone/Lung inference.','No patient-level inferential tests performed.'],'input_sha256':hashes}
for source,g in patients.groupby('source'):
    summary['source_counts'][source]={'author_resolved_patients':len(g),'patients_with_reference_pass':int(g.reference_pass_libraries.gt(0).sum()),'patients_with_reviewed_candidates':int(g.reviewed_candidates.gt(0).sum()),'patients_with_ambient_supported_candidates':int(g.ambient_supported_candidates.gt(0).sum()),'patients_with_legacy_strict_supported_candidates':int(g.legacy_strict_supported_candidates.gt(0).sum())}
(out/'VERIFIED_COUNTS.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
assert len(pd.read_csv(out/'patient_counts.tsv',sep='\t',encoding='utf-8'))==len(patients)
print(json.dumps(summary,indent=2),flush=True)
