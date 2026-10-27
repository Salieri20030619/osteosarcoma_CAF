from pathlib import Path
import json,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012'
old=pd.read_csv(P/'phase27_identity_review_20260926/tables/SCPCS000525_review.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
base=P/'phase26_client_reanalysis_20260922/resource_recovery_SCPCS000525_20260925'
print('RECOVERY_FILES',[p.name for p in base.iterdir()])
new=pd.read_csv(base/'tables/SCPCS000525_identity_evidence.tsv.gz',sep='\t',low_memory=False).set_index('cell_id')
assert old.index.is_unique and new.index.is_unique and set(old.index)==set(new.index)
new=new.loc[old.index];fib=old.provisional_lineage.eq('Fibroblast');db=old.known_doublet.astype(str).str.lower().eq('true');gate=fib&~db&new.CNV_evidence.isin(['concordant_no_alteration','single_method_no_alteration'])
out=dict(sample='SCPCS000525',patient='SARC0088',old_cells=len(old),old_CAF_candidates=int(old.CAF_candidate.sum()),recovery_rule_candidates=int(gate.sum()),current_patient_already_represented_by='SCPCS000526',new_patient_added=False,new_CNV_evidence=new.CNV_evidence.value_counts().to_dict(),recoverable_provisional_candidates_need_final_lineage_review=True)
(O/'resources/LATE_S525_RECOVERY_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
for f in base.glob('**/*QA.json'):
    if f.stat().st_size<10000:print(f.relative_to(base),f.read_text()[:3000])
