from pathlib import Path
import json
import pandas as pd,numpy as np
O=Path(__file__).resolve().parent;T=O/'tables';rd=lambda n:pd.read_csv(T/n,sep='\t')
c=rd('01_clinical_sample_audit.tsv');cells=rd('02_customer_panel_cell_evidence.tsv.gz');s=rd('03_surface_screen_sample_counts.tsv.gz');rank=rd('03_surface_candidate_ranking.tsv')
assert int((O/'exit_code').read_text())==0 and len(c)==40 and c['sample'].is_unique
assert c.unified_strict.sum()==3976 and c.retained_cells.sum()==255241
assert cells.cell_id.is_unique and cells.unified_strict.sum()==3976 and cells.customer_myCAF_target.sum()==1926
for k in ['target_in_gate','strict_CAF_in_gate','mural_in_gate','CNV_risk_in_gate','other_or_unresolved_in_gate']:assert (s[k]<=s.gate_cells).all()
assert (s.target_in_gate<=s.target_cells).all()
sizes=s.groupby('sample').all_QC_cells.first();assert sizes.sum()==255241
assert sizes.reindex(c['sample']).to_numpy().tolist()==c.retained_cells.tolist()
target=cells.groupby('sample').customer_myCAF_target.sum()
assert s.groupby('sample').target_cells.first().reindex(c['sample'],fill_value=0).tolist()==target.reindex(c['sample'],fill_value=0).tolist()
my=rd('02_panel_sample_counts.tsv');my=my[(my['set']=='strict')&(my.panel=='myCAF_customer')].set_index('sample')
complete=c[['sample','dataset','patient_verified','clinical_group','site_verified','treatment_evidence','unified_strict']].copy();complete['myCAF_supported2']=complete['sample'].map(my.supported2).fillna(0).astype(int);complete['myCAF_supported3']=complete['sample'].map(my.supported3).fillna(0).astype(int)
complete.to_csv(T/'02_all40_customer_myCAF_availability.tsv',sep='\t',index=False)
qa={'all_checks_passed':True,'clinical_rows':40,'strict_CAF':3976,'myCAF_operational_candidates':1926,'all_QC_screen_denominator':255241,'RNA_fractions_not_FACS_purity':True,'formal_A_B1_B2_DE_not_run':True}
(O/'REVIEW_QA.json').write_text(json.dumps(qa,indent=2));print(complete.groupby('clinical_group')[['unified_strict','myCAF_supported2','myCAF_supported3']].sum().to_string())
print(rank[(rank.scope=='B2_target_records')&rank.markers.isin(['FAP+THY1','PDGFRA+FAP','FAP'])].to_string(index=False))
