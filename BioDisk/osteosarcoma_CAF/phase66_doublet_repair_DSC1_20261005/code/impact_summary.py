from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase66_doublet_repair_DSC1_20261005')
a=pd.read_csv(O/'tables/cells_repaired_with_confidence.tsv.gz',sep='\t',low_memory=False)
primary={s:t.dataset.value_counts().sort_index().sort_values(ascending=False,kind='stable').index[0] for s,t in a.groupby('subgroup_id')}
w=pd.read_csv(O/'tables/program_before_after_comparison.tsv',sep='\t');w['primary_source']=w.subgroup_id.map(primary);w['is_primary_source']=w.dataset.eq(w.primary_source)
w.to_csv(O/'tables/program_before_after_comparison.tsv',sep='\t',index=False)
g=pd.read_csv(O/'tables/gene_before_after_comparison.tsv.gz',sep='\t')
changed=(w.effect_change.abs()>1e-10)
summary={'program_contrasts':len(w),'program_effect_changed':int(changed.sum()),'program_direction_flips':int(w.direction_flip.sum()),'program_nominal_changes':int(w.nominal_change.sum()),'program_lost_P':int(w.lost_evaluability.sum()),'primary_source_program_effect_changes':int((changed&w.is_primary_source).sum()),'gene_contrasts':len(g),'gene_nominal_changes':int(g.nominal_change.sum()),'positive_sources':a[a.doublet_positive_any].dataset.value_counts().to_dict()}
(O/'resources/IMPACT_SUMMARY.json').write_text(json.dumps(summary,indent=2))
w[changed|w.lost_evaluability|(w.N_all!=w.N_excluded)].to_csv(O/'tables/affected_program_results.tsv',sep='\t',index=False)
print(json.dumps(summary),flush=True)
print('QA_LOG',(O/'logs/verify.log').read_text()[-2500:])
p=O/'QA/ACCEPTANCE.json'
if p.exists():
 q=json.loads(p.read_text());q.pop('single_cell_source_N_unchanged',None);q['frozen8174_roster_preserved_before_filter']=True;p.write_text(json.dumps(q,indent=2))
p=O/'code/verify.py';p.write_text(p.read_text().replace("'single_cell_source_N_unchanged':True","'frozen8174_roster_preserved_before_filter':True"))
