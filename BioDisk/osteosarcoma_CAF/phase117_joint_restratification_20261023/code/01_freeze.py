from caf117_common import *
import importlib.metadata,shutil
event('01_freeze','started')
assert (O/'PLAN_ZH.txt').exists()
primary=[B/'objects/log1p_cp10k_shared.npz',B/'tables/cell_registry.tsv.gz',B/'tables/shared_genes.tsv',B/'tables/cell_function_scores.tsv.gz',B/'tables/program_coverage.tsv',B/'resources/PROGRAMS.json',R/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',R/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv',R/'phase62_F1F5_internal_20261005/resources/libraries_canonical.json',R/'phase116_original33_analogs_20261023/tables/inherited_pathway_cell_scores.tsv.gz',R/'phase116_original33_analogs_20261023/tables/inherited_pathway_coverage.tsv']
for h in json.loads((B/'resources/INPUT_HASHES.json').read_text()):assert sha(h['path'])==h['sha256'],h['path']
m=pd.read_csv(primary[1],sep='\t',low_memory=False);old=pd.read_csv(primary[6],sep='\t',low_memory=False).set_index('cell_id');assert len(m)==10510 and m.cell_id.is_unique and (~m.doublet).sum()==10274
assert set(m.loc[m.cohort.eq('baseline'),'cell_id'])==set(old.index)
assert m[m.source.eq('NCC')].donor.nunique()==7 and m[m.source.eq('NCI')].donor.nunique()==17
m['original_class']=m.class_label.fillna('Added_unassigned');m['original_state']=m.state_label.fillna('Added_unassigned');m['old_lineage']=m.cell_id.map(old.review_lineage_code).fillna('Added_unverified')
m['old_lineage_confidence']=m.cell_id.map(old.lineage_confidence).fillna('not_assessed');m['old_mixed']=m.cell_id.map(old.mixed_program_flag).fillna(False).astype(bool)
m['truth_lineage']=np.where(m.old_lineage.isin(['Fibro','Osteo','Mural'])&m.old_lineage_confidence.isin(['高','中'])&~m.old_mixed,m.old_lineage,'Unassessed')
save(m,'frozen_registry.tsv.gz');save(m.groupby(['source','cohort']).agg(cells=('cell_id','size'),source_patients=('donor','nunique'),doublets=('doublet','sum')).reset_index(),'frozen_source_counts.tsv')
shutil.copy2(primary[7],O/'tables/original33_dictionary.tsv')
dump('resources/INPUT_HASHES.json',[dict(path=str(p),sha256=sha(p)) for p in primary]);dump('resources/ENVIRONMENT.json',{k:importlib.metadata.version(k) for k in ['numpy','pandas','scanpy','scipy','anndata','harmonypy','scikit-learn','leidenalg']})
protocol=dict(version=1,frozen_at=now(),plan_sha256=sha(O/'PLAN_ZH.txt'),primary_cells=10274,all_cells=10510,original_cells=8174,new_cells=2336,seed=20261023,threads=4,primary_resolution=.5,sensitivity_resolutions=[.3,.8],pilot_theta=[1,2],pilot_macro_purity_loss_max=.05,pilot_each_lineage_loss_max=.10,unknown_labels_preserved=True,figure_queue='paused_by_user',flow=['01_freeze','02_structure','03_correspondence','04_function','05_holdout','06_acceptance_report'])
dump('PROTOCOL_V1.json',protocol);dump('FREEZE_COMPLETE.json',dict(status='PASS',cells=len(m),old_states=old.subgroup_id.nunique(),input_files=len(primary),upstream38_hashes_checked=True,plan_sha256=protocol['plan_sha256']))
event('01_freeze','completed',cells=len(m),fit_cells=int((~m.doublet).sum()))
