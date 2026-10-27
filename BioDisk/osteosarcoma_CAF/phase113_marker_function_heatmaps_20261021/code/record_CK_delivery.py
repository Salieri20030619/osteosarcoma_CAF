from pathlib import Path
import json,hashlib,zipfile,shutil
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase113_marker_function_heatmaps_20261021'
r=json.loads((O/'QA/FINAL_STATUS.json').read_text(encoding='utf-8'));archive=O/'delivery/CAF_Fig1c_Fig1k_600dpi.zip'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==r['archive_sha256']
with zipfile.ZipFile(archive) as z:assert z.testzip() is None and len(z.namelist())==r['archive_members']
paths=['phase27_identity_review_20260926/resources/lineage_panels.json','phase66_doublet_repair_DSC1_20261005/resources/PROGRAMS.json','phase112_allcell_3D_landscape_20261021/tables/all_cell_metadata.tsv.gz','phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv','phase111_reference_figure_adaptation_20261021/tables/Fig1i_cell_values.tsv.gz','phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz','phase53_de_novo_reannotation_20261002/tables/SCPCP000017/input_manifest.tsv','phase53_de_novo_reannotation_20261002/tables/SCPCP000023/input_manifest.tsv']
(O/'resources/INPUT_HASHES.json').write_text(json.dumps([dict(path=str(R/p),sha256=hashlib.sha256((R/p).read_bytes()).hexdigest()) for p in paths],indent=2),encoding='utf-8')
r.update(server_archive=str(archive),server_archive_verified=True)
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
dump(O/'FINAL_STATUS.json',r);dump(O/'STATUS.json',r)
current=R/'CAF_REFERENCE_FIGURES_CURRENT.json'
if current.exists():shutil.copy2(current,O/'resources/FIGURES_CURRENT_before_phase113.json')
dump(current,r)
resume=R/'CAF_RESUME_STATUS_CURRENT.json';shutil.copy2(resume,O/'resources/RESUME_before_phase113.json');d=json.loads(resume.read_text(encoding='utf-8'))
d.update(authoritative_current_file=str(current),phase=O.name,current_phase=str(O),status='Fig1c_Fig1k_batch_complete',current_state='completed_verified',task_status='completed_verified',active_task='Full-cell marker dotplots and CAF candidate functional heatmaps complete',updated_at=r['completed_at'],next_action='Continue authorized reference-figure queue; retain frozen labels, source scopes and six candidates.')
dump(resume,d);print(json.dumps(dict(status='PASS',archive_sha256=r['archive_sha256'],figures=6,progress_updated=True)))
