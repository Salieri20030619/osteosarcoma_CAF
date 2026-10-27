from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile,shutil
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase112_allcell_3D_landscape_20261021'
receipt=json.loads((O/'resources/FINAL_STATUS.json').read_text(encoding='utf-8'))
archive=O/'delivery/CAF_Fig1b_1d_600dpi.zip'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==receipt['archive_sha256']
with zipfile.ZipFile(archive) as z:assert z.testzip() is None and len(z.namelist())==receipt['zip_members']
receipt.update(server_archive=str(archive),server_archive_hash_verified=True,server_delivery_checked_at=datetime.now(timezone.utc).isoformat())
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
dump(O/'FINAL_STATUS.json',receipt);dump(O/'STATUS.json',receipt);dump(R/'CAF_REFERENCE_FIGURES_CURRENT.json',receipt)
p=R/'CAF_RESUME_STATUS_CURRENT.json';shutil.copy2(p,O/'resources/CAF_RESUME_STATUS_before_phase112.json')
current=json.loads(p.read_text(encoding='utf-8'))
current.update(authoritative_current_file=str(R/'CAF_REFERENCE_FIGURES_CURRENT.json'),phase=O.name,current_phase=str(O),status='allcell_3D_and_composition_figures_complete',current_state='completed_verified',task_status='completed_verified',active_task='Full-cell 3D UMAP and cell-type composition figure batch complete',updated_at=receipt['server_delivery_checked_at'],next_action='Continue authorized paper-figure adaptations from preserved queue; prior Fig1m boxplot v2 accepted. No original analysis overwritten.')
dump(p,current)
print(json.dumps(dict(status='PASS',archive_sha256=receipt['archive_sha256'],static_figures=receipt['static_figures'],interactive_html=receipt['interactive_html'],progress_updated=True)))
