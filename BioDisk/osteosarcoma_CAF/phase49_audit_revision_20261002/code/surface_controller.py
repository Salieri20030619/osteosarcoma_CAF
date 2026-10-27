from pathlib import Path
import subprocess,sys
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase49_audit_revision_20261002')
with (R/'logs/surface_revision.log').open('w') as log: rc=subprocess.call([sys.executable,'-u',str(R/'code/revise_surface.py')],stdout=log,stderr=subprocess.STDOUT)
(R/'surface_revision.exit_code').write_text(str(rc))
raise SystemExit(rc)
