from pathlib import Path
import subprocess
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
r=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python",str(O/"code/investigate_S23.py")]);(O/"logs/S23_investigation.exit_code").write_text(str(r.returncode)+"\n")
