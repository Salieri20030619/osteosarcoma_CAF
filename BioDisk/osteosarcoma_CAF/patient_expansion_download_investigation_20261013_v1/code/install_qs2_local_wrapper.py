from pathlib import Path
import subprocess,os
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
env=os.environ.copy();env["MAKEFLAGS"]="-j2"
r=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript","--vanilla",str(O/"code/install_qs2_local_sources.R")],env=env);(O/"logs/qs2_install_local.exit_code").write_text(str(r.returncode)+"\n")
