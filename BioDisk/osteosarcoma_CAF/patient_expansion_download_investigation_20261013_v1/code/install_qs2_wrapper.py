from pathlib import Path
import subprocess,os
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
env=os.environ.copy();env.update(http_proxy="http://127.0.0.1:18797",https_proxy="http://127.0.0.1:18797",MAKEFLAGS="-j2")
r=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript","--vanilla",str(O/"code/install_qs2.R")],env=env)
(O/"logs/qs2_install.exit_code").write_text(str(r.returncode)+"\n")
