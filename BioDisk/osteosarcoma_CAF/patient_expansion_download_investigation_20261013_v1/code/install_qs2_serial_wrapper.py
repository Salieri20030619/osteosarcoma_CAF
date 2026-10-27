from pathlib import Path
import subprocess,os
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
E="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2"
env=os.environ.copy();env.update(PATH=E+"/bin:"+env.get("PATH",""),CONDA_PREFIX=E,LD_LIBRARY_PATH=E+"/lib:"+env.get("LD_LIBRARY_PATH",""),MAKEFLAGS="-j1",CMAKE_BUILD_PARALLEL_LEVEL="2")
r=subprocess.run([E+"/bin/Rscript","--vanilla",str(O/"code/install_qs2_local_sources.R")],env=env);(O/"logs/qs2_install_fixed.exit_code").write_text(str(r.returncode)+"\n")
