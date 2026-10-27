from pathlib import Path
import subprocess,os
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1");E="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2"
env=os.environ.copy();env.update(PATH=E+"/bin:"+env.get("PATH",""),LD_LIBRARY_PATH=E+"/lib:"+env.get("LD_LIBRARY_PATH",""))
r=subprocess.run([E+"/bin/Rscript","--vanilla",str(O/"code/resolve_primary_gene_aliases.R")],env=env);(O/"logs/primary_alias_resolution.exit_code").write_text(str(r.returncode)+"\n")
