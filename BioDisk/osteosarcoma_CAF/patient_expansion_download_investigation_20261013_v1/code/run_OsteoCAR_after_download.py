from pathlib import Path
import subprocess
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
for job in ["OsteoCAR_download","qs2_install_fixed"]:
 pid=(O/"resources"/(job+".pid")).read_text().strip();subprocess.run(["tail","--pid="+pid,"-f","/dev/null"],check=True)
assert (O/"logs/OsteoCAR_download.exit_code").read_text().strip()=="0", "Download has not passed"
assert (O/"resources/OSTEOCAR_DOWNLOAD_COMPLETE.json").exists()
assert (O/"logs/qs2_install_fixed.exit_code").read_text().strip()=="0", "Reader installation has not passed"
r=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript","--vanilla",str(O/"code/extract_OsteoCAR.R")]);(O/"logs/OsteoCAR_extraction.exit_code").write_text(str(r.returncode)+"\n")
