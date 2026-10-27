from pathlib import Path
import subprocess
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
pid=(O/"resources/OsteoCAR_extraction.pid").read_text().strip();subprocess.run(["tail","--pid="+pid,"-f","/dev/null"])
assert (O/"logs/OsteoCAR_extraction.exit_code").read_text().strip()=="0"
r=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python",str(O/"code/investigate_all_human_objects.py")]);(O/"logs/all_human_investigation.exit_code").write_text(str(r.returncode)+"\n")
