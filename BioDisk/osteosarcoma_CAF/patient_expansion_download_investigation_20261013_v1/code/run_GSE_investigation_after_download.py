from pathlib import Path
import subprocess,json
O=Path("/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/patient_expansion_download_investigation_20261013_v1")
pid=(O/"resources/GSE234187_download.pid").read_text().strip()
subprocess.run(["tail","--pid="+pid,"-f","/dev/null"],check=True)
assert (O/"logs/GSE234187_download.exit_code").read_text().strip()=="0"
assert (O/"resources/GSE234187_DOWNLOAD_COMPLETE.json").exists()
r=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python",str(O/"code/investigate_GSE234187.py")]);(O/"logs/GSE234187_investigation.exit_code").write_text(str(r.returncode)+"\n")
