from pathlib import Path
import subprocess,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927')
p=subprocess.Popen(["python3","-u",'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927/code/reconstruct_sra_worker_v1.py'],cwd=O)
(O/"reconstruction.pid").write_text(str(p.pid)+"\n")
r=p.wait()
(O/"reconstruction_native.exit_code").write_text(str(r)+"\n")
sys.exit(r)
