from pathlib import Path
import subprocess,sys,time,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927')
p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',"-u",'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927/code/prepare_reference_worker_v1.py'],cwd=O)
(O/"reference.pid").write_text(str(p.pid)+"\n")
r=p.wait()
(O/"reference_native.exit_code").write_text(str(r)+"\n")
sys.exit(r)
