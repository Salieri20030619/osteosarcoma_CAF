from pathlib import Path
import subprocess,json,time,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927')
p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',"-u",'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase37_requested_source_search_20260927/code/sra_download_worker_v2.py'],cwd=O)
(O/"download.pid").write_text(str(p.pid)+"\n")
r=p.wait()
(O/"download_native.exit_code").write_text(str(r)+"\n")
(O/"DOWNLOAD_NATIVE_COMPLETION.json").write_text(json.dumps({"worker_pid":p.pid,"exit_code":r,"finished_at":time.time()}))
sys.exit(r)
