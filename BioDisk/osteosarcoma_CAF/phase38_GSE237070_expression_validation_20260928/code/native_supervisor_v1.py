from pathlib import Path
import subprocess,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase38_GSE237070_expression_validation_20260928')
p=subprocess.Popen(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python',"-u",'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase38_GSE237070_expression_validation_20260928/code/external_expression_validation_v1.py'],cwd=O)
(O/"analysis.pid").write_text(str(p.pid)+"\n")
r=p.wait()
(O/"analysis_native.exit_code").write_text(str(r)+"\n")
sys.exit(r)
