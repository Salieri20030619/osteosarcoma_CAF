from pathlib import Path
import os,signal,subprocess
target=b'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026/resources/pdf_runtime';found=[]
for d in Path('/proc').iterdir():
 if not d.name.isdigit():continue
 try:a=(d/'cmdline').read_bytes().split(b'\0')
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if target in a and b'pip' in a and b'install' in a:found.append(int(d.name))
assert len(found)<=1,found
for pid in found:os.kill(pid,signal.SIGTERM)
print('STOPPED_OWN_STALLED_MIRROR_INSTALL',found,flush=True)
subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python','-m','pip','install','--no-deps','--timeout','12','--retries','0','--index-url','https://pypi.org/simple','--target',target.decode(),'PyMuPDF'],timeout=100,check=True)
