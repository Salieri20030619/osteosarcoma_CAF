from pathlib import Path
import subprocess,sys
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
T=O/'interim_GSE21257';T.mkdir(exist_ok=True)
for name in ['sources','resources']:
 p=T/name
 if not p.exists():p.symlink_to(O/name,target_is_directory=True)
for name in ['tables','objects','QA','figures']:(T/name).mkdir(exist_ok=True)
s=(O/'code/analyze_response.py').read_text()
s=s.replace("O=B/'phase64_response_followup_20261005'","O=B/'phase64_response_followup_20261005/interim_GSE21257'")
s=s.replace("[('GSE39055','GPL14951'),('GSE21257','GPL10295')]","[('GSE21257','GPL10295')]")
s=s.replace("zip(axs,['GSE39055','GSE21257'])","zip(axs,['GSE21257'])")
p=O/'code/analyze_response_interim21257.py';p.write_text(s)
with (O/'logs/interim21257.log').open('w') as log:
 r=subprocess.run([sys.executable,str(p)],stdout=log,stderr=subprocess.STDOUT)
print('exit',r.returncode,flush=True)
print((O/'logs/interim21257.log').read_text()[-15000:])
sys.exit(r.returncode)
