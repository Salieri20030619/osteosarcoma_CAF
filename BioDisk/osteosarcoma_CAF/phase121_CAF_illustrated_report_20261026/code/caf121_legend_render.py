from pathlib import Path
import subprocess,json,hashlib,zipfile
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase121_CAF_illustrated_report_20261026');p=O/'delivery/CAF_illustrated_report_legends.pdf';assert hashlib.sha256(p.read_bytes()).hexdigest()=='4c22e02df33f01f93aaa5943eab873cc00cec50398c89d3d6e5f6ac11dd57ead'
d=O/'resources/legend_render';d.mkdir(exist_ok=True)
for n in [3,4,11]:subprocess.run(['pdftoppm','-f',str(n),'-l',str(n),'-singlefile','-r','130','-png',str(p),str(d/('page'+str(n)))],check=True,capture_output=True)
with zipfile.ZipFile(O/'delivery/legend_render.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in d.glob('*.png'):z.write(f,f.name)
print(hashlib.sha256((O/'delivery/legend_render.zip').read_bytes()).hexdigest())
