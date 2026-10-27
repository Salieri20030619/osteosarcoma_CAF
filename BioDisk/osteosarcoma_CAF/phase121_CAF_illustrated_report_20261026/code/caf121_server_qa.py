from pathlib import Path
import subprocess,json,hashlib,zipfile,shutil
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase121_CAF_illustrated_report_20261026';p=O/'delivery/CAF_illustrated_report.pdf';expected=json.loads((O/'resources/PDF_BUILD_RECEIPT.json').read_text())['sha256']
assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
render=O/'resources/pdf_render';render.mkdir(exist_ok=True)
r=subprocess.run(['pdftoppm','-r','110','-png',str(p),str(render/'page')],capture_output=True)
(O/'logs/pdf_render.stderr.txt').write_bytes(r.stderr);assert r.returncode==0
pages=sorted(render.glob('page-*.png'));assert len(pages)==19
with zipfile.ZipFile(O/'delivery/render_review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in pages:z.write(f,f.name)
shutil.copyfile(R/'caf121_figures.py',O/'code/01_figures.py')
rec=dict(status='PASS_RENDER',pdf_sha256=expected,pages=19,pdftoppm_exit=0,render_archive_sha256=hashlib.sha256((O/'delivery/render_review.zip').read_bytes()).hexdigest())
(O/'resources/RENDER_RECEIPT.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec),flush=True)
