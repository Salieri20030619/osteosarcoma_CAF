#!/usr/bin/env python3
import argparse,pathlib,zipfile
p=argparse.ArgumentParser();p.add_argument('--results',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
assert (a.results/'SUCCESS.txt').exists(),'Analysis has not completed'
a.output.parent.mkdir(parents=True,exist_ok=True)
temp=a.output.with_suffix('.tmp')
with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED) as z:
 for f in sorted(a.results.rglob('*')):
  if f.is_file() and f.suffix!='.h5ad':z.write(f,str(f.relative_to(a.results)))
 docs=pathlib.Path(__file__).resolve().parent.parent/'docs'
 if docs.exists():
  for f in docs.glob('*.txt'):z.write(f,'docs/'+f.name)
temp.replace(a.output)
print('REVIEW_PACKAGE_COMPLETE',a.output,flush=True)
