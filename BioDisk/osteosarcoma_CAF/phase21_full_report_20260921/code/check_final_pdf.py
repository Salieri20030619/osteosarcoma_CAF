from pathlib import Path
import hashlib,json
from PIL import Image,ImageChops
from pypdf import PdfReader
O=Path('C:/codex_test/CAF_full_report_20260921')
p=next((O/'delivery').glob('*.pdf'));r=PdfReader(str(p))
assert len(r.pages)==15,len(r.pages)
same=[]
for n in range(1,14):
 a=Image.open(O/f'qa/fixed-{n:02}.png').convert('RGB')
 b=Image.open(O/f'qa/final-{n:02}.png').convert('RGB')
 assert a.size==b.size and ImageChops.difference(a,b).getbbox() is None,n
 same.append(n)
t='\n'.join(x.extract_text() for x in r.pages)
assert '305,054' in t and '6,566' in t and '留出' in t
assert '\ufffd' not in t
qa={'pages':len(r.pages),'prior_visually_reviewed_pages_unchanged':same,'pages_for_final_visual_review':[14,15],'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
(O/'qa/final_integrity.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa))
