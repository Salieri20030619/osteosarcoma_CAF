from pathlib import Path
import sys
sys.path=[q for q in sys.path if q and Path(q).resolve()!=Path(__file__).parent.resolve()]
from reportlab.pdfbase.ttfonts import TTFont
import pymupdf,json
p=Path(r'C:\Users\31937\Downloads\CAF_F1F5第一轮交付_20261004')
f=TTFont('check',r'C:\Windows\Fonts\msyh.ttc')
s=(p/'CAF_F1F5报告纯文本.txt').read_text(encoding='utf-8');missing=sorted(set(c for c in s if ord(c)>32 and ord(c) not in f.face.charToGlyph));print('MISSING_GLYPHS',missing)
d=pymupdf.open(p/'CAF_F1F5第一轮全转录组与GSEA报告.pdf');bounds=[]
for i,page in enumerate(d):
 blocks=page.get_text('blocks');bad=[b[:4] for b in blocks if b[0]<45 or b[2]>550 or b[1]<35 or b[3]>825];bounds.append({'page':i+1,'text_characters':len(page.get_text()),'out_of_bounds':bad})
assert len(d)==17 and not missing and all(not z['out_of_bounds'] for z in bounds)
(Path(r'C:\Users\31937\CAF_F1F5_round1_20261004')/'PDF_layout_QA.json').write_text(json.dumps({'pages':len(d),'missing_glyphs':missing,'bounds':bounds},indent=2));print('PDF_LAYOUT_PASS',len(d))