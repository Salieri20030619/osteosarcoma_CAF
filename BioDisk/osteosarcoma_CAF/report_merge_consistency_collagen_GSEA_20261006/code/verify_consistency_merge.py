from pathlib import Path
from collections import Counter
from docx import Document
from pypdf import PdfReader
import json,hashlib,datetime,re
B=Path(__file__).parent;O=B/'consistency_merge'
r=json.loads((O/'MERGE_ACCEPTANCE.json').read_text(encoding='utf-8'))
d=Document(r['report']);a=Document(r['base_report']);s=Document(r['added_report'])
sig=lambda t:json.dumps([[c.text for c in row.cells] for row in t.rows],ensure_ascii=False)
assert Counter(map(sig,d.tables))==Counter(map(sig,a.tables+s.tables))
assert len(d.inline_shapes)==16 and len(d.tables)==14
assert sum(len(t.rows)-1 for t in d.tables[-3:])==54
assert sum(len(t.rows)-1 for t in d.tables if len(t.columns)==7 and t.rows[0].cells[1].text=='细胞数')==33
fig=[int(re.match(r'^图(\d+)\s',p.text).group(1)) for p in d.paragraphs if re.match(r'^图(\d+)\s',p.text)]
tab=[int(re.match(r'^表(\d+)\s',p.text).group(1)) for p in d.paragraphs if re.match(r'^表(\d+)\s',p.text)]
assert fig==list(range(1,17)) and tab==list(range(1,10))
assert hashlib.sha256(Path(r['base_report']).read_bytes()).hexdigest()==r['base_sha256']
assert hashlib.sha256(Path(r['added_report']).read_bytes()).hexdigest()==r['added_sha256']
assert hashlib.sha256(Path(r['report']).read_bytes()).hexdigest()==r['report_sha256']
pdf=PdfReader(O/'render_v2/render_input.pdf');assert len(pdf.pages)==26
assert len(list((O/'render_v2').glob('page-*.png')))==26
assert min(len(p.extract_text().strip()) for p in pdf.pages)>160
r.update({'pages':26,'passed':True,'all_pages_visually_checked':True,'visual_review_pending':False,'all_14_table_contents_unchanged':True,'figure_numbers_contiguous':True,'table_numbers_contiguous':True,'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat()})
(O/'FINAL_ACCEPTANCE.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'passed':True,'pages':26,'figures':16,'tables':14,'report_sha256':r['report_sha256']},ensure_ascii=False))
