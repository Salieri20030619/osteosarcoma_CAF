from pathlib import Path
import pymupdf as fitz,json,hashlib,re,shutil,io,os,sys
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.lib import colors
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026/expanded_QA_revision/q9_denominator_fix');B.mkdir(exist_ok=True);target=Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答_扩充版.pdf');backup=B/'previous_expanded_QA.pdf';staged=B/'review.pdf';record=B/'Q9_COMPLETION.json'
if '--commit' in sys.argv:
 r=json.loads(record.read_text(encoding='utf-8'));assert hashlib.sha256(staged.read_bytes()).hexdigest()==r['sha256'];tmp=target.with_name(target.name+'.tmp');shutil.copyfile(staged,tmp)
 try:os.replace(tmp,target)
 except PermissionError:
  with target.open('r+b') as f:f.write(staged.read_bytes());f.truncate();f.flush();os.fsync(f.fileno())
  tmp.unlink()
 assert hashlib.sha256(target.read_bytes()).hexdigest()==r['sha256'];r.update(status='COMPLETE',visual_review='Changed paragraph inspected; no clipping; other pages render identically');record.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(r,ensure_ascii=False));raise SystemExit
if not backup.exists():shutil.copyfile(target,backup)
d=fitz.open(backup);assert len(d)==4;p=d[2];bs=p.get_text('blocks');ids=[19,20,21,22];old=''.join(bs[i][4] for i in ids);assert old.startswith('估计精度') and old.count('且 HLCA 成纤维比例为 100.0%')==1
new=old.replace('\n','').replace('且 HLCA 成纤维比例为 100.0%','且在每位至少5个候选细胞的2位可评估患者中，患者等权平均的成纤维细胞样比例为100.0%')
r=fitz.Rect(bs[ids[0]][:4])
for i in ids[1:]:r|=fitz.Rect(bs[i][:4])
r.x1=550;r.y0-=.7;r.y1+=3
pdfmetrics.registerFont(TTFont('YH','C:/Windows/Fonts/msyh.ttc'));style=ParagraphStyle('body',fontName='YH',fontSize=9.5,leading=14.5,wordWrap='CJK',textColor=colors.HexColor('#243b49'));para=Paragraph(new,style);_,height=para.wrap(r.width,1000);assert height<=r.height,(height,r.height)
p.add_redact_annot(r,fill=None);p.apply_redactions(images=0,graphics=0);buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(p.rect.width,p.rect.height));para.drawOn(c,r.x0,p.rect.height-r.y0-height);c.save();overlay=fitz.open(stream=buf.getvalue(),filetype='pdf');p.show_pdf_page(p.rect,overlay,0);d.save(staged,garbage=4,deflate=True);d.close()
a=fitz.open(backup);b=fitz.open(staged)
def questions(doc):
 t='\n'.join(p.get_text(sort=True) for p in doc);matches=list(re.finditer(r'Q(\d+)\s*｜',t));return {int(m.group(1)):re.sub(r'\s+','',t[m.start():matches[i+1].start() if i+1<len(matches) else len(t)]) for i,m in enumerate(matches)}
aa=questions(a);bb=questions(b);assert set(aa)==set(bb)==set(range(1,12));unchanged=[i for i in range(1,12) if i!=9]
assert all(aa[i]==bb[i] for i in unchanged);assert bb[9]==aa[9].replace('且HLCA成纤维比例为100.0%','且在每位至少5个候选细胞的2位可评估患者中，患者等权平均的成纤维细胞样比例为100.0%')
assert all(a[i].get_pixmap().samples==b[i].get_pixmap().samples for i in [0,1,3]);b[2].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(B/'page3.png')
rec=dict(status='STAGED_PASS',file=str(target),pages=4,changed_question='Q9 only',change='100% denominator clarified as equal-patient mean across2 evaluable patients each with>=5 candidate cells',unchanged_questions_verified=unchanged,other_pages_pixel_identical=True,sha256=hashlib.sha256(staged.read_bytes()).hexdigest(),previous_sha256=hashlib.sha256(backup.read_bytes()).hexdigest(),backup=str(backup));record.write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(rec,ensure_ascii=False))
