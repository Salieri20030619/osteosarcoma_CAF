from pathlib import Path
from zipfile import ZipFile
from lxml import etree
import json,hashlib,shutil
src=Path(r'D:\电脑管家迁移文件\xwechat_files\wxid_913cjq0lvmfd32_3e57\msg\file\2026-10\20261005骨肉瘤CAF干实验生信分析避雷清单.docx')
out=Path(r'C:\Users\31937\CAF_client_checklist_review_20261005');out.mkdir(exist_ok=True)
shutil.copy2(src,out/'client_checklist_original.docx')
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with ZipFile(src) as z:
 root=etree.fromstring(z.read('word/document.xml'));lines=[]
 for i,el in enumerate(root.find('w:body',ns)):
  text=''.join(el.xpath('.//w:t/text()',namespaces=ns))
  if el.tag.endswith('}tbl'):
   rows=[]
   for tr in el.findall('w:tr',ns):rows.append(' | '.join(''.join(c.xpath('.//w:t/text()',namespaces=ns)) for c in tr.findall('w:tc',ns)))
   text='\n'.join(rows)
  if text:lines.append(f'[{i}] {text}')
 extra={}
 for name in z.namelist():
  if name in ['word/comments.xml','word/footnotes.xml','word/endnotes.xml']:
   extra[name]='\n'.join(etree.fromstring(z.read(name)).xpath('//w:t/text()',namespaces=ns))
 info={'source':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'bytes':src.stat().st_size,'inline_images':len([n for n in z.namelist() if n.startswith('word/media/')]),'tracked_insertions':len(root.findall('.//w:ins',ns)),'tracked_deletions':len(root.findall('.//w:del',ns)),'extra_parts':extra}
(out/'source_receipt.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8')
(out/'extracted_text.txt').write_text('\n'.join(lines),encoding='utf-8')
print(json.dumps(info,ensure_ascii=False));print('\n'.join(lines))
