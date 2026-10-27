doc=Document(original)
for child in list(doc._element.body):
 if child.tag!=qn('w:sectPr'):doc._element.body.remove(child)
sec=doc.sections[0];sec.top_margin=sec.bottom_margin=Inches(.64);sec.left_margin=sec.right_margin=Inches(.68)
sec.footer_distance=Inches(.27)
for name in ['Normal','Title','Subtitle','Heading 1','Heading 2','Caption']:
 st=doc.styles[name];st.font.name='Calibri';st.font.color.rgb=RGBColor(0,0,0);st.font.italic=False
 st._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体' if name in ['Normal','Caption'] else '微软雅黑')
st=doc.styles['Normal'];st.font.size=Pt(10.5);st.paragraph_format.line_spacing=1.13;st.paragraph_format.space_after=Pt(6)
for name,size in [('Title',20),('Subtitle',11),('Heading 1',14),('Heading 2',11.5)]:
 st=doc.styles[name];st.font.size=Pt(size);st.font.bold=name!='Subtitle';st.paragraph_format.space_before=Pt(8);st.paragraph_format.space_after=Pt(6)
doc.styles['Title'].paragraph_format.space_before=Pt(0)
doc.styles['Caption'].font.size=Pt(8.8);doc.styles['Caption'].paragraph_format.line_spacing=1.08;doc.styles['Caption'].paragraph_format.space_after=Pt(6)
def p(text,style=None):return doc.add_paragraph(text,style)
def h(text,level=1):return doc.add_heading(text,level)
def page():doc.add_page_break()
def fig(name,caption,width=6.9):
 par=doc.add_paragraph();par.alignment=WD_ALIGN_PARAGRAPH.CENTER;par.paragraph_format.space_after=Pt(3);par.paragraph_format.keep_with_next=True
 r=par.add_run();shape=r.add_picture(str(F/(name+'.png')),width=Inches(width));shape._inline.docPr.set('descr',caption)
 p(caption,'Caption')
def table(headers,rows,widths,font=9.2):
 t=doc.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
 for col,w in zip(t.columns,widths):col.width=Inches(w)
 for cell,txt in zip(t.rows[0].cells,headers):cell.text=str(txt)
 for row in rows:
  for cell,txt in zip(t.add_row().cells,row):cell.text=str(txt)
 bd=OxmlElement('w:tblBorders')
 for edge in ['top','left','bottom','right','insideH','insideV']:
  x=OxmlElement('w:'+edge);x.set(qn('w:val'),'single');x.set(qn('w:sz'),'4');x.set(qn('w:color'),'D9D9D9');bd.append(x)
 t._tbl.tblPr.append(bd)
 for i,row in enumerate(t.rows):
  tr=row._tr.get_or_add_trPr();tr.append(OxmlElement('w:cantSplit'))
  if i==0:tr.append(OxmlElement('w:tblHeader'))
  for j,cell in enumerate(row.cells):
   cell.width=Inches(widths[j]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER;pr=cell._tc.get_or_add_tcPr()
   shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'284B63' if i==0 else ('F1F5F8' if i%2 else 'FFFFFF'));pr.append(shade)
   m=OxmlElement('w:tcMar')
   for side in ['top','bottom','left','right']:
    x=OxmlElement('w:'+side);x.set(qn('w:w'),'75');x.set(qn('w:type'),'dxa');m.append(x)
   pr.append(m)
   for par in cell.paragraphs:
    par.paragraph_format.space_after=Pt(0);par.paragraph_format.line_spacing=1.05
    if j==0 or len(str(cell.text))<9:par.alignment=WD_ALIGN_PARAGRAPH.CENTER
    for r in par.runs:r.font.size=Pt(font);r.font.bold=i==0;r.font.color.rgb=RGBColor.from_string('FFFFFF' if i==0 else '000000')
 return t
