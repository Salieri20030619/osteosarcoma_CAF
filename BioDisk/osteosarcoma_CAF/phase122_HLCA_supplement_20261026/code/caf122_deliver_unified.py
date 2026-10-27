from pathlib import Path
import json,hashlib,os,shutil,pymupdf,base64,zlib,zipfile,io
B=Path('C:/Users/31937/caf122_hlca_supplement_20261026');W=B/'unified_full_report';v=json.loads((W/'UNIFIED_FINAL_QA.json').read_text(encoding='utf-8'));assert v['status']=='PASS'
paths=[(W/'full_report_review.pdf',Path('C:/Users/31937/Downloads/骨肺CAF亚群_图文分析与Lung04推荐报告.pdf'),v['report_sha256']),(W/'reviewer_QA_review.pdf',Path('C:/Users/31937/Downloads/骨肺CAF研究_审稿人问答.pdf'),v['QA_sha256'])]
for src,dst,expected in paths:
 assert hashlib.sha256(src.read_bytes()).hexdigest()==expected
 if dst.exists() and dst==paths[1][1] and not (W/'previous_reviewer_QA.pdf').exists():shutil.copyfile(dst,W/'previous_reviewer_QA.pdf')
 tmp=dst.with_name(dst.name+'.tmp');shutil.copyfile(src,tmp)
 try:os.replace(tmp,dst)
 except PermissionError:
  with dst.open('r+b') as f:f.write(src.read_bytes());f.truncate();f.flush();os.fsync(f.fileno())
  tmp.unlink()
 assert hashlib.sha256(dst.read_bytes()).hexdigest()==expected
record=dict(status='COMPLETE',report=str(paths[0][1]),report_sha256=v['report_sha256'],report_pages=17,QA=str(paths[1][1]),QA_sha256=v['QA_sha256'],QA_pages=2,QA_questions=5,main_filename_unchanged=True,normal_reference_page=11,recommendation_page=13,checks=78,passed=78,previous_report_backup=str(W/'previous_17page_report.pdf'))
(W/'UNIFIED_DELIVERY_COMPLETE.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
# Reuse verified identical image pixels already archived on the configured server.
new=pymupdf.open(paths[0][0]);old=pymupdf.open(W/'previous_17page_report.pdf');pixels={}
for pg in old:
 for it in pg.get_images():
  pix=pymupdf.Pixmap(old,it[0]);pixels[hashlib.sha256(pix.samples).hexdigest()]=True
blob=paths[0][0].read_bytes();items=[];seen=set()
for pg in new:
 for it in pg.get_images():
  x=it[0]
  if x in seen:continue
  seen.add(x);raw=new.xref_stream_raw(x);pix=pymupdf.Pixmap(new,x);key=hashlib.sha256(pix.samples).hexdigest();assert key in pixels
  encoded=base64.a85encode(zlib.compress(pix.samples))+b'~>';assert raw==encoded,('stream_encoding',x)
  pos=blob.find(raw);assert pos>=0 and blob.find(raw,pos+1)<0;items.append(dict(start=pos,length=len(raw),pixel_sha256=key,stream_sha256=hashlib.sha256(raw).hexdigest()))
items.sort(key=lambda x:x['start']);cursor=0;manifest=[]
with zipfile.ZipFile(W/'unified_server_transport.zip','w',zipfile.ZIP_DEFLATED) as out:
 for i,it in enumerate(items):
  name=f'literal_{i}.bin';out.writestr(name,blob[cursor:it['start']]);manifest.append(dict(literal=name,**it));cursor=it['start']+it['length']
 out.writestr('tail.bin',blob[cursor:]);out.writestr('manifest.json',json.dumps(dict(report_sha256=v['report_sha256'],items=manifest,tail='tail.bin',source_server_pdf='CAF_full_customer_report_integrated_20261027.pdf',output_server_pdf='CAF_full_customer_report_unified_20261027.pdf')))
 out.write(paths[1][0],'reviewer_QA.pdf');out.write(W/'UNIFIED_DELIVERY_COMPLETE.json','UNIFIED_DELIVERY_COMPLETE.json');out.write(W/'UNIFIED_FINAL_QA.json','UNIFIED_FINAL_QA.json')
print(json.dumps(record,ensure_ascii=False));print('TRANSPORT_BYTES', (W/'unified_server_transport.zip').stat().st_size)
