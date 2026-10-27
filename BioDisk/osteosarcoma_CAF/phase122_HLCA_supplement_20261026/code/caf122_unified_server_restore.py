from pathlib import Path
import json,hashlib,base64,zlib,zipfile,io,datetime
import sys
sys.path.insert(0,'/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026/resources/pdf_runtime')
import pymupdf as fitz
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';P=O/'delivery'
with zipfile.ZipFile(P/'unified_server_transport.zip') as z:
 assert z.testzip() is None;m=json.loads(z.read('manifest.json'));source=fitz.open(P/m['source_server_pdf']);pixels={}
 for pg in source:
  for it in pg.get_images():
   pix=fitz.Pixmap(source,it[0]);pixels[hashlib.sha256(pix.samples).hexdigest()]=pix.samples
 out=io.BytesIO()
 for it in m['items']:
  out.write(z.read(it['literal']));data=base64.a85encode(zlib.compress(pixels[it['pixel_sha256']]))+b'~>';assert len(data)==it['length'] and hashlib.sha256(data).hexdigest()==it['stream_sha256'];out.write(data)
 out.write(z.read(m['tail']));pdf=out.getvalue();assert hashlib.sha256(pdf).hexdigest()==m['report_sha256'];final=P/m['output_server_pdf'];final.write_bytes(pdf)
 rec=json.loads(z.read('UNIFIED_DELIVERY_COMPLETE.json'));q=z.read('reviewer_QA.pdf');assert hashlib.sha256(q).hexdigest()==rec['QA_sha256'];qp=P/'CAF_reviewer_QA_20261027.pdf';qp.write_bytes(q);assert len(fitz.open(final))==17 and len(fitz.open(qp))==2
 (P/'UNIFIED_FINAL_QA.json').write_bytes(z.read('UNIFIED_FINAL_QA.json'));rec.update(server_report=str(final),server_QA=str(qp),server_hash_verified=True,restored_identical_image_streams=len(m['items']),server_verified_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
 (P/'UNIFIED_CUSTOMER_DELIVERY_COMPLETE.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8');(R/'LATEST_CAF_INTEGRATED_CUSTOMER_REPORT.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8')
 for pointer in [O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:
  state=json.loads(pointer.read_text(encoding='utf-8'));state['integrated_customer_report']=dict(path=str(final),sha256=rec['report_sha256'],pages=17,normal_reference_page=11,reviewer_QA=str(qp),reviewer_QA_pages=2);pointer.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(dict(status='PASS',report_sha256=rec['report_sha256'],QA_sha256=rec['QA_sha256'],reused_verified_images=len(m['items']))))
