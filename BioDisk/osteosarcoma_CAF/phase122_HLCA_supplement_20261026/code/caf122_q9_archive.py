from pathlib import Path
import json,hashlib,difflib,base64,gzip,sys
if '--server' not in sys.argv:
 W=Path('C:/Users/31937/caf122_hlca_supplement_20261026/expanded_QA_revision/q9_denominator_fix');a=(W/'previous_expanded_QA.pdf').read_bytes();b=(W/'review.pdf').read_bytes();aa=a.splitlines(keepends=True);bb=b.splitlines(keepends=True);offset=[0]
 for line in aa:offset.append(offset[-1]+len(line))
 ops=[]
 for tag,i,j,k,l in difflib.SequenceMatcher(None,aa,bb).get_opcodes():
  if tag=='equal':ops.append(['copy',offset[i],offset[j]-offset[i]])
  elif tag in ['insert','replace']:ops.append(['data',base64.b64encode(b''.join(bb[k:l])).decode()])
 payload=dict(source_sha256=hashlib.sha256(a).hexdigest(),target_sha256=hashlib.sha256(b).hexdigest(),operations=ops,receipt=json.loads((W/'Q9_COMPLETION.json').read_text(encoding='utf-8')));p=W/'q9_archive_patch.json.gz';p.write_bytes(gzip.compress(json.dumps(payload,ensure_ascii=False).encode('utf-8')));print('PATCH_BYTES',p.stat().st_size)
else:
 R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase122_HLCA_supplement_20261026';P=O/'delivery';m=json.loads(gzip.decompress((P/'q9_archive_patch.json.gz').read_bytes()));old=(P/'CAF_reviewer_QA_expanded_reviewed_20261027.pdf').read_bytes();assert hashlib.sha256(old).hexdigest()==m['source_sha256'];new=b''.join(old[op[1]:op[1]+op[2]] if op[0]=='copy' else base64.b64decode(op[1]) for op in m['operations']);assert hashlib.sha256(new).hexdigest()==m['target_sha256'];dest=P/'CAF_reviewer_QA_expanded_Q9_20261027.pdf';dest.write_bytes(new);r=m['receipt'];r.update(server_file=str(dest),server_sha256_verified=True);(P/'Q9_FINAL_RECEIPT.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
 for pointer in [R/'LATEST_CAF_INTEGRATED_CUSTOMER_REPORT.json',O/'STATUS.json',R/'LATEST_HLCA_SUPPLEMENT.json']:
  d=json.loads(pointer.read_text(encoding='utf-8'));d['reviewer_QA_expanded']=r;pointer.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(dict(status='PASS',sha256=m['target_sha256'],changed_question='Q9 only')))
