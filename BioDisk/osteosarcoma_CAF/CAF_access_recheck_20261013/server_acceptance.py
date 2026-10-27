from pathlib import Path
import json,hashlib,pandas as pd,gzip,sys
ROOT=Path(__file__).resolve().parent
P=ROOT.parent;old=P/'patient_expansion_download_investigation_20261013_v1'
d=json.loads((ROOT/'NCC_DOWNLOAD_COMPLETE.json').read_text())
for f in ([] if '--barcodes-only' in sys.argv else d['files']):
 p=ROOT/'downloads/NCC'/f['file'];h=hashlib.sha256()
 with p.open('rb') as inp:
  for b in iter(lambda:inp.read(1024*1024),b''):h.update(b)
 assert p.stat().st_size==f['bytes'] and h.hexdigest()==f['sha256']
query=pd.read_csv(ROOT/'tables/NCC_BARCODES.tsv.gz',sep='\t');idx=pd.read_csv(old/'tables/ORIGINAL_SOURCE_BARCODE_INDEX.tsv',sep='\t')
refs={r.sample:set(gzip.open(r.barcode_list,'rt').read().splitlines()) for r in idx.itertuples()}
rows=[]
for name,group in query.groupby('sample'):
 bc=set(group.barcode16)
 for r in idx.itertuples():
  overlap=len(bc&refs[r.sample]);rows.append(dict(NCC_sample=name,reference_sample=r.sample,reference_dataset=r.dataset,shared_barcodes=overlap,query_unique_barcodes=len(bc),reference_unique_barcodes=len(refs[r.sample]),query_overlap_fraction=overlap/len(bc)))
df=pd.DataFrame(rows);df.to_csv(ROOT/'tables/NCC_VS_ORIGINAL_BARCODES.tsv',sep='\t',index=False)
best=df.sort_values('query_overlap_fraction',ascending=False).groupby('NCC_sample').head(1).sort_values('NCC_sample')
best.to_csv(ROOT/'tables/NCC_BEST_ORIGINAL_BARCODE_OVERLAP.tsv',sep='\t',index=False)
receipt={'status':'PASS','files_verified':len(d['files']),'bytes':d['bytes'],'reference_samples':len(idx),'NCC_samples_barcode_comparable':query['sample'].nunique(),'BS6_excluded_missing_labels':True,'best_overlap_min':float(best.query_overlap_fraction.min()),'best_overlap_max':float(best.query_overlap_fraction.max()),'barcode_overlap_does_not_prove_patient_independence':True,'scientific_baseline_unchanged':True}
if '--barcodes-only' in sys.argv:
 receipt['files_verified']=0;receipt['status']='barcode_comparison_complete_file_transfer_pending'
out='NCC_BARCODE_COMPARISON.json' if '--barcodes-only' in sys.argv else 'SERVER_ACCEPTANCE.json'
(ROOT/out).write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt));print(best.to_string(index=False))
