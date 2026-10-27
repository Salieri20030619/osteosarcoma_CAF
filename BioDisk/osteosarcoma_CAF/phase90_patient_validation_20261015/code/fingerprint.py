from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,re,hashlib,datetime
import pandas as pd,numpy as np,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase90_patient_validation_20261015'
pat=re.compile(r'([ACGT]{16})')
def bc(s):
 m=pat.search(str(s));return m.group(1) if m else None
manifest=pd.read_csv(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',sep='\t')
ncc=pd.read_csv(O/'resources/manifest.tsv',sep='\t');cache={}
for sm in ncc['sample']:
 a=ad.read_h5ad(P/'phase87_NCC_preprocessing_20261014/inputs'/sm/'counts.h5ad')
 ids=[bc(s) for s in a.obs_names];assert None not in ids and len(set(ids))==len(ids)
 cache[sm]=(a,dict(zip(ids,range(len(ids)))))
rows=[];matches=[];sources=[]
for r in manifest.itertuples():
 f=P/'phase26_client_reanalysis_20260922/inputs'/r.sample/'full_called_counts.h5ad'
 stat=f.stat();a=ad.read_h5ad(f,backed='r');bcs=[bc(s) for s in a.obs_names]
 counts=pd.Series(bcs).value_counts();ambig=set(counts[counts>1].index)
 ix={b:i for i,b in enumerate(bcs) if b is not None and b not in ambig}
 sources.append({'sample':r.sample,'path':str(f),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'duplicate_16nt_barcodes_excluded':len(ambig)})
 for sm,(b,bi) in cache.items():
  common=sorted(set(a.var_names)&set(b.var_names));shared=sorted(set(ix)&set(bi))
  row={'NCC_sample':sm,'baseline_sample':r.sample,'baseline_dataset':r.dataset,'shared_genes':len(common),'shared_unambiguous_barcodes':len(shared),'exact_count_vectors':0,'nonzero_exact_count_vectors':0,'max_nonzero_gene_overlap':0}
  if common and shared:
   # Read a small union of matched barcode rows; never densify full matrices.
   ai=np.array([ix[s] for s in shared]);order=np.argsort(ai);inv=np.argsort(order)
   x=a[ai[order],:].X.tocsr()[inv][:,a.var_names.get_indexer(common)]
   y=b.X[[bi[s] for s in shared],:][:,b.var_names.get_indexer(common)].tocsr()
   diff=(x-y).tocsr();diff.eliminate_zeros()
   equal=np.asarray(diff.getnnz(1)==0);nz=np.asarray((x>0).sum(1)).ravel()
   row['exact_count_vectors']=int(equal.sum());row['nonzero_exact_count_vectors']=int((equal&(nz>0)).sum());row['max_nonzero_gene_overlap']=int(nz.max())
   for j in np.flatnonzero(equal): matches.append({'NCC_sample':sm,'baseline_sample':r.sample,'barcode':shared[j],'shared_genes':len(common),'nonzero_shared_genes':int(nz[j]),'baseline_cell_id':str(a.obs_names[ix[shared[j]]])})
  rows.append(row)
 a.file.close()
 pd.DataFrame(rows).to_csv(O/'tables/NCC_baseline_count_fingerprints.tsv',sep='\t',index=False,encoding='utf-8')
 print('FINGERPRINT',r.sample,flush=True)
pd.DataFrame(matches,columns=['NCC_sample','baseline_sample','barcode','shared_genes','nonzero_shared_genes','baseline_cell_id']).to_csv(O/'tables/NCC_exact_count_matches.tsv',sep='\t',index=False,encoding='utf-8')
out={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'comparisons':len(rows),'NCC_samples':len(cache),'baseline_samples':len(manifest),'matched_vectors':sum(r['exact_count_vectors'] for r in rows),'nonzero_matched_vectors':sum(r['nonzero_exact_count_vectors'] for r in rows),'scope':'Detect exact shared-library cell vectors on common stable gene IDs and unambiguous 16nt barcodes; does not establish absence of same-patient different-library overlap.','sources':sources}
(O/'resources/FINGERPRINT_COMPLETE.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
