
from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1'
m=pd.read_csv(O/'tables/OsteoCAR/patient_prim_author_metadata.tsv.gz',sep='\t',low_memory=False);scores=pd.read_csv(O/'tables/OsteoCAR/patient_prim_cell_scores.tsv.gz',sep='\t');m=m.merge(scores,on='cell_id',validate='one_to_one');mp=pd.read_csv(O/'tables/OsteoCAR/patient_prim_best_overlap_candidates.tsv',sep='\t');genes=[c.split('__',1)[1] for c in scores.columns if c.startswith('gene_log1pCP10k__')];rows=[]
for r in mp.itertuples():
 g=m[m.sample_name.eq(r.author_sample)].copy();g['barcode16']=g.cell_id.str.extract(r'([ACGT]{16})')[0];ev=pd.read_csv(P/'phase27_identity_review_20260926/tables'/f'{r.project_sample}_review.tsv.gz',sep='\t',usecols=['cell_id','n_counts']);ev['barcode16']=ev.cell_id.str.extract(r'([ACGT]{16})')[0];g=g.merge(ev,on='barcode16',suffixes=('_author','_project'),validate='one_to_one')
 library_error=g.nCount_RNA-g.n_counts
 selected=g.sort_values('cell_id_author').iloc[np.linspace(0,len(g)-1,min(40,len(g)),dtype=int)].copy()
 f=P/'phase26_client_reanalysis_20260922/inputs'/r.project_sample/'full_called_counts.h5ad';a=ad.read_h5ad(f,backed='r');mapping={}
 for j,sy in enumerate(a.var.symbol.astype(str)):
  for gene in sy.split(';'):
   if gene in genes:mapping.setdefault(gene,[]).append(j)
 js=sorted(set(j for v in mapping.values() for j in v));ix=a.obs_names.get_indexer(selected.cell_id_project);assert (ix>=0).all();order=np.argsort(ix);rev=np.argsort(order);x=a[ix[order],js].X.tocsr()[rev];pos={v:j for j,v in enumerate(js)};a.file.close();matches=0;comparisons=0;nonzero=0;errors=[]
 for gene in genes:
  if gene not in mapping:continue
  transformed=selected['gene_log1pCP10k__'+gene].to_numpy(float);observed=selected.full_library_counts.to_numpy()*np.expm1(transformed)/1e4;valid=np.isfinite(observed);observed=np.rint(observed[valid]);actual=np.asarray(x[:,[pos[j] for j in mapping[gene]]].sum(1)).ravel()[valid];err=np.abs(observed-actual);matches+=int((err==0).sum());comparisons+=len(err);nonzero+=int(((observed>0)|(actual>0)).sum());errors.extend(err)
 rows.append(dict(author_sample=r.author_sample,project_sample=r.project_sample,matched_barcodes=len(g),fingerprint_cells=len(selected),gene_cell_comparisons=comparisons,nonzero_comparisons=nonzero,exact_count_matches=matches,count_match_fraction=matches/comparisons if comparisons else None,max_gene_count_error=max(errors) if errors else None,nCount_RNA_matches_full_source=int(np.isclose(library_error,0).sum()),library_comparisons=len(g),max_library_count_error=float(np.abs(library_error).max())))
 print('FINGERPRINT',r.author_sample,matches,comparisons,flush=True)
out=pd.DataFrame(rows);out.to_csv(O/'tables/OsteoCAR/primary_raw_count_fingerprint_verification.tsv',sep='\t',index=False);print(out.to_string(index=False));(O/'resources/PRIMARY_COUNT_FINGERPRINT_COMPLETE.json').write_text(json.dumps({'status':'completed','samples':len(out),'gene_cell_comparisons':int(out.gene_cell_comparisons.sum()),'exact_matches':int(out.exact_count_matches.sum()),'all_source_library_counts_identical':bool((out.nCount_RNA_matches_full_source==out.library_comparisons).all())},indent=2))
