
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,gzip,json,shutil
from scipy import io
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1';D=O/'downloads/GSE234187';meta_file=O/'tables/GSE234187_cell_investigation.tsv.gz';meta=pd.read_csv(meta_file,sep='\t');features=pd.read_csv(D/'GSM7453693_features.tsv.gz',sep='\t',header=None);stable=features[0].astype(str).str.replace(r'\.\d+$','',regex=True)
with gzip.open(D/'GSM7453693_matrix.mtx.gz','rb') as f:x=io.mmread(f).T.tocsr()
barcodes=pd.read_csv(D/'GSM7453693_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str);assert np.array_equal(meta.barcode.astype(str),barcodes)
lib=np.asarray(x.sum(1)).ravel();assert np.array_equal(lib,meta.n_counts)
reg=pd.read_csv(P/'phase82_allCAF_crossdataset_markers_20261012/tables/primary_same_layer_crossdataset_markers.tsv',sep='\t');rows=[];summ=[];original=pd.read_csv(O/'tables/GSE234187_lineage_gene_evidence.tsv',sep='\t');gene_names=set(reg.symbol)
shutil.copy2(meta_file,O/'resources/GSE234187_cell_investigation_before_gene_id_resolution.tsv.gz')
for r in reg.itertuples():
 ids=str(r.feature_ids).split(';');js=np.flatnonzero(stable.isin(ids));raw=np.asarray(x[:,js].sum(1)).ravel() if len(js) else np.full(len(meta),np.nan);v=np.log1p(raw/np.maximum(lib,1)*1e4);meta['gene_log1pCP10k__'+r.symbol]=v
 rows.append(dict(canonical_gene=r.symbol,stable_feature_ids=r.feature_ids,matched_rows=len(js),source_symbols=';'.join(features.iloc[js,1].astype(str)),assayed=bool(len(js))))
 for label,g in meta.groupby('provisional_expression_lineage'):
  ix=g.index;summ.append(dict(provisional_lineage=label,gene=r.symbol,assayed=bool(len(js)),cells=len(g),detected_cells=int(np.sum(raw[ix]>0)) if len(js) else None,detection_fraction=float(np.mean(raw[ix]>0)) if len(js) else None,mean_log1pCP10k=float(np.mean(v[ix])) if len(js) else None))
meta.to_csv(meta_file,sep='\t',index=False);pd.concat([original[~original.gene.isin(gene_names)],pd.DataFrame(summ)],ignore_index=True).to_csv(O/'tables/GSE234187_lineage_gene_evidence.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(O/'tables/GSE234187_candidate_stable_gene_mapping.tsv',sep='\t',index=False)
receipt={'status':'completed','canonical_candidates_retained':14,'assayed_candidates':sum(r['assayed'] for r in rows),'mapping':'frozen stable Ensembl gene IDs; missing IDs are not filled with zero','cell_labels_unchanged':True,'source_counts_unchanged':True};(O/'resources/GSE234187_GENE_ID_RESOLUTION_COMPLETE.json').write_text(json.dumps(receipt,indent=2));print(pd.DataFrame(rows).to_string(index=False));print(json.dumps(receipt))
