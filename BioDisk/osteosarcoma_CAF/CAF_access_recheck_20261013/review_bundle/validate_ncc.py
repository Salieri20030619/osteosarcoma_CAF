from pathlib import Path
import h5py,numpy as np,pandas as pd,json,re,gzip,zipfile,io
from scipy import sparse
ROOT=Path(__file__).resolve().parent;D=ROOT/'downloads/NCC';T=ROOT/'tables';T.mkdir(exist_ok=True)
metadata=pd.read_csv(D/'Patient metadata',sep='\t');rows=[];bcrows=[]
for p in sorted(D.glob('*')):
 if p.suffix not in ('.h5','.h5ad'):continue
 sample=p.stem.split('-')[0]
 with h5py.File(p) as h:
  if 'matrix' in h:
   g=h['matrix'];shape=tuple(g['shape'][:]);genes=g['features/name'].asstr()[:];ids=g['features/id'].asstr()[:];bc=g['barcodes'].asstr()[:]
   x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=shape)
   encoding='10x_genes_by_cells_csc'
  else:
   g=h['X'];shape=tuple(g.attrs['shape']);genes=h['obs/_index'].asstr()[:];bc=h['var/_index'].asstr()[:];ids=genes
   x=sparse.csr_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=shape)
   encoding='h5ad_36601_by_7291_numeric_indices_only_axes_unverified'
  assert x.shape==(len(genes),len(bc))
  assert len(set(bc))==len(bc)
  assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
  x.check_format(full_check=True)
  # Check axis identity independently of shape using barcode patterns and a gene marker panel.
  barcode_pattern=sum(bool(re.search('[ACGT]{16}',b)) for b in bc)
  gene_markers=[g for g in ['COL1A1','DCN','LUM','PTPRC','EPCAM','RUNX2'] if g in set(genes)]
  labeled_axes=len(gene_markers)>=4 and barcode_pattern==len(bc)
  if p.suffix=='.h5':assert labeled_axes
  totals=np.asarray(x.sum(axis=0)).ravel()
  mr=metadata.loc[metadata.Sample==sample].iloc[0].to_dict()
  row=dict(file=p.name,sample=sample,matrix_rows=x.shape[0],matrix_columns=x.shape[1],genes=len(genes) if labeled_axes else None,barcodes=len(bc) if labeled_axes else None,nonzero_entries=int(x.nnz),zero_sum_columns=int((totals==0).sum()),median_column_sum=float(np.median(totals)),encoding=encoding,axes_labeled=labeled_axes,finite_nonnegative_integer=True,**{k:str(v) for k,v in mr.items()})
  rows.append(row)
  if labeled_axes:bcrows.extend(dict(sample=sample,barcode=b,barcode16=re.search('[ACGT]{16}',b).group()) for b in bc)
  print(json.dumps(row),flush=True)
pd.DataFrame(rows).to_csv(T/'NCC_MATRIX_AND_PATIENT_INVENTORY.tsv',sep='\t',index=False)
pd.DataFrame(bcrows).to_csv(T/'NCC_BARCODES.tsv.gz',sep='\t',index=False)
spatial={}
with zipfile.ZipFile(D/'SP_BS3.zip') as z:
 for name in z.namelist():
  if name.endswith('.gz'):
   raw=gzip.decompress(z.read(name))
   if name.endswith('matrix.mtx.gz'):
    from scipy.io import mmread
    mx=mmread(io.BytesIO(raw));assert np.isfinite(mx.data).all() and (mx.data>=0).all() and np.equal(mx.data,np.floor(mx.data)).all()
    spatial.update(genes=int(mx.shape[0]),spots=int(mx.shape[1]),nonzeros=int(mx.nnz))
   elif name.endswith('barcodes.tsv.gz'):spatial['barcode_rows']=len(raw.splitlines())
   elif name.endswith('features.tsv.gz'):spatial['feature_rows']=len(raw.splitlines())
 assert spatial['spots']==spatial['barcode_rows'] and spatial['genes']==spatial['feature_rows']
receipt={'status':'file_integrity_pass_BS6_semantic_labels_missing','single_cell_files':len(rows),'labeled_single_cell_files':sum(r['axes_labeled'] for r in rows),'sample_records':len(metadata),'unique_patient_IDs':metadata.ID.nunique(),'confirmed_barcodes_total':sum(r['barcodes'] or 0 for r in rows),'BS6_unlabeled_dimensions':[36601,7291],'BS6_not_admitted_to_gene_or_barcode_analysis':True,'author_postQC_cell_count':49356,'no_QC_or_CAF_admission_performed':True,'spatial':spatial,'patient_ID_not_global_deduplication':True,'paired_prepost_patient':'Patient4: BS4 / BS14','lung_specimens_documented':0}
(ROOT/'NCC_MATRIX_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2),encoding='utf8');print(json.dumps(receipt))
