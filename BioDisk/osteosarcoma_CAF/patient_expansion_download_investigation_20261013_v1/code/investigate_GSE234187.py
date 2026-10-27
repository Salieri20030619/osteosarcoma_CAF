
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,json,gzip,hashlib
from scipy import io,sparse
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1';D=O/'downloads/GSE234187'
assert json.loads((O/'resources/GSE234187_DOWNLOAD_COMPLETE.json').read_text())['status']=='complete_verified'
with gzip.open(D/'GSM7453693_matrix.mtx.gz','rb') as f:x=io.mmread(f).T.tocsr().astype(np.float64)
features=pd.read_csv(D/'GSM7453693_features.tsv.gz',sep='\t',header=None);barcodes=pd.read_csv(D/'GSM7453693_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str);assert x.shape==(len(barcodes),len(features));assert len(barcodes)==13684 and barcodes.is_unique;assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data))
symbols=features[1].astype(str);counts=np.asarray(x.sum(1)).ravel();nfeat=np.diff(x.indptr);mito=np.asarray(x[:,symbols.str.upper().str.startswith('MT-').to_numpy()].sum(1)).ravel()/np.maximum(counts,1)
lineage=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text());programs=json.loads((P/'phase86_functional_bone_lung_20261012/resources/PROGRAMS.json').read_text());maps={g:np.flatnonzero(symbols.eq(g)) for g in set(sum(lineage.values(),[])+sum(programs.values(),[]))};coverage=[]
def score_panels(panels,prefix):
 result={}
 for name,genes in panels.items():
  measured=[g for g in genes if len(maps.get(g,[]))];val=np.zeros(x.shape[0])
  for g in measured:
   ix=maps[g]
   # Multiple Ensembl rows for one symbol are summed on counts before transformation.
   v=np.log1p(np.asarray(x[:,ix].sum(1)).ravel()/np.maximum(counts,1)*1e4);val+=v
  result[prefix+name]=val/len(measured) if measured else np.full(x.shape[0],np.nan)
  coverage.append(dict(type=prefix,program=name,defined=len(genes),measured=len(measured),fraction=len(measured)/len(genes),missing=';'.join(g for g in genes if g not in measured)))
 return result
scores=score_panels(lineage,'lineage__');scores.update(score_panels(programs,'program__'))
L=np.column_stack([scores['lineage__'+k] for k in lineage]);top=np.nanmax(L,axis=1);ix=np.nanargmax(L,axis=1);winner=np.array(list(lineage))[ix].astype(object);ties=(np.abs(L-top[:,None])<1e-12).sum(1)!=1;winner[ties|np.isnan(top)|(top<=0)]='Unresolved'
meta=pd.DataFrame(dict(cell_id='GSE234187__'+barcodes,barcode=barcodes,patient_id='GSE234187_patient1',site='Lung',n_counts=counts,n_features=nfeat,mitochondrial_fraction=mito,provisional_expression_lineage=winner,CNV_status='not_assessed',doublet_status='not_assessed',formal_CAF_admission=False))
for name,v in scores.items():meta[name]=v
candidates=json.loads((P/'CAF_ALL14_RESEARCH_CURRENT.json').read_text())['retained_genes'];genes=sorted(set(sum(lineage.values(),[])+candidates+['FAP','PDPN','CD248','GPC1','THY1','PTPRC','EPCAM','COL11A1','POSTN']))
gene_summary=[]
for g in genes:
 js=np.flatnonzero(symbols.eq(g));v=np.asarray(x[:,js].sum(1)).ravel() if len(js) else np.full(x.shape[0],np.nan)
 if g in candidates:meta['gene_log1pCP10k__'+g]=np.log1p(v/np.maximum(counts,1)*1e4)
 for label in np.unique(winner):
  keep=winner==label
  gene_summary.append(dict(provisional_lineage=label,gene=g,assayed=bool(len(js)),cells=int(keep.sum()),detected_cells=int(np.sum(v[keep]>0)) if len(js) else None,detection_fraction=float(np.mean(v[keep]>0)) if len(js) else None,mean_log1pCP10k=float(np.mean(np.log1p(v[keep]/np.maximum(counts[keep],1)*1e4))) if len(js) else None))
meta.to_csv(O/'tables/GSE234187_cell_investigation.tsv.gz',sep='\t',index=False);pd.DataFrame(gene_summary).to_csv(O/'tables/GSE234187_lineage_gene_evidence.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/GSE234187_program_coverage.tsv',sep='\t',index=False)
summary=meta.groupby('provisional_expression_lineage').agg(cells=('cell_id','size'),median_counts=('n_counts','median'),median_features=('n_features','median'),median_mito=('mitochondrial_fraction','median')).reset_index();summary.to_csv(O/'tables/GSE234187_lineage_QC_summary.tsv',sep='\t',index=False)
receipt={'status':'completed_descriptive_investigation','cells':int(x.shape[0]),'features':int(x.shape[1]),'nonzeros':int(x.nnz),'independent_patients':1,'zero_library_cells':int((counts==0).sum()),'provisional_fibroblast_expression_cells':int((winner=='Fibroblast').sum()),'formal_new_CAF_patients':0,'reason_not_formally_admitted':'Single-case source; competitive expression screen is provisional; doublet and malignant-lineage adjudication not yet established','no_new_QC_thresholds_applied':True,'source_counts_unchanged':True,'no_cell_level_inferential_P':True}
(O/'resources/GSE234187_INVESTIGATION_COMPLETE.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True);print(summary.to_string(index=False),flush=True)
