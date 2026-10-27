import pathlib,json,concurrent.futures
import numpy as np,pandas as pd,anndata as ad
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');q=root/'phase95_NCI_expanded_reanalysis_20261016';p=root/'phase100_NCI_completion_20261020'
markers=['DCN','LUM','COL1A1','COL1A2','COL3A1','COL6A3','FN1','PDGFRA','FAP','PDPN','PTPRC','LST1','TYROBP','LYZ','CSF1R','PECAM1','VWF','KDR','EMCN','RGS5','MYH11','ACTA2','PDGFRB','RUNX2','SATB2','ALPL','IBSP','SP7','ACP5','CTSK','MMP9','EPCAM','SFTPB','KRT8','KRT18','CD3D','CD3E','CD79A']
run=json.loads((q/'COMPUTATION_FINISHED.json').read_text(encoding='utf-8'))['samples']
def one(r):
 sm=r['sample']
 if not r['reference_pass']:return []
 c=pd.read_csv(q/'tables'/(sm+'_all_cell_evidence.tsv.gz'),sep='\t',encoding='utf-8')
 ids=c[c.mechanical_expanded_candidate].copy();a=ad.read_h5ad(q/'inputs'/sm/'available_counts.h5ad',backed='r')
 ci=a.obs_names.get_indexer(ids.cell_id);assert (ci>=0).all()
 sy=a.var.symbol.astype(str).to_numpy();idx=np.flatnonzero(np.isin(sy,markers))
 x=a[ci,:].X[:,idx].toarray();den=ids.n_counts.to_numpy();assert (den>0).all();y=np.log1p(x/den[:,None]*1e4)
 rows=[]
 for cl in sorted(ids.cluster.unique()):
  mask=ids.cluster.eq(cl).to_numpy()
  for j,col in enumerate(idx):rows.append(dict(sample=sm,cluster=str(cl),gene=sy[col],candidate_cells=int(mask.sum()),detected_fraction=float((x[mask,j]>0).mean()),mean_log1pCP10k=float(y[mask,j].mean())))
 a.file.close();return rows
with concurrent.futures.ThreadPoolExecutor(2) as e:parts=list(e.map(one,run))
pd.DataFrame([r for z in parts for r in z]).to_csv(p/'tables/candidate_identity_marker_audit.tsv',sep='\t',encoding='utf-8',index=False)
print('MARKER_AUDIT_COMPLETE',sum(map(len,parts)),flush=True)
