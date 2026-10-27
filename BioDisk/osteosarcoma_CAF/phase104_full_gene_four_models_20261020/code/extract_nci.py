import pathlib,json,hashlib
import numpy as np,pandas as pd,anndata as ad
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase104_full_gene_four_models_20261020'
for d in ['code','resources','tables','objects','figures','report','logs']: (P/d).mkdir(parents=True,exist_ok=True)
q=R/'phase95_NCI_expanded_reanalysis_20261016'
run=json.loads((q/'COMPUTATION_FINISHED.json').read_text())['samples'];features=None;counts={};totals={};cells={}
for rr in run:
 sm=rr['sample'];m=pd.read_csv(q/'tables'/(sm+'_all_cell_evidence.tsv.gz'),sep='\t',usecols=['cell_id','author_CAF','Ann_Level0','Ann_Level2','n_counts'])
 m=m[m.author_CAF];assert m.Ann_Level0.eq('Host').all() and m.Ann_Level2.eq('CAFs').all()
 a=ad.read_h5ad(q/'inputs'/sm/'available_counts.h5ad',backed='r')
 ix=a.obs_names.get_indexer(m.cell_id);assert (ix>=0).all();z=a[ix,:].to_memory();a.file.close()
 if features is None:features=z.var_names.to_numpy();symbols=z.var.symbol.to_numpy()
 assert np.array_equal(features,z.var_names)
 raw=np.asarray(z.X.sum(0)).ravel();den=float(z.X.sum());assert den==m.n_counts.sum()
 pid='NCI_PATIENT_'+sm.rsplit('_',1)[-1].removesuffix('b')
 counts[pid]=counts.get(pid,np.zeros(len(features)))+raw;totals[pid]=totals.get(pid,0)+den;cells[pid]=cells.get(pid,0)+len(m)
 print(sm,len(m),flush=True)
ids=sorted(counts);mat=np.vstack([counts[i] for i in ids]);lib=np.array([totals[i] for i in ids])
np.savez_compressed(P/'objects/NCI_allgene_author_CAF.npz',features=features.astype(str),symbols=symbols.astype(str),donors=np.array(ids),raw=mat,library_UMIs=lib,cells=np.array([cells[i] for i in ids]),expression=np.log2(mat/lib[:,None]*1e4+.1))
assert len(ids)==17 and sum(cells.values())==1530
print('NCI_COMPLETE',mat.shape,flush=True)
