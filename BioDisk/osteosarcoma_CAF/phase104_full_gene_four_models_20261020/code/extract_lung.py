import pathlib,gzip,json
import numpy as np,pandas as pd
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');P=R/'phase104_full_gene_four_models_20261020';Q=R/'phase97_normal_background_20261020'
m=pd.read_csv(Q/'tables/lung_selected_cell_metadata.tsv',sep='\t',dtype=str).set_index('cellId')
groups={k:z.index for k,z in m.groupby('patient')};assert len(groups)==3 and min(map(len,groups.values()))>=20
ids=sorted(groups);counts=[];features=[]
with gzip.open(Q/'inputs/lung_counts_verified.tsv.gz','rt',encoding='utf-8') as f:
 cells=pd.Index(f.readline().rstrip('\r\n').split('\t')[1:]);pos=[cells.get_indexer(groups[i]) for i in ids];assert all((ix>=0).all() for ix in pos)
 for line in f:
  label,values=line.rstrip('\r\n').split('\t',1);a=np.fromstring(values,sep='\t')
  assert len(a)==len(cells) and np.isfinite(a).all() and (a>=0).all() and (a==np.floor(a)).all()
  features.append(label);counts.append([a[ix].sum() for ix in pos])
  if len(features)%5000==0:print('LUNG',len(features),flush=True)
mat=np.asarray(counts).T;lib=mat.sum(1)
for i,donor in enumerate(ids):assert lib[i]==m.loc[groups[donor],'nUMI'].astype(float).sum()
assert mat.shape==(3,26485)
np.savez_compressed(P/'objects/normal_lung_allgene.npz',features=np.array(features),donors=np.array(ids),raw=mat,library_UMIs=lib,cells=np.array([len(groups[i]) for i in ids]),expression=np.log2(mat/lib[:,None]*1e4+.1))
print('LUNG_COMPLETE',mat.shape,flush=True)
