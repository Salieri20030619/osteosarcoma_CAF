from common import *
import anndata as ad
samples=['PT1','MR1','MR10A'];raws=[];edges=[];z=pd.read_csv(O/'tables/all_spot_RNA_and_components.tsv.gz',sep='\t')
for sm in samples:
 a=ad.read_h5ad(P/'phase35_external_validation_prep_20260927/objects/GSE293065'/f'{sm}_counts_spatial.h5ad');gm=feature_map(a);raw=pd.DataFrame({'spot_id':a.obs_names,'sample':sm,'library_UMI':np.asarray(a.X.sum(1)).ravel()})
 for g in ['COL1A1','LOX','MMP14','CD8A','CD8B']:raw[g]=np.asarray(a.X[:,gm[g]].toarray()).ravel()
 raws.append(raw);npz=np.load(O/'objects'/f'{sm}_adjacency.npz');edges.append(pd.DataFrame({'sample':sm,'source':npz['row']+1,'target':npz['col']+1}))
pd.concat(raws).to_csv(O/'tables/independent_QA_raw_gene_counts.tsv',sep='\t',index=False);pd.concat(edges).to_csv(O/'tables/independent_QA_adjacency.tsv',sep='\t',index=False)
print('QA_INPUTS_PREPARED',samples)
