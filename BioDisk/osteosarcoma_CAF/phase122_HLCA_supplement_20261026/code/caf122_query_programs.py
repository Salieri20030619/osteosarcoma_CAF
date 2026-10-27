from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json,anndata as ad
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');I=R/'phase119_Bone_Lung_functional_subgroups_20261024';O=R/'phase122_HLCA_supplement_20261026'
p=json.loads((I/'resources/PROGRAMS.json').read_text());desired=set(sum(p.values(),[]))|set(['CCL2','COL1A1','CTHRC1','MMP19','SERPINE1','HIF1A','SERPINH1','P4HA1','P4HA2','P4HA3','PLOD1','PLOD2','PLOD3','PCOLCE','PCOLCE2','BMP1','LOX','LOXL1','LOXL2','FNDC1','COMP','MMP2','COL11A1','COL10A1'])
genes=sorted(desired&set(pd.read_csv(I/'tables/genes.tsv',sep='\t').symbol));m=pd.read_csv(O/'tables/frozen_registry.tsv.gz',sep='\t');fit=~m.doublet
can=pd.read_csv(R/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');parts=[];order=[];coverage=[]
for item in json.loads((I/'resources/OBJECTS.json').read_text()):
 a=ad.read_h5ad(item['path'],backed='r');src=item['source'];ix=m.index[fit&m.source.eq(src)&m.cell_id.isin(a.obs_names)].to_numpy()
 if not len(ix):a.file.close();continue
 if src!='NCI':mp={r.symbol:a.var_names.get_indexer([r.feature])[0] for r in can[can.dataset.eq(src)].itertuples()}
 else:
  sy=a.var.symbol.astype(str);ok=~sy.duplicated(keep=False);mp={s:j for j,s in enumerate(sy) if ok.iloc[j]}
 ids=a.obs_names.get_indexer(m.loc[ix,'cell_id']);raw=a.X[ids,:].tocsr();den=np.asarray(raw.sum(1)).ravel();assert np.allclose(den,m.loc[ix,'library_counts'])
 assert all(mp.get(g,-1)>=0 for g in genes),'Shared feature missing'
 parts.append(raw[:,[mp[g] for g in genes]]);order.extend(ix);coverage.append(dict(source=src,cells=len(ix),genes=len(genes)));a.file.close()
 allix=m.index[fit].to_numpy()
X=sparse.vstack(parts).tocsr()[pd.Index(order).get_indexer(allix)];obs=m.loc[allix].copy().set_index('cell_id',drop=False);ad.AnnData(X,obs=obs,var=pd.DataFrame(index=genes)).write_h5ad(O/'objects/query_program_counts.h5ad',compression='gzip');(O/'resources/QUERY_PROGRAM_COUNTS_COMPLETE.json').write_text(json.dumps(dict(cells=X.shape[0],genes=X.shape[1],coverage=coverage),indent=2),encoding='utf-8');print('COMPLETE',X.shape)
