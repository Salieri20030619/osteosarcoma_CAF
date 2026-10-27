from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import sys,json,re
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase95_NCI_expanded_reanalysis_20261016';sm=sys.argv[1];d=O/'inputs'/sm
assert (d/'EXPORT_COMPLETE').exists()
c=pd.read_csv(d/'author_cells.tsv',sep='\t',encoding='utf-8');g=pd.read_csv(d/'genes.tsv',sep='\t',encoding='utf-8')
x=mmread(d/'counts.mtx.gz').T.tocsr();assert x.shape==(len(c),len(g)) and np.array_equal(np.asarray(x.sum(1)).ravel(),c.n_counts)
assert c.cell_id.is_unique and g.feature_id.is_unique
panels=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text(encoding='utf-8'))
den=np.asarray(x.sum(1)).ravel();y=x.multiply((1e4/np.maximum(den,1))[:,None]).tocsr();y.data=np.log1p(y.data)
coverage=[]
for name,genes in panels.items():
 idx=np.flatnonzero(g.symbol.isin(genes));c['phase27_'+name]=np.asarray(y[:,idx].mean(1)).ravel() if len(idx) else np.nan
 c['det_'+name]=np.asarray((x[:,idx]>0).sum(1)).ravel() if len(idx) else 0
 coverage.append({'panel':name,'measured_object_features':len(idx),'defined':len(genes),'missing_features':';'.join(sorted(set(genes)-set(g.symbol))),'source_specific_measurement_mask':'unavailable_in_merged_author_object'})
values=c[['phase27_'+n for n in panels]].to_numpy();mx=np.nanmax(values,axis=1);label=np.array(list(panels),dtype=object)[np.nanargmax(values,axis=1)];label[(mx<=0)|(np.sum(values==mx[:,None],axis=1)!=1)]='Unresolved'
c['provisional_lineage']=label;c['author_CAF']=c.Ann_Level0.eq('Host')&c.Ann_Level2.eq('CAFs')
c['patient_id']='unverified';c['putative_patient_group']='NCI_PUTATIVE__'+re.sub(r'b$','',str(c.sample_name.iloc[0]).split('_')[-1]);c['patient_group_status']='suffix_hypothesis_not_verified_independent_N'
a=ad.AnnData(y,obs=pd.DataFrame(index=c.cell_id),var=pd.DataFrame(index=g.feature_id))
sc.settings.n_jobs=4;sc.pp.highly_variable_genes(a,n_top_genes=min(2000,a.n_vars),flavor='seurat')
z=a[:,a.var.highly_variable].copy();sc.pp.scale(z,max_value=10);sc.tl.pca(z,n_comps=30,svd_solver='arpack',random_state=20261016);sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=20261016);sc.tl.leiden(z,resolution=.8,key_added='cluster',random_state=20261016);sc.tl.umap(z,random_state=20261016)
c['cluster']=z.obs.cluster.astype(str).to_numpy();c['UMAP1']=z.obsm['X_umap'][:,0];c['UMAP2']=z.obsm['X_umap'][:,1]
c['reference']=False;c['heldout']=False;c['reference_used']=False;c['heldout_clean']=False;c['reference_role']='query'
c.to_csv(d/'cells.tsv',sep='\t',index=False,encoding='utf-8');pd.DataFrame(coverage).to_csv(d/'panel_coverage.tsv',sep='\t',index=False,encoding='utf-8')
markers=[]
global_sum=np.asarray(y.sum(0)).ravel()
for cluster,ix in c.groupby('cluster').groups.items():
 ix=np.array(list(ix));inside=np.asarray(y[ix].mean(0)).ravel();outside=(global_sum-inside*len(ix))/max(len(c)-len(ix),1);delta=inside-outside;top=np.argsort(-delta)[:40];det=np.asarray((x[ix][:,top]>0).mean(0)).ravel()
 for rank,j in enumerate(top):markers.append({'sample':sm,'cluster':cluster,'rank':rank+1,'symbol':g.symbol.iloc[j],'mean_log1pCP10k':inside[j],'difference_from_other_cells':delta[j],'detected_fraction':det[rank],'method':'descriptive_ranking_no_inferential_P'})
pd.DataFrame(markers).to_csv(O/'tables'/(sm+'_cluster_markers.tsv'),sep='\t',index=False,encoding='utf-8')
c.groupby(['cluster','provisional_lineage','Ann_Level0','Ann_Level2']).size().rename('cells').reset_index().to_csv(O/'tables'/(sm+'_cluster_lineage.tsv'),sep='\t',index=False,encoding='utf-8')
obs=c.set_index('cell_id').copy()
for col in obs:
 if obs[col].dtype==object:obs[col]=obs[col].fillna('').astype(str)
raw=ad.AnnData(x,obs=obs,var=g.set_index('feature_id'));raw.uns['input_scope']='All available NCI cells in author-processed counts; no additional global QC filter';raw.write_h5ad(d/'available_counts.h5ad',compression='gzip')
(d/'PREPARATION_COMPLETE').write_text('complete\n',encoding='utf-8')
print(json.dumps({'sample':sm,'cells':len(c),'clusters':c.cluster.nunique(),'fibroblast_winner':int(c.provisional_lineage.eq('Fibroblast').sum()),'author_CAF':int(c.author_CAF.sum())}),flush=True)
