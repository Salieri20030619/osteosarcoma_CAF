from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,gzip,hashlib,subprocess
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase91_lung_patient_validation_20261015'
S=P/'patient_expansion_download_investigation_20261013_v1';D=S/'downloads/GSE234187';sm='GSE234187__GSM7453693';out=O/'inputs'/sm
out.mkdir(parents=True,exist_ok=True)
def dump(f,v):f.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8388608),b''):h.update(b)
 return h.hexdigest()
receipt=json.loads((S/'resources/GSE234187_DOWNLOAD_COMPLETE.json').read_text(encoding='utf-8'))
assert all(sha(Path(r['file']))==r['sha256'] for r in receipt['files'])
feat=pd.read_csv(D/'GSM7453693_features.tsv.gz',sep='\t',header=None)
bar=pd.read_csv(D/'GSM7453693_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
x=mmread(D/'GSM7453693_matrix.mtx.gz').T.tocsr();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
assert x.shape==(13684,33694) and feat[0].is_unique and bar.is_unique
assert feat[0].str.startswith('ENSG').all() and (x.data==np.floor(x.data)).all() and (x.data>=0).all()
den=np.asarray(x.sum(1)).ravel();assert (den>0).all()
panels=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text(encoding='utf-8'))
target=sorted(set(sum(panels.values(),[])));symbols=feat[1].astype(str).to_numpy();gm={g:np.flatnonzero(symbols==g) for g in target}
raw=np.column_stack([np.asarray(x[:,gm[g]].sum(1)).ravel() for g in target]).astype(np.float32)
norm=np.log1p(raw/den[:,None]*1e4);scores={};det={};pct={}
for key,gs in panels.items():
 ii=[target.index(g) for g in gs if len(gm[g])]
 scores[key]=norm[:,ii].mean(1);det[key]=(raw[:,ii]>0).sum(1);pct[key]=raw[:,ii].sum(1)/den
v=np.column_stack(list(scores.values()));mx=v.max(1);winner=np.array(list(panels),dtype=object)[v.argmax(1)]
winner[(mx<=0)|(np.sum(v==mx[:,None],axis=1)!=1)]='Unresolved'
cells=pd.DataFrame({'cell_id':'GSE234187__'+bar,'barcode':bar,'tool_id':[sm+'_c'+str(i).zfill(6) for i in range(len(bar))],'sample':sm,'patient_id':'GSE234187_patient1','dataset':'GSE234187','site':'Lung','n_counts':den,'n_genes':np.asarray((x>0).sum(1)).ravel(),'provisional_lineage':winner,'reference':False,'heldout':False})
for k in panels:cells['phase27_'+k]=scores[k];cells['det_'+k]=det[k]
genes=pd.DataFrame({'ensembl':feat[0],'symbol':feat[1]})
genes.to_csv(out/'genes.tsv',sep='\t',index=False,encoding='utf-8');cells.to_csv(out/'cells.tsv',sep='\t',index=False,encoding='utf-8')
if not (out/'counts.mtx.gz').exists():(out/'counts.mtx.gz').symlink_to(D/'GSM7453693_matrix.mtx.gz')
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
with (O/'logs/doublets.log').open('w',encoding='utf-8') as f:r=subprocess.run([R,'--vanilla',str(O/'code/doublets_full.R'),str(out)],stdout=f,stderr=subprocess.STDOUT)
(O/'logs/doublets.exit_code').write_text(str(r.returncode)+'\n',encoding='utf-8');assert r.returncode==0 and (out/'DOUBLETS_COMPLETE').exists()
db=pd.read_csv(out/'full_scDblFinder_calls.tsv',sep='\t').set_index('cell_id')
assert set(db.index)==set(cells.cell_id)
cells['full_scDblFinder_class']=cells.cell_id.map(db.scDblFinder_class)
cells['full_scDblFinder_score']=cells.cell_id.map(db.scDblFinder_score)
# Fixed clean-reference strategy inherited from Phase88; diagnostic minimums are not CAF cell-entry cutoffs.
mes=set(sum([panels[k] for k in ['Fibroblast','Osteogenic','Mural','Myogenic','Chondrocyte']],[]))
mesidx=[target.index(g) for g in mes if g in target];mesfrac=raw[:,mesidx].sum(1)/den
cells['reference_role']='query';rng=np.random.default_rng(20261015)
for group in ['T_NK','B_Plasma','Myeloid','Endothelial','Mast']:
 ix=np.flatnonzero((winner==group)&(det[group]>=2)&(mesfrac<=.002)&cells.full_scDblFinder_class.eq('singlet').to_numpy())
 if len(ix)<50:continue
 rng.shuffle(ix);nh=max(10,int(round(len(ix)*.2)));hold=ix[:nh];train=ix[nh:]
 cells.loc[hold,'heldout']=True;cells.loc[train,'reference']=True
 cells.loc[hold,'reference_role']='holdout_'+group;cells.loc[train,'reference_role']='ref_'+group
assert cells.reference.sum()>=30 and cells.heldout.sum()>=20
cells.to_csv(out/'cells.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame({'cell':cells.tool_id,'group':np.where(cells.reference,'reference','query')}).to_csv(out/'annotations.tsv',sep='\t',index=False,header=False,encoding='utf-8')
print('INPUT_READY',len(cells),int(cells.reference.sum()),int(cells.heldout.sum()),flush=True)
var=genes.rename(columns={'ensembl':'gene_id'}).set_index('gene_id');obs=cells.set_index('cell_id')
a=ad.AnnData(x.copy(),obs=obs,var=var);a.write_h5ad(O/'objects/full_called_counts.h5ad',compression='gzip')
# All cells retained for neighborhood review; software doublets are flagged, never hidden.
b=a.copy();sc.pp.normalize_total(b,target_sum=1e4);sc.pp.log1p(b)
sc.pp.highly_variable_genes(b,n_top_genes=2000,flavor='seurat');b=b[:,b.var.highly_variable].copy()
sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=20261015)
sc.pp.neighbors(b,n_neighbors=15,n_pcs=30,random_state=20261015)
sc.tl.leiden(b,resolution=.8,key_added='cluster',random_state=20261015)
sc.tl.umap(b,random_state=20261015)
cells['cluster']=b.obs['cluster'].astype(str).to_numpy();cells['UMAP1']=b.obsm['X_umap'][:,0];cells['UMAP2']=b.obsm['X_umap'][:,1]
cells.to_csv(out/'cells.tsv',sep='\t',index=False,encoding='utf-8')
b.write_h5ad(O/'objects/review_embedding.h5ad',compression='gzip')
cp=x.multiply(1e4/den[:,None]).tocsr();total=np.asarray(cp.sum(0)).ravel();totaldet=np.asarray((x>0).sum(0)).ravel();markers=[]
for cl,g in cells.groupby('cluster'):
 ii=g.index.to_numpy();n=len(ii);m=len(cells)-n
 pos=np.asarray(cp[ii].sum(0)).ravel();detg=np.asarray((x[ii]>0).sum(0)).ravel()
 logfc=np.log2((pos/n+1)/((total-pos)/m+1));delta=detg/n-(totaldet-detg)/m;rank=logfc*np.maximum(delta,0)
 for j in np.argsort(-rank)[:50]:markers.append({'cluster':cl,'cells':n,'gene_id':feat.iloc[j,0],'symbol':feat.iloc[j,1],'log2FC_CP10k_plus1':float(logfc[j]),'detection_fraction':float(detg[j]/n),'detection_difference':float(delta[j]),'descriptive_rank':float(rank[j])})
pd.DataFrame(markers).to_csv(O/'tables/cluster_markers.tsv',sep='\t',index=False,encoding='utf-8')
pd.crosstab(cells.cluster,cells.provisional_lineage).to_csv(O/'tables/cluster_lineage_counts.tsv',sep='\t',encoding='utf-8')
dump(O/'PREPARATION_COMPLETE.json',{'status':'PASS','cells':len(cells),'genes':len(feat),'all_original_counts_preserved':True,'source_hashes':receipt['files'],'reference_cells':int(cells.reference.sum()),'heldout_cells':int(cells.heldout.sum()),'doublets':int(cells.full_scDblFinder_class.eq('doublet').sum()),'clusters':int(cells.cluster.nunique()),'independent_source_patients':1,'new_formal_CAF_patients':0,'new_inferential_tests':False})
print('PREPARATION_COMPLETE',flush=True)
