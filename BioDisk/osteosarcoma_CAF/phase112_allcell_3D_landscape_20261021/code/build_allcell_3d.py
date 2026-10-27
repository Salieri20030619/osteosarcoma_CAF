from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,traceback,sys,gc
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
from sklearn.manifold import trustworthiness
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase112_allcell_3D_landscape_20261021'
for s in ['code','tables','objects','resources','logs','delivery']:(O/s).mkdir(parents=True,exist_ok=True)
def dump(name,value):(O/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def stage(s):
 dump('STATUS.json',dict(status=s,updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat()));print(s,flush=True)
try:
 stage('loading_frozen_annotations')
 ann=R/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz';clin=R/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv'
 cols=['original_cell_id','sample','dataset','original_donor_record','review_lineage_code','review_expression_lineage','lineage_confidence','mixed_program_flag']
 d=pd.read_csv(ann,sep='\t',usecols=cols);d=d[d.dataset.isin(['SCPCP000017','SCPCP000023'])].copy()
 cm=pd.read_csv(clin,sep='\t');d=d.merge(cm[['sample','donor_record','strict_site_group','full_called_cells']],on='sample',validate='many_to_one')
 assert d.original_cell_id.is_unique and (d.original_donor_record==d.donor_record).all()
 assert d.groupby('sample').size().equals(d.groupby('sample').full_called_cells.first().astype('int64'))
 mapping={'Chondro':'Chondrocyte-like','Endo':'Endothelial','Epi':'Epithelial','Ery':'Erythroid','Fibro':'Fibroblast-like','Lymph':'Lymphatic endothelium','Macro':'Macrophage','Mast':'Mast','Mesench':'Mesenchymal unresolved','MixBPlasma':'Mixed B / plasma','MixEndoEpi':'Mixed endothelial / epithelial','MonoDC':'Monocyte / DC unresolved','Mural':'Mural / pericyte','Myogenic':'Myogenic','OC':'Osteoclast-like','Osteo':'Osteogenic-like','Plasma':'Plasma-like','Renal':'Renal epithelial-like','T':'T cell','U':'Unresolved lineage','pDC':'pDC-like'}
 assert set(d.review_lineage_code)<=set(mapping)
 d['cell_type']=d.review_lineage_code.map(mapping);d=d.rename(columns={'original_cell_id':'cell_id','strict_site_group':'site'})
 d.to_csv(O/'tables/all_cell_metadata.tsv.gz',sep='\t',index=False)
 pd.DataFrame([dict(code=k,display=v) for k,v in mapping.items()]).to_csv(O/'resources/label_mapping.tsv',sep='\t',index=False)
 counts=d.groupby(['dataset','site','donor_record','cell_type']).size().rename('cells').reset_index()
 # Explicit zero categories for each patient-site, retaining all 21 labels.
 rows=[]
 for (src,site,donor),q in d.groupby(['dataset','site','donor_record']):
  vc=q.cell_type.value_counts()
  for label in sorted(mapping.values()):rows.append(dict(dataset=src,site=site,donor=donor,cell_type=label,cells=int(vc.get(label,0)),denominator=len(q),fraction_percent=float(vc.get(label,0)/len(q)*100)))
 comp=pd.DataFrame(rows);assert np.allclose(comp.groupby(['dataset','site','donor']).fraction_percent.sum(),100)
 comp.to_csv(O/'tables/patient_site_composition.tsv',sep='\t',index=False)
 means=comp.groupby(['dataset','site','cell_type']).agg(mean_patient_percent=('fraction_percent','mean'),patients=('donor','nunique'),pooled_cells=('cells','sum')).reset_index()
 means.to_csv(O/'tables/site_patient_balanced_composition.tsv',sep='\t',index=False)
 cohort=d.groupby(['dataset','site']).agg(cells=('cell_id','size'),patients=('donor_record','nunique'),samples=('sample','nunique')).reset_index();cohort.to_csv(O/'tables/cohort_counts.tsv',sep='\t',index=False)
 overlaps=[]
 for src,q in d.groupby('dataset'):
  both=set(q[q.site.eq('Bone')].donor_record)&set(q[q.site.eq('Lung')].donor_record)
  overlaps.append(dict(dataset=src,total_unique_donors=q.donor_record.nunique(),paired_Bone_Lung_donors=sorted(both)))
 dump('resources/PATIENT_OVERLAP.json',overlaps)
 paths=[ann,clin,R/'phase53_de_novo_reannotation_20261002/code/build_structure.py'];metrics=[]
 protocol=dict(scope='All upstream called cells in two original sources; no CAF-only denominator and no newly applied gate',annotation='Phase55 frozen reviewed lineages, including unresolved and mixed labels',clinical='Phase81 strict anatomical site; unknown/mixed sites retained as separate context',comparison='Descriptive only, no group inference; shared Bone/Lung donors retained and documented',composition='Per-patient/site fractions then equal patient weighting for site-average bars; missing lineages have zero fraction',embedding='Reuse complete Phase53 source-specific integrated 50D representation and saved neighborhood graph; optimize true 3D UMAP on all source cells. Bone/Lung views are subsets of fixed source coordinates.',parameters=dict(n_components=3,n_neighbors=15,min_dist=.5,spread=1.,epochs=200,random_state=20261002,init='spectral'),reference_2D='Frozen original Phase53 2D embedding with same representation/graph and parameter settings',inherited_gene_exclusion=['COL1A1','LOX','MMP14'],no_reclustering=True,no_label_reassignment=True,source_coordinates_not_comparable=True,NCI='Not included in this first all-cell landscape batch; only original sources with Bone/Lung cohorts',scanpy=sc.__version__)
 dump('resources/PROTOCOL.json',protocol)
 for src in ['SCPCP000017','SCPCP000023']:
  stage('embedding_'+src)
  base=R/'phase53_de_novo_reannotation_20261002/objects'/src;emb=base/'embeddings.npz';graph=base/'connectivities.npz';paths.extend([emb,graph]);z=np.load(emb)
  ids=pd.Index(z['cell_id']);q=d[d.dataset.eq(src)].set_index('cell_id').loc[ids].copy();assert len(q)==len(ids) and set(q.index)==set(d.loc[d.dataset.eq(src),'cell_id'])
  x=z['integrated'].astype(np.float32);g=sparse.load_npz(graph).tocsr();assert g.shape==(len(q),len(q)) and np.isfinite(g.data).all()
  a=ad.AnnData(sparse.csr_matrix((len(q),0)),obs=pd.DataFrame(index=ids));a.obsm['X_integrated']=x;a.obsp['connectivities']=g
  a.uns['neighbors']=dict(connectivities_key='connectivities',distances_key='distances',params=dict(n_neighbors=15,method='umap',metric='euclidean',use_rep='X_integrated'))
  out=O/'objects'/(src+'_UMAP3D.npz')
  if out.exists():
   cached=np.load(out);assert np.array_equal(cached['cell_id'],ids);u=cached['UMAP3D']
  else:
   sc.tl.umap(a,n_components=3,min_dist=.5,spread=1.,maxiter=200,init_pos='spectral',random_state=20261002)
   u=a.obsm['X_umap'];assert u.shape==(len(q),3) and np.isfinite(u).all()
   np.savez_compressed(out,cell_id=ids.to_numpy(str),UMAP3D=u,UMAP2D=z['UMAP'])
  for j in range(3):q['UMAP'+str(j+1)]=u[:,j]
  q['legacy_UMAP1']=z['UMAP'][:,0];q['legacy_UMAP2']=z['UMAP'][:,1]
  q.reset_index(names='cell_id').to_csv(O/'tables'/(src+'_plot_cells.tsv.gz'),sep='\t',index=False)
  stage('checking_'+src)
  ix=np.sort(np.random.default_rng(20261021).choice(len(q),2000,replace=False))
  # Sample-restricted diagnostic, not a whole-dataset neighborhood-fidelity estimate.
  for dim,uu in [(2,z['UMAP']),(3,u)]:metrics.append(dict(dataset=src,dimensions=dim,diagnostic='trustworthiness within a fixed random 2000-cell subset',sampled_cells=2000,k=15,value=float(trustworthiness(x[ix],uu[ix],n_neighbors=15))))
  del a,g,x,z,u,q;gc.collect()
 pd.DataFrame(metrics).to_csv(O/'tables/embedding_diagnostic.tsv',sep='\t',index=False)
 dump('resources/INPUT_HASHES.json',[dict(path=str(p),sha256=sha(p)) for p in paths])
 dump('COMPUTATION_COMPLETE.json',dict(status='PASS',cells=len(d),sources=2,annotations_unique=d.cell_id.is_unique,full_sample_denominators_verified=True,patient_fraction_sums_verified=True,source_3d_embeddings_complete=True,original_files_modified=False,figures_pending=True))
 stage('computation_complete')
except Exception:
 dump('FAILURE.json',dict(traceback=traceback.format_exc()));raise
