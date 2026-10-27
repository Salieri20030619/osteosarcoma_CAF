import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json,itertools,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'audit_20261002_process_results';O.mkdir(exist_ok=True)
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad');m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='	',keep_default_na=False)
genes=['COL1A1','LOX','MMP14','MRC2','FGFR2','CA12','RUNX2','ANGPT1','SMAD3','ATP1A1','BACE1','LAMP2'];cw=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='	',keep_default_na=False);cw['prio']=~cw.mapping_method.eq('HGNC_Ensembl_approved');cw=cw.sort_values(['prio','feature']).drop_duplicates('canonical_symbol').set_index('canonical_symbol');features=cw.loc[genes,'feature'].tolist();out=[];native={};samplechecks=[]
for sm in m.loc[m.dataset.eq('SCPCP000017'),'sample']:
 f=P/'phase26_client_reanalysis_20260922/inputs'/sm/'full_called_counts.h5ad';b=ad.read_h5ad(f,backed='r');want=a.obs_names[a.obs['sample'].astype(str).eq(sm)];idx=b.obs_names.get_indexer(want);ix=np.argsort(idx);want=want[ix];idx=idx[ix];j=b.var_names.get_indexer(features);assert (j>=0).all()
 if len(idx):
  z=b.X[idx,:].tocsr();lib=np.asarray(z.sum(axis=1)).ravel();counts=z[:,j].toarray();ai=a.obs_names.get_indexer(want);re=a.X[ai][:,a.var_names.get_indexer(features)].toarray();e=float(np.max(np.abs(re-counts)));le=float(np.max(np.abs(a.obs.iloc[ai].full_library_counts.to_numpy(float)-lib)));assert e==0 and le==0
  for cid,vals,l in zip(want,counts,lib):native[str(cid)]=(vals,l)
  samplechecks.append(dict(sample=sm,cells=len(want),selected_count_error=e,library_error=le))
 b.file.close()
ob=a.obs[a.obs.dataset.astype(str).eq('SCPCP000017')].copy()
for layer,cl in [('supported','F3'),('expanded','F4')]:
 sub=ob[ob.in_supported.astype(bool)] if layer=='supported' else ob
 pat=[];don=[]
 for donor,d in sub.groupby('donor_record',observed=True):
  site=[]
  for st,z in d.groupby('site',observed=True):
   masks=[z[layer+'_class_id'].astype(str).eq(cl),~z[layer+'_class_id'].astype(str).eq(cl)]
   if not all(x.any() for x in masks):continue
   vv=[]
   for mask in masks:
    pp=[native[str(i)] for i in z.index[mask]];tot=np.sum([r[0] for r in pp],axis=0);den=sum(r[1] for r in pp);vv.append(np.log2(tot/den*10000+.1))
   site.append(vv[0]-vv[1])
  if site:pat.append(np.mean(site,axis=0));don.append(str(donor))
 mat=np.asarray(pat);sg=np.asarray(list(itertools.product([-1,1],repeat=len(don))));pv=np.mean(np.abs(sg@mat/len(don))>=np.abs(mat.mean(0))-1e-12,axis=0)
 table=pd.read_csv(C/'tables'/f'{layer}__SCPCP000017__full_gene_statistics.tsv.gz',sep='	',low_memory=False);table=table[table.contrast.eq(cl+'_vs_rest')].set_index('feature')
 for k,(g,f) in enumerate(zip(genes,features)):
  r=table.loc[f];out.append(dict(layer=layer,gene=g,N=len(don),recomputed_mean=float(mat[:,k].mean()),stored_mean=float(r.mean_delta),recomputed_P=float(pv[k]),stored_P=float(r.P),max_mean_error=float(abs(mat[:,k].mean()-r.mean_delta)),P_error=float(abs(pv[k]-r.P))))
 pd.DataFrame(mat,index=don,columns=genes).to_csv(O/f'{layer}_native_patient_effects.tsv',sep='	',index_label='patient_id')
pd.DataFrame(out).to_csv(O/'raw_count_reconstruction.tsv',sep='	',index=False);pd.DataFrame(samplechecks).to_csv(O/'native_alignment_checks.tsv',sep='	',index=False)
result={'source':'SCPCP000017','native_samples_checked':len(samplechecks),'candidate_cells_checked':sum(x['cells'] for x in samplechecks),'genes':genes,'comparisons':len(out),'max_mean_error':max(x['max_mean_error'] for x in out),'max_P_error':max(x['P_error'] for x in out),'production_changed':False};(O/'RAW_RECONSTRUCTION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result));print('AUDIT_RAW_RECONSTRUCTION_COMPLETE')
