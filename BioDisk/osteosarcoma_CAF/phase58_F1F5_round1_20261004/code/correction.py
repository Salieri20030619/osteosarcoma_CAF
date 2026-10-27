from common import *
import anndata as ad,gzip,gc
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
assert (O/'FREEZE_COMPLETE').exists()
raw=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad');ob=read(O/'tables/frozen_expanded_cells.tsv.gz').set_index('cell_id').loc[raw.obs_names];features=raw.var_names.astype(str).to_numpy();pieces=[];ids=[];libs=[];qc=[]
for i,(sample,g) in enumerate(ob.groupby('sample',sort=True,observed=True)):
 base=P/'phase26_client_reanalysis_20260922/inputs'/sample;de=P/'phase31_customer_completion_20260926/decontx'/sample;assert (de/'COMPLETE').exists() and (de/'exit_code').read_text().strip()=='0'
 cells=read(base/'cells.tsv');genes=read(base/'genes.tsv');cont=read(de/'cell_contamination.tsv');assert cells.cell_id.tolist()==cont.cell_id.tolist()
 with gzip.open(de/'corrected_counts.mtx.gz','rb') as h:
  with threadpool_limits(limits=4):x=mmread(h).tocsr().T.tocsr()
 ix=pd.Index(cells.cell_id).get_indexer(g.index);assert (ix>=0).all();jx=pd.Index(features).get_indexer(genes.ensembl);assert (jx>=0).all();part=x[ix].tocoo();full=sparse.csr_matrix((part.data,(part.row,jx[part.col])),shape=(len(g),len(features)));den=cont.corrected_UMI.to_numpy()[ix];assert np.allclose(np.asarray(full.sum(1)).ravel(),den,rtol=1e-5,atol=1e-4)
 pieces.append(full);ids.extend(g.index);libs.extend(den);qc.append({'sample':sample,'candidate_cells':len(g),'corrected_library_verified':True});del x,full,part;gc.collect();print('CORRECTED_SAMPLE',i+1,sample,len(g),flush=True)
 allx=None
x=sparse.vstack(pieces,format='csr');order=pd.Index(ids).get_indexer(ob.index);assert (order>=0).all();x=x[order];den=np.asarray(libs)[order];cb=ob.copy();cb['corrected_full_library_counts']=den;corrected=ad.AnnData(x,obs=cb,var=raw.var.copy());corrected.write_h5ad(O/'objects/expanded_corrected_full_genes.h5ad',compression='gzip');pd.DataFrame(qc).to_csv(O/'tables/correction_input_QA.tsv',sep='\t',index=False)
del pieces,raw;gc.collect();presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict();manifest=read(A/'tables/clinical_sample_registry.tsv');cw=read(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',keep_default_na=False).set_index('feature').loc[features];symbols=cw.canonical_symbol.to_numpy();quality=cw.mapping_method.ne('HGNC_Ensembl_approved').astype(int).to_numpy();sets=json.loads((O/'resources/programs_frozen.json').read_text())['sets'];programs=list(sets);allprog=[];summary=[]
for ds in manifest.dataset.unique():
 samples=manifest.loc[manifest.dataset.eq(ds),'sample'].tolist();gi=np.flatnonzero(np.asarray(presence[[sidx[s] for s in samples]].sum(0)).ravel()>0);fe=features[gi];sy=symbols[gi];mapping={}
 for j,(f,s,q) in enumerate(zip(fe,sy,quality[gi])):
  if s and (s not in mapping or (q,f)<mapping[s][0]):mapping[s]=((q,f),j)
 mem={name:[mapping[s][1] for s in members if s in mapping] for name,members in sets.items()};sg=cb[cb.dataset.eq(ds)];gf=[];pf=[]
 for cls in ['F1','F2','F3','F4','F5']:
  rows=[];donors=[]
  for donor,dg in sg.groupby('donor_record',observed=True):
   vv=[]
   for site,g in dg.groupby('site',observed=True):
    left=g.class_id.eq(cls).to_numpy();right=~left
    if not left.any() or not right.any():continue
    il=cb.index.get_indexer(g.index[left]);ir=cb.index.get_indexer(g.index[right]);ll=den[il].sum();rr=den[ir].sum()
    if ll<=0 or rr<=0:continue
    lc=np.asarray(x[il][:,gi].sum(0)).ravel();rc=np.asarray(x[ir][:,gi].sum(0)).ravel();required=g['sample'].unique();valid=np.asarray(presence[[sidx[s] for s in required]][:,gi].sum(0)).ravel()==len(required);dd=np.log2(lc/ll*1e4+.1)-np.log2(rc/rr*1e4+.1);dd[~valid]=np.nan;vv.append(dd)
   if vv:rows.append(mean_valid(np.vstack(vv)));donors.append(str(donor))
  d=np.vstack(rows) if rows else np.empty((0,len(gi)));pv,N,method=tests(d);context=f'expanded__{ds}__{cls}_vs_rest';np.savez_compressed(O/'objects'/f'{context}__corrected_patient_effects.npz',effects=d,feature_indices=gi,donors=np.asarray(donors,dtype=str))
  gf.append(pd.DataFrame({'dataset':ds,'class_id':cls,'contrast':cls+'_vs_rest','feature':fe,'symbol':sy,'N':N,'mean_delta':mean_valid(d),'positive_patients':np.sum(d>0,axis=0),'negative_patients':np.sum(d<0,axis=0),'P':pv,'test':method,'family':f'expanded|{ds}|corrected_all_F1_F5_genes'}))
  score=np.column_stack([mean_valid(d[:,ix],axis=1) if ix else np.full(len(donors),np.nan) for ix in mem.values()]);pN,nn,mm=tests(score);pf.append(pd.DataFrame({'dataset':ds,'class_id':cls,'contrast':cls+'_vs_rest','program':programs,'N':nn,'mean_delta':mean_valid(score),'positive_patients':np.sum(score>0,axis=0),'negative_patients':np.sum(score<0,axis=0),'P':pN,'test':mm,'family':f'expanded|{ds}|corrected_all_F1_F5_programs'}));np.savez_compressed(O/'objects'/f'{context}__corrected_program_effects.npz',effects=score,programs=np.asarray(programs),donors=np.asarray(donors,dtype=str));print('CORRECTED_STATS',ds,cls,len(donors),flush=True)
 full=fdr(pd.concat(gf,ignore_index=True));full.to_csv(O/'tables/genes'/f'expanded__{ds}__corrected_statistics.tsv.gz',sep='\t',index=False);pr=fdr(pd.concat(pf,ignore_index=True));allprog.append(pr)
pd.concat(allprog,ignore_index=True).to_csv(O/'tables/corrected_program_statistics.tsv.gz',sep='\t',index=False);dump(O/'delivery/CORRECTION_ACCEPTANCE.json',{'cells':len(cb),'samples_with_candidates':len(qc),'raw_source_called_cells_not_recomputed':True,'DecontX_reused':True,'no_new_cell_exclusions':True,'full_feature_presence_retained':True});(O/'CORRECTION_COMPLETE').write_text('SUCCESS\n');print('CORRECTION_COMPLETE',flush=True)