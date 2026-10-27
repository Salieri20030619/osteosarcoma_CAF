from pathlib import Path
import os,sys,gc,json,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase63_matrix_resistance_20261005';B=P/'phase58_F1F5_round1_20261004';A=P/'phase41_full_transcriptome_setup_20261001'
sys.path.insert(0,str(B/'code'))
from common import tests,fdr,mean_valid
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
assert (O/'FREEZE_COMPLETE').exists()
config=json.loads((O/'resources/PROGRAMS.json').read_text());sets=config['sets'];roles=config['roles'];genes=sorted(set().union(*map(set,sets.values())));gi={g:i for i,g in enumerate(genes)};mem={p:np.array([gi[g] for g in gs],int) for p,gs in sets.items()}
ob=read(P/'phase62_F1F5_internal_20261005/tables/F1F5_internal_cell_labels.tsv.gz').set_index('cell_id');assert len(ob)==8174 and ob.subgroup_id.nunique()==33
mp=read(B/'resources/source_canonical_feature_map.tsv');presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict();assay=np.zeros((len(ob),len(genes)),bool);counts={};den={};featureledger=[]
for mode,path in [('raw',A/'objects/all_candidate_full_genes.h5ad'),('corrected',B/'objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==ob.index.tolist();v=np.zeros((len(ob),len(genes)),float);den[mode]=np.asarray(a.X.sum(1)).ravel();assert (den[mode]>0).all()
 for ds in sorted(ob.dataset.unique()):
  rr=np.flatnonzero(ob.dataset.eq(ds));mapping=mp[mp.dataset.eq(ds)].set_index('symbol').feature.to_dict();ix=a.var_names.get_indexer([mapping.get(g,'__MISSING__') for g in genes]);ok=ix>=0;v[np.ix_(rr,np.flatnonzero(ok))]=a.X[rr][:,ix[ok]].toarray()
  if mode=='raw':
   for sm in ob.iloc[rr]['sample'].unique():
    ss=np.flatnonzero(ob['sample'].eq(sm));assay[np.ix_(ss,np.flatnonzero(ok))]=presence[sidx[sm],ix[ok]].toarray().ravel()>0
   featureledger.extend({'dataset':ds,'gene':g,'feature':mapping.get(g,''),'assayed_samples':int(sum(assay[ob.index.get_indexer(ob[ob['sample'].eq(s)].index)[0],j] for s in ob.iloc[rr]['sample'].unique())),'samples':ob.iloc[rr]['sample'].nunique()} for j,g in enumerate(genes))
 counts[mode]=v;del a;gc.collect()
np.savez_compressed(O/'objects/mechanism_counts.npz',raw=counts['raw'],corrected=counts['corrected'],assayed=assay,raw_library=den['raw'],corrected_library=den['corrected'],genes=np.array(genes),cells=ob.index.to_numpy(str))
read(P/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv').to_csv(O/'tables/subgroup_baseline.tsv',sep='\t',index=False)
pd.DataFrame(featureledger).to_csv(O/'tables/mechanism_feature_coverage.tsv',sep='\t',index=False)
cell=ob.copy();profiles=[];coverage=[]
for mode in counts:
 log=np.log1p(counts[mode]/den[mode][:,None]*1e4)
 for name,ii in mem.items():
  cell[mode+'__'+name]=np.nan
  for ds in sorted(ob.dataset.unique()):
   rr=np.flatnonzero(ob.dataset.eq(ds));jj=ii[assay[rr][:,ii].all(0)];ok=len(jj)>=min(3,len(ii)) and len(jj)/len(ii)>=.8
   coverage.append({'dataset':ds,'matrix':mode,'program':name,'measured':len(jj),'defined':len(ii),'score_evaluable':ok})
   if ok:cell.iloc[rr,cell.columns.get_loc(mode+'__'+name)]=log[np.ix_(rr,jj)].mean(1)
  for (ds,st),g in cell.groupby(['dataset','subgroup_id']):
   values=g[mode+'__'+name].to_numpy(float)
   profiles.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'program':name,'cells':len(g),'donors':g.donor_record.nunique(),'mean_cell_score':np.nanmean(values) if np.isfinite(values).any() else np.nan,'median_cell_score':np.nanmedian(values) if np.isfinite(values).any() else np.nan,'evidence':'descriptive_mean_log1p_full_library_CP10k'})
cell.to_csv(O/'tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t')
pd.DataFrame(profiles).to_csv(O/'tables/descriptive_subgroup_scores.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/cell_score_coverage.tsv',sep='\t',index=False)
effects=[];stats=[];genestats=[];roster=[];genepat=[];lineageshares=[]
scopes=[('raw','all',1),('corrected','all',1),('corrected','all',5),('corrected','Fibro',1)]
for ds in sorted(ob.dataset.unique()):
 source=ob.dataset.eq(ds).to_numpy();stlist=sorted(ob.loc[source,'subgroup_id'].unique());bydonor=list(ob.loc[source].groupby('donor_record',sort=True))
 for st in stlist:
  for mode,scope,minimum in scopes:
   l=source&ob.subgroup_id.eq(st).to_numpy();r=source&ob.subgroup_id.ne(st).to_numpy()
   if scope=='Fibro':
    fib=ob.lineage_stratum.eq('Fibro').to_numpy();l &=fib;r &=fib
   if not l.any() or not r.any():continue
   de=[];pe=[];donors=[]
   for donor,dg in bydonor:
    dsitelist=[];psitelist=[]
    for site,sg in dg.groupby('site',dropna=False):
     ix=ob.index.get_indexer(sg.index);il=ix[l[ix]];ir=ix[r[ix]]
     if min(len(il),len(ir))<minimum:continue
     x=counts[mode];lib=den[mode];present=assay[np.r_[il,ir]].all(0);dv=np.log2(x[il].sum(0)/lib[il].sum()*1e4+.1)-np.log2(x[ir].sum(0)/lib[ir].sum()*1e4+.1);dv[~present]=np.nan;dsitelist.append(dv)
     pv=[]
     for name,jj in mem.items():
      measured=jj[present[jj]];ok=len(measured)>=min(3,len(jj)) and len(measured)/len(jj)>=.8;pv.append(float(dv[measured].mean()) if ok else np.nan)
     psitelist.append(pv);roster.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'donor_record':donor,'site':site,'target_cells':len(il),'reference_cells':len(ir),'assayed_genes':int(present.sum())})
    if dsitelist:de.append(mean_valid(np.array(dsitelist),axis=0));pe.append(mean_valid(np.array(psitelist),axis=0));donors.append(donor)
   gd=np.array(de).reshape((-1,len(genes)));pdiff=np.array(pe).reshape((-1,len(mem)));pvals,nn,method=tests(pdiff)
   for j,name in enumerate(mem):
    v=pdiff[:,j];valid=np.isfinite(v);vec=v[valid];fam=f'{ds}|{mode}|{scope}|min{minimum}|{roles[name]}'
    stats.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'program':name,'N':int(nn[j]),'mean_effect':float(vec.mean()) if len(vec) else np.nan,'median_effect':float(np.median(vec)) if len(vec) else np.nan,'positive_patients':int((vec>0).sum()),'negative_patients':int((vec<0).sum()),'P':pvals[j],'test':method[j],'family':fam,'evidence':'patient_paired_internal8174_mechanism_contrast'})
    for donor,value in zip(donors,v):effects.append({'dataset':ds,'subgroup_id':st,'matrix':mode,'scope':scope,'minimum_cells':minimum,'program':name,'donor_record':donor,'effect':value})
   if scope=='all' and minimum==1:
    gp,gn,gmeth=tests(gd);fr=pd.DataFrame({'gene':genes,'dataset':ds,'subgroup_id':st,'matrix':mode,'N':gn,'mean_effect':mean_valid(gd,axis=0),'positive_patients':(gd>0).sum(0),'negative_patients':(gd<0).sum(0),'P':gp,'test':gmeth,'family':f'{ds}|{mode}|all_subgroup_targeted_mechanism_genes'});genestats.append(fr)
    np.savez_compressed(O/'objects'/f'{ds}__{st}__{mode}__patient_gene_effects.npz',effects=gd,genes=np.array(genes),donors=np.array(donors,str))
  print('MECHANISM',ds,st,flush=True)
 # RNA contribution is descriptive, inside the original pool, never whole-tissue source purity.
 for donor,dg in bydonor:
  rr=ob.index.get_indexer(dg.index)
  for mode in ['raw','corrected']:
   for name,jj in mem.items():
    jj=jj[assay[rr][:,jj].all(0)]
    if not len(jj):continue
    denominator=float(counts[mode][np.ix_(rr,jj)].sum())
    for (st,lineage),gg in dg.groupby(['subgroup_id','lineage_stratum']):
     ii=ob.index.get_indexer(gg.index);num=float(counts[mode][np.ix_(ii,jj)].sum());lineageshares.append({'dataset':ds,'donor_record':donor,'matrix':mode,'program':name,'subgroup_id':st,'lineage_stratum':lineage,'cells':len(ii),'RNA_counts':num,'original_pool_RNA_share':num/denominator if denominator>0 else np.nan,'denominator_original_pool_RNA':denominator})
s=fdr(pd.DataFrame(stats));s.to_csv(O/'tables/patient_program_statistics.tsv.gz',sep='\t',index=False)
g=fdr(pd.concat(genestats,ignore_index=True));g.to_csv(O/'tables/targeted_mechanism_gene_statistics.tsv.gz',sep='\t',index=False);g[g.nominal].to_csv(O/'tables/targeted_nominal_genes.tsv.gz',sep='\t',index=False)
pd.DataFrame(effects).to_csv(O/'tables/patient_program_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(roster).to_csv(O/'tables/paired_patient_site_roster.tsv.gz',sep='\t',index=False);pd.DataFrame(lineageshares).to_csv(O/'tables/within_pool_RNA_contribution.tsv.gz',sep='\t',index=False)
(O/'MECHANISM_COMPLETE').write_text('SUCCESS\n');print('MECHANISM_COMPLETE',len(s),len(g),len(genes),flush=True)
