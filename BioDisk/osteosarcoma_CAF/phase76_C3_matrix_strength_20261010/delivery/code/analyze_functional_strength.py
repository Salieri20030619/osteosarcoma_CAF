from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
import json,sys,itertools
import numpy as np,pandas as pd,anndata as ad
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase76_C3_matrix_strength_20261010';O.mkdir(exist_ok=True)
sys.path.insert(0,str(P/'phase66_doublet_repair_DSC1_20261005/code'));from current_metadata import load_current_cells
m=load_current_cells().set_index('cell_id');mp=pd.read_csv(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');old=pd.read_csv(P/'phase75_F2C3_novelty_20261009/tables/cell_signature_scores.tsv.gz',sep='\t');refs=['other_Fibro','other_F2','F2_C0','F2_C1','F3_C1','F4_C8','F4_C10'];panels={'Published_COL5A1_MSC4':['COL5A1','COL16A1','EMILIN1','ADAMTS2'],'C3_assembly5':['COL3A1','COL5A1','COL5A2','PCOLCE','LOXL2'],'processing2':['PCOLCE','LOXL2']};genes=sorted(set(sum(panels.values(),[]))|{'COL1A1','COL1A2','COL3A1','COL4A1','COL4A2','FN1','POSTN','VCAN','COL11A1','LOX','PXDN'})
protocol={'scope':'C3 collagen/matrix strength and functional novelty; no myCAF identity reclassification','source':'SCPCP000017; captured frozen-pool Fibro cells, current doublets excluded','new_reference':'Zheng et al. J Cancer Res Clin Oncol 2024, DOI10.1007/s00432-024-05787-2; four genes explicitly named in Results/FigS3F, not a complete published signature','panels':panels,'comparison':'same donor/site target vs stated comparator, site means equally averaged per donor; min1 and min5 per arm','tests':'Two-sided exact donor signflip N>=2 and nonzero; raw P<.05 nominal; all 3 panels x 7 comparisons grouped by matrix/min-cell threshold for BH; residual comparisons form one separate family','residual':'processing2 ~ published COL5A1_MSC4 + log1p(library); train on non-C3 Fibro of OTHER donors, equal total weight per training donor; compare residuals within held-out donor/site','ranking':'donor-equal descriptive mean across subgroups, no cross-source pooling; a score built from previously selected C3 genes is not independent validation','contribution':'RNA count share in captured Fibro candidates only; not tissue abundance, protein secretion or matrix deposition','source_fields_unchanged':True}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2));effects=[];stats=[];res_effects=[];res_stats=[];contributions=[];rankrows=[];models=[];scoresout=[]
def test(v):
 v=np.asarray(v,float);n=len(v)
 if n<2 or not np.any(v!=0):return np.nan
 return float((abs(np.array(list(itertools.product([-1.,1.],repeat=n)))@v/n)>=abs(v.mean())-1e-12).mean())
def compare(d,value,target='F2_C3',minimum=1,refname='other_Fibro'):
 vals={}
 for donor,dg in d.groupby('donor_record'):
  es=[]
  for site,z in dg.groupby('site'):
   left=z[z.subgroup_id==target];right=z[z.subgroup_id!=target]
   if refname=='other_F2':right=right[right.class_id=='F2']
   elif refname!='other_Fibro':right=right[right.subgroup_id==refname]
   if min(len(left),len(right))>=minimum:es.append(left[value].mean()-right[value].mean())
  if es:vals[donor]=float(np.mean(es))
 return vals
for mode,path in [('raw',P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad'),('corrected',P/'phase58_F1F5_round1_20261004/objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==m.index.tolist();ii=np.flatnonzero(m.dataset.eq('SCPCP000017')&m.lineage_stratum.eq('Fibro')&~m.doublet_positive_any);meta=m.iloc[ii].copy();maps=mp[(mp.dataset=='SCPCP000017')&mp.symbol.isin(genes)].set_index('symbol').loc[genes];jj=a.var_names.get_indexer(maps.feature);x=a.X[ii][:,jj].toarray();lib=np.asarray(a.X.sum(1)).ravel()[ii];y=np.log1p(x/lib[:,None]*1e4);meta['library']=lib
 for name,gs in panels.items():meta[name]=y[:,[genes.index(g) for g in gs]].mean(1)
 orig=old[(old.matrix==mode)].set_index('cell_id').loc[meta.index];assert np.allclose(meta.C3_assembly5,orig.F2C3_assembly5) and np.allclose(meta.processing2,orig.F2C3_processing2)
 scoresout.append(meta[['dataset','sample','donor_record','site','class_id','subgroup_id','library']+list(panels)].assign(matrix=mode).reset_index())
 for st,d in meta.groupby('subgroup_id'):
  for panel in panels:rankrows.append(dict(matrix=mode,subgroup=st,panel=panel,cells=len(d),donors=d.donor_record.nunique(),donor_equal_mean=d.groupby('donor_record')[panel].mean().mean()))
 for j,g in enumerate(genes):
  mask=meta.subgroup_id.eq('F2_C3').to_numpy();contributions.append(dict(matrix=mode,gene=g,target_cells=int(mask.sum()),background_cells=len(meta),target_count_sum=x[mask,j].sum(),total_count_sum=x[:,j].sum(),target_RNA_fraction=x[mask,j].sum()/x[:,j].sum() if x[:,j].sum()>0 else np.nan,cell_fraction=mask.mean()))
 for minimum in [1,5]:
  for ref in refs:
   for panel in panels:
    ds=compare(meta,panel,minimum=minimum,refname=ref);v=list(ds.values());n=len(v);stats.append(dict(matrix=mode,min_cells=minimum,reference=ref,panel=panel,N=n,mean_effect=np.mean(v) if n else np.nan,positive=sum(z>0 for z in v),P=test(v),family=f'{mode}|min{minimum}|all_functional_panels_refs'))
    for donor,val in ds.items():effects.append(dict(matrix=mode,min_cells=minimum,reference=ref,panel=panel,donor_record=donor,effect=val))
 residuals=[]
 for donor in meta.loc[meta.subgroup_id.eq('F2_C3'),'donor_record'].unique():
  tr=meta[(meta.donor_record!=donor)&~meta.subgroup_id.eq('F2_C3')];ho=meta[meta.donor_record==donor].copy();w=1/tr.groupby('donor_record').donor_record.transform('size').to_numpy();X=np.column_stack([np.ones(len(tr)),tr.Published_COL5A1_MSC4,np.log1p(tr.library)]);coef=np.linalg.lstsq(X*np.sqrt(w)[:,None],tr.processing2*np.sqrt(w),rcond=None)[0];ho['residual']=ho.processing2-np.column_stack([np.ones(len(ho)),ho.Published_COL5A1_MSC4,np.log1p(ho.library)])@coef;residuals.append(ho);models.append(dict(matrix=mode,heldout_donor=donor,training_donors=tr.donor_record.nunique(),training_cells=len(tr),intercept=coef[0],beta_COL5_program=coef[1],beta_library=coef[2]))
 rd=pd.concat(residuals)
 for minimum in [1,5]:
  for ref in refs:
   ds=compare(rd,'residual',minimum=minimum,refname=ref);v=list(ds.values());n=len(v);res_stats.append(dict(matrix=mode,min_cells=minimum,reference=ref,N=n,mean_effect=np.mean(v) if n else np.nan,positive=sum(z>0 for z in v),P=test(v),family='COL5_reference_residual_all_matrices_refs_thresholds'))
   for donor,val in ds.items():res_effects.append(dict(matrix=mode,min_cells=minimum,reference=ref,donor_record=donor,effect=val))
for name,records in [('functional_patient_statistics',stats),('residual_processing_statistics',res_stats)]:
 d=pd.DataFrame(records);d['FDR']=np.nan
 for fam,g in d[d.P.notna()].groupby('family'):d.loc[g.index,'FDR']=multipletests(g.P,method='fdr_bh')[1]
 d.to_csv(O/f'{name}.tsv',sep='\t',index=False)
for name,records in [('patient_effects',effects),('residual_patient_effects',res_effects),('captured_RNA_contribution',contributions),('descriptive_subgroup_ranks',rankrows),('LODO_models',models)]:pd.DataFrame(records).to_csv(O/f'{name}.tsv',sep='\t',index=False)
pd.concat(scoresout).to_csv(O/'cell_functional_scores.tsv.gz',sep='\t',index=False)
(O/'COMPLETE.json').write_text(json.dumps({'status':'completed','cells_per_matrix':len(meta),'target_cells':int(meta.subgroup_id.eq('F2_C3').sum()),'functional_records':len(stats),'residual_records':len(res_stats),'prior_two_panel_scores_reproduced':True},indent=2));print('COMPLETE',O)
