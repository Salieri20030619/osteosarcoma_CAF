from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
p=Path(r'C:\Users\31937\CAF_F2C3_novelty_20261009\deep_review\results');t=p/'tables'
c=pd.read_csv(t/'cell_signature_scores.tsv.gz',sep='\t');c=c[(c.dataset=='SCPCP000017')&(c.lineage_stratum=='Fibro')&(~c.doublet_positive_any)].copy();rows=[];models=[];stats=[]
for mode,d in c.groupby('matrix'):
 targets=d.loc[d.subgroup_id.eq('F2_C3'),'donor_record'].unique()
 for nuisance in ['Zheng_mCAF6','Atlas_mCAF4']:
  for donor in targets:
   train=d[(d.donor_record!=donor)&(d.subgroup_id!='F2_C3')].copy();hold=d[d.donor_record==donor].copy();assert not train.donor_record.eq(donor).any() and not train.subgroup_id.eq('F2_C3').any()
   X=np.column_stack([np.ones(len(train)),train[nuisance],np.log1p(train.library)]);Y=train.F2C3_processing2.to_numpy();weights=1/train.groupby('donor_record').donor_record.transform('size').to_numpy();root=np.sqrt(weights);coef=np.linalg.lstsq(X*root[:,None],Y*root,rcond=None)[0]
   hold['residual']=hold.F2C3_processing2-np.column_stack([np.ones(len(hold)),hold[nuisance],np.log1p(hold.library)])@coef
   models.append(dict(matrix=mode,nuisance=nuisance,heldout_donor=donor,training_donors=train.donor_record.nunique(),training_cells=len(train),intercept=coef[0],beta_mCAF=coef[1],beta_log_library=coef[2]))
   for minimum in [1,5]:
    for refname in ['other_Fibro','other_F2']:
     diffs=[]
     for site,sg in hold.groupby('site'):
      target=sg[sg.subgroup_id=='F2_C3'];ref=sg[sg.subgroup_id!='F2_C3']
      if refname=='other_F2':ref=ref[ref.class_id=='F2']
      if min(len(target),len(ref))>=minimum:diffs.append(target.residual.mean()-ref.residual.mean())
     if diffs:rows.append(dict(matrix=mode,nuisance=nuisance,min_cells=minimum,comparison=refname,donor_record=donor,effect=np.mean(diffs)))
e=pd.DataFrame(rows)
for key,g in e.groupby(['matrix','nuisance','min_cells','comparison']):
 n=len(g);v=g.effect.to_numpy();pv=np.nan
 if n>=2 and np.any(v!=0):pv=float((np.abs(np.array(list(itertools.product([-1.,1.],repeat=n)))@v/n)>=abs(v.mean())-1e-12).mean())
 stats.append(dict(zip(['matrix','nuisance','min_cells','comparison'],key),N=n,mean_effect=v.mean(),positive=int((v>0).sum()),P=pv,LOPO_min=min((v.sum()-z)/(n-1) for z in v) if n>1 else np.nan))
s=pd.DataFrame(stats);s['family']='posthoc_residual_processing_all_models_contrasts_matrices_thresholds';s['FDR']=np.nan;ix=s.P.notna();s.loc[ix,'FDR']=multipletests(s.loc[ix,'P'],method='fdr_bh')[1]
e.to_csv(t/'residual_processing_patient_effects.tsv',sep='\t',index=False);s.to_csv(t/'residual_processing_statistics.tsv',sep='\t',index=False);pd.DataFrame(models).to_csv(t/'residual_processing_LODO_models.tsv',sep='\t',index=False)
protocol={'purpose':'Posthoc targeted check of processing2 signal beyond published matrix signature and library depth','model':'PCOLCE-LOXL2 score ~ intercept + mCAF score + log1p(full_library_counts)','training':'Other Fibro only, no F2_C3; held-out donor entirely excluded; equal total weight per training donor','inference':'Two-sided exact signflip on source-local donor/site residual differences; no cell-level or regression-coefficient P values; entire 16 residual comparisons one BH family','limitations':'Exploratory linear nuisance adjustment, conditional on captured candidate cells and available panels; not proof of novel subtype or causal mechanism'}
(p/'resources/RESIDUAL_PROTOCOL.json').write_text(json.dumps(protocol,indent=2));print(s.round(5).to_string(index=False))
d=pd.read_csv(t/'F2C3_donor_gene_detection.tsv',sep='\t');m=d[d.gene=='MMP9'].copy();m['count_total']=m.cells*m.mean_corrected_count;print('MMP9 dominant fraction',m.count_total.max()/m.count_total.sum())
