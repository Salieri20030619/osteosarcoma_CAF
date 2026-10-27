from pathlib import Path
import sys,os,json,itertools
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,anndata as ad
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'audit_myCAF_report_20261010';O.mkdir(exist_ok=True)
sys.path.insert(0,str(P/'phase66_doublet_repair_DSC1_20261005/code'));from current_metadata import load_current_cells
m=load_current_cells().set_index('cell_id');genes=['ACTA2','TAGLN','MYH11','MYL9','TPM2','CNN1'];mp=pd.read_csv(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');effects=[];stats=[];dets=[];scoreout=[]
for mode,path in [('raw',P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad'),('corrected',P/'phase58_F1F5_round1_20261004/objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==m.index.tolist();rr=np.flatnonzero(m.dataset.eq('SCPCP000017'));md=m.iloc[rr].copy();mapping=mp[(mp.dataset=='SCPCP000017')&mp.symbol.isin(genes)].set_index('symbol').loc[genes];jj=a.var_names.get_indexer(mapping.feature);x=a.X[rr][:,jj].toarray();lib=np.asarray(a.X.sum(1)).ravel()[rr];v=np.log1p(x/lib[:,None]*1e4);md['exact_myCAF6']=v.mean(1)
 md[['dataset','sample','donor_record','site','class_id','subgroup_id','lineage_stratum','candidate_tier','CNV_evidence_current','doublet_positive_any','program_myCAF_customer','exact_myCAF6']].assign(matrix=mode).reset_index().to_csv(O/f'cell_six_gene_scores_{mode}.tsv',sep='\t',index=False)
 for scope in ['frozen_all','Fibro_no_doublets']:
  allowed=np.ones(len(md),bool) if scope=='frozen_all' else md.lineage_stratum.eq('Fibro').to_numpy()&~md.doublet_positive_any.to_numpy(bool)
  base=md.iloc[np.flatnonzero(allowed)]
  for target in ['F2','F2_C0','F2_C1','F2_C2','F2_C3','F2_C4','F2_C5']:
   for refname in (['other_candidates'] if target=='F2' else ['other_F2','other_candidates']):
    ds=[]
    for donor,dg in base.groupby('donor_record'):
     siteeff=[]
     for site,sg in dg.groupby('site'):
      left=sg[sg.class_id.eq('F2')] if target=='F2' else sg[sg.subgroup_id.eq(target)]
      right=sg[~sg.class_id.eq('F2')] if target=='F2' else sg[~sg.subgroup_id.eq(target)]
      if refname=='other_F2':right=right[right.class_id.eq('F2')]
      if len(left) and len(right):siteeff.append(left.exact_myCAF6.mean()-right.exact_myCAF6.mean())
     if siteeff:
      e=np.mean(siteeff);ds.append(e);effects.append(dict(matrix=mode,scope=scope,target=target,reference=refname,donor_record=donor,effect=e))
    n=len(ds);pv=np.nan
    if n>=2 and np.any(np.array(ds)!=0):pv=float((np.abs(np.array(list(itertools.product([-1.,1.],repeat=n)))@ds/n)>=abs(np.mean(ds))-1e-12).mean())
    stats.append(dict(matrix=mode,scope=scope,target=target,reference=refname,N=n,mean_effect=np.mean(ds) if n else np.nan,median_effect=np.median(ds) if n else np.nan,positive=sum(z>0 for z in ds),negative=sum(z<0 for z in ds),P=pv,family=f'{mode}|{scope}|six_gene_identity_audit_all_targets_references'))
   mask=(md.class_id.eq('F2') if target=='F2' else md.subgroup_id.eq(target)).to_numpy()&allowed;ids=np.flatnonzero(mask)
   for j,g in enumerate(genes):dets.append(dict(matrix=mode,scope=scope,target=target,gene=g,cells=len(ids),donors=md.iloc[ids].donor_record.nunique(),detected=int((x[ids,j]>=(1 if mode=='raw' else .5)).sum()),mean_logCP10k=float(v[ids,j].mean()) if len(ids) else np.nan))
 if mode=='raw':print('ORIGINAL_SCORE_DIFF',float(np.max(np.abs(md.exact_myCAF6-md.program_myCAF_customer))))
s=pd.DataFrame(stats);s['FDR']=np.nan
for fam,g in s[s.P.notna()].groupby('family'):s.loc[g.index,'FDR']=multipletests(g.P,method='fdr_bh')[1]
s.to_csv(O/'exact_six_gene_patient_statistics.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(O/'exact_six_gene_patient_effects.tsv',sep='\t',index=False);pd.DataFrame(dets).to_csv(O/'exact_six_gene_detection.tsv',sep='\t',index=False)
m[m.class_id.eq('F2')].groupby(['subgroup_id','dataset','lineage_stratum','candidate_tier','CNV_evidence_current'],dropna=False).size().rename('cells').reset_index().to_csv(O/'F2_current_identity_composition.tsv',sep='\t',index=False)
(O/'COMPLETE.json').write_text(json.dumps({'status':'completed','panel':genes,'source':'SCPCP000017','scope':'Read-only audit recomputation, no label changes','mean_logCP10k_score':True,'source_files_unchanged':True,'new_FDR_families_not_substitute_historical':True},indent=2));print(s[(s.matrix=='raw')&s.target.isin(['F2','F2_C1','F2_C3'])].round(6).to_string(index=False))
