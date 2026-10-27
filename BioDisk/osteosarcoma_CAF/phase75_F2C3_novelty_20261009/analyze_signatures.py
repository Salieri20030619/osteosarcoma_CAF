from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,json,hashlib,itertools,datetime
import pandas as pd,numpy as np,anndata as ad
from scipy import sparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase75_F2C3_novelty_20261009';O.mkdir(exist_ok=True)
for d in ['tables','resources','figures']:(O/d).mkdir(exist_ok=True)
sys.path.insert(0,str(P/'phase66_doublet_repair_DSC1_20261005/code'))
from current_metadata import load_current_cells
panels={
'F2C3_identity3':['CD248','POSTN','VCAN'],
'F2C3_assembly5':['COL3A1','COL5A1','COL5A2','PCOLCE','LOXL2'],
'F2C3_processing2':['PCOLCE','LOXL2'],
'Zheng_mCAF6':['FAP','MMP11','COL1A1','COL11A1','POSTN','FN1'],
'Zheng_mCAF_nonstructural2':['FAP','MMP11'],
'Zheng_myoCAF7':['ACTA2','MYH11','TAGLN','ACTG2','MYL9','TPM1','TPM2'],
'Zheng_iCAF12':['VCAN','COMP','CXCL14','ADAMDEC1','MFAP5','CFD','APOD','C3','IGF1','ITM2A','PLA2G2A','FBLN2'],
'Atlas_mCAF4':['SPARC','THBS2','FN1','POSTN'],
'Atlas_myCAF4':['TAGLN','ACTA2','MYLK','MYL9'],
'Atlas_iCAF5':['CXCL12','CXCL14','IL6','CXCL1','CXCL2'],
'Atlas_pericyte5':['RGS5','ABCC9','ENPEP','PLXDC1','PDGFRB'],
'FNDC1_CAF6':['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3'],
'FNDC1_nonstructural3':['FNDC1','MMP2','P4HA3'],
'CD36_CAF5':['CD36','FABP4','APOE','STEAP4','APOLD1'],
'Mineralization_control3':['SPP1','IBSP','BGLAP']}
protocol={'date':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Signature projection onto frozen 8174 candidate cells; not reference label transfer or independent replication','panels':panels,'signatures':'Short published figure/code panels; not complete unbiased marker sets. Zheng mCAF6 follows Supplementary Fig S2I; VCAN/COMP higher in iCAF. Atlas from describing_markers.qmd finalized block. Residual processing2/FAP-MMP11/FNDC1-MMP2-P4HA3 are explicitly posthoc mechanistic contrasts, not new validated signatures.','normalization':'Per-cell mean log1p(CP10k) using original full-library denominator and source-specific canonical gene map; missing features never set to zero','biological_unit':'source-local donor, within donor/site mean target minus reference, equal-weight sites then donors','filters':'current doublet positives excluded; comparisons restrict both sides to Fibro; thresholds 1 and 5 cells per donor-site-arm','contrasts':['F2_C3 vs other F2 Fibro','F2_C3 vs all other Fibro','F2_C3 vs F2_C1 Fibro','F2_C3 vs F4_C8 Fibro'],'statistics':'Two-sided exact donor signflip, N>=2 and nonzero; raw P<.05 nominal; BH across all 15 panels x 4 contrasts within each source/matrix/min-cell threshold; all tests retained; N<10 exploratory; no cell-level P','correlation':'Descriptive Spearman only within source Fibro, including gene-disjoint comparisons and per-donor correlations; no causal or identity-equivalence inference','cross_source':'No pooling. Target absent or not Fibro is recorded; no pseudo-validation','selection':'Published panels and comparisons fixed before score output; no re-clustering or label changes'}
(O/'resources/PROTOCOL.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
m=load_current_cells().set_index('cell_id');assert len(m)==8174 and m.index.is_unique
m[m.subgroup_id.eq('F2_C3')].to_csv(O/'tables/F2C3_metadata.tsv',sep='\t')
m[m.subgroup_id.eq('F2_C3')].groupby(['dataset','donor_record','site','lineage_stratum','candidate_tier'],dropna=False).size().rename('cells').reset_index().to_csv(O/'tables/F2C3_composition.tsv',sep='\t',index=False)
mp=pd.read_csv(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t')
pres=sparse.load_npz(P/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz').tocsr()
sidx=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample').presence_row.to_dict()
genes=sorted(set(sum(panels.values(),[]))|{'MMP9','OLFML3','COL5A3','HGF','PRSS23','DCN','BGN','COL1A2'})
allstats=[];effects=[];summaries=[];gene_summ=[];cellframes=[];cov=[];corr=[];rosters=[]
def perm(d):
 d=np.asarray(d,float);n=len(d)
 if n<2 or not np.any(d!=0):return np.nan
 signs=np.array(list(itertools.product([-1.,1.],repeat=n)))
 return float((np.abs(signs@d/n)>=abs(d.mean())-1e-12).mean())
for mode,path in [('raw',P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad'),('corrected',P/'phase58_F1F5_round1_20261004/objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==m.index.tolist();X=a.X.tocsr();lib=np.asarray(X.sum(1)).ravel();feat=a.var_names
 for ds in sorted(m.dataset.unique()):
  ii=np.flatnonzero(m.dataset.eq(ds));meta=m.iloc[ii].copy();meta['library']=lib[ii];mapping=mp[mp.dataset.eq(ds)&mp.symbol.isin(genes)].copy();jj=feat.get_indexer(mapping.feature);assert (jj>=0).all()
  present=np.asarray(pres[[sidx[s] for s in meta['sample'].unique()]][:,jj].sum(0)).ravel()==meta['sample'].nunique();mapping=mapping.iloc[np.flatnonzero(present)];jj=jj[present];syms=mapping.symbol.tolist();counts=X[ii][:,jj].toarray();norm=np.log1p(np.divide(counts,lib[ii,None],out=np.zeros_like(counts,dtype=float),where=lib[ii,None]>0)*1e4)
  sx=pd.DataFrame(norm,index=meta.index,columns=syms);sc=pd.DataFrame(index=meta.index)
  for name,gs in panels.items():
   avail=[g for g in gs if g in syms];cov.append(dict(dataset=ds,matrix=mode,panel=name,expected=len(gs),measured=len(avail),missing=';'.join(set(gs)-set(avail))))
   sc[name]=sx[avail].mean(axis=1) if len(avail)==len(gs) else np.nan
  keep=~meta.doublet_positive_any.astype(bool);fib=meta.lineage_stratum.eq('Fibro')&keep
  cf=meta[['dataset','sample','donor_record','site','class_id','subgroup_id','lineage_stratum','doublet_positive_any','library']].join(sc);cf['matrix']=mode;cellframes.append(cf.reset_index())
  for st,group in meta[keep].groupby('subgroup_id'):
   for scope in ['all','Fibro']:
    ids=group.index if scope=='all' else group.index[group.lineage_stratum.eq('Fibro')]
    if not len(ids):continue
    donor_means=sc.loc[ids].groupby(meta.loc[ids,'donor_record']).mean()
    for name in panels:summaries.append(dict(dataset=ds,matrix=mode,subgroup_id=st,scope=scope,panel=name,cells=len(ids),donors=len(donor_means),cell_mean=sc.loc[ids,name].mean(),donor_equal_mean=donor_means[name].mean()))
    if ds=='SCPCP000017':
     for g in syms:
      pos=meta.index.get_indexer(ids);gene_summ.append(dict(dataset=ds,matrix=mode,subgroup_id=st,scope=scope,gene=g,cells=len(ids),donors=meta.loc[ids].donor_record.nunique(),detection=float((counts[pos,syms.index(g)]>=(1 if mode=='raw' else .5)).mean()),mean_log1p_cp10k=sx.loc[ids,g].mean(),donor_equal_mean=sx.loc[ids,g].groupby(meta.loc[ids,'donor_record']).mean().mean()))
  if fib.sum()>2:
   for left,right in [('F2C3_identity3','Zheng_mCAF6'),('F2C3_assembly5','Zheng_mCAF6'),('F2C3_processing2','Zheng_mCAF_nonstructural2'),('F2C3_assembly5','Atlas_mCAF4'),('F2C3_assembly5','FNDC1_CAF6')]:
    for donor,ids in [('ALL_CELLS_descriptive',meta.index[fib])]+[(donor,g.index) for donor,g in meta[fib].groupby('donor_record') if len(g)>=5]:
     v=sc.loc[ids,[left,right]].dropna();rho=spearmanr(v[left],v[right]).statistic if len(v)>2 and v[left].nunique()>1 and v[right].nunique()>1 else np.nan
     corr.append(dict(dataset=ds,matrix=mode,donor=donor,left=left,right=right,cells=len(v),rho=rho))
  for threshold in [1,5]:
   for comp in ['other_F2','other_Fibro','F2_C1','F4_C8']:
    patient_diffs={}
    for donor,dg in meta[fib].groupby('donor_record'):
     site_diffs=[]
     for site,sg in dg.groupby('site',dropna=False):
      target=sg.index[sg.subgroup_id.eq('F2_C3')]
      if comp=='other_F2':ref=sg.index[sg.class_id.eq('F2')&~sg.subgroup_id.eq('F2_C3')]
      elif comp=='other_Fibro':ref=sg.index[~sg.subgroup_id.eq('F2_C3')]
      else:ref=sg.index[sg.subgroup_id.eq(comp)]
      if min(len(target),len(ref))<threshold:continue
      site_diffs.append(sc.loc[target].mean()-sc.loc[ref].mean());rosters.append(dict(dataset=ds,matrix=mode,min_cells=threshold,comparison=comp,donor_record=donor,site=site,target_cells=len(target),reference_cells=len(ref)))
     if site_diffs:patient_diffs[donor]=pd.DataFrame(site_diffs).mean()
    dif=pd.DataFrame.from_dict(patient_diffs,orient='index',columns=list(panels))
    for name in panels:
     v=dif[name].dropna();n=len(v)
     allstats.append(dict(dataset=ds,matrix=mode,min_cells=threshold,comparison=comp,panel=name,N=n,mean_effect=v.mean() if n else np.nan,positive=int((v>0).sum()),negative=int((v<0).sum()),P=perm(v),LOPO_min=(min((v.sum()-x)/(n-1) for x in v) if n>1 else np.nan),family=f'{ds}|{mode}|min{threshold}|all_panels_and_contrasts'))
     for donor,e in v.items():effects.append(dict(dataset=ds,matrix=mode,min_cells=threshold,comparison=comp,panel=name,donor_record=donor,effect=e))
  print('SCORED',mode,ds,'cells',len(meta),'F2C3',int(meta.subgroup_id.eq('F2_C3').sum()),flush=True)
 del a,X
stats=pd.DataFrame(allstats);stats['FDR']=np.nan
for fam,g in stats.groupby('family'):
 ix=g.index[g.P.notna()]
 if len(ix):stats.loc[ix,'FDR']=multipletests(stats.loc[ix,'P'],method='fdr_bh')[1]
stats['nominal']=stats.P.lt(.05);stats['FDR_significant']=stats.FDR.lt(.05)
for name,df in [('signature_patient_statistics',stats),('patient_effects',pd.DataFrame(effects)),('signature_subgroup_summary',pd.DataFrame(summaries)),('diagnostic_gene_summary',pd.DataFrame(gene_summ)),('signature_coverage',pd.DataFrame(cov)),('descriptive_signature_correlations',pd.DataFrame(corr)),('patient_site_contrast_roster',pd.DataFrame(rosters))]:df.to_csv(O/'tables'/f'{name}.tsv',sep='\t',index=False)
pd.concat(cellframes).to_csv(O/'tables/cell_signature_scores.tsv.gz',sep='\t',index=False)
# Full-transcriptome prior contrast is reused without altering statistics.
prior=P/'phase62_F1F5_internal_20261005/tables/genes/SCPCP000017__F2_C3__corrected__all.tsv.gz'
g=pd.read_csv(prior,sep='\t');r=pd.read_csv(str(prior).replace('__corrected__all','__raw__all'),sep='\t');g=g.merge(r[['symbol','mean_delta']],on='symbol',suffixes=('','_raw'))
g['display_eligible']=(g.target_detection>=.1)&(g.detection_difference>=.1)&(g.mean_delta>=.5)&(g.mean_delta_raw>0);g['display_score']=g.mean_delta*g.detection_difference
g.sort_values('display_score',ascending=False).to_csv(O/'tables/F2C3_full_transcriptome_prior_contrast.tsv.gz',sep='\t',index=False)
g[g.display_eligible].sort_values('display_score',ascending=False).head(50).to_csv(O/'tables/F2C3_top50_descriptive_markers.tsv',sep='\t',index=False)
source_files=[P/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',prior,P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv']
manifest=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in source_files]
(O/'resources/INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
assert len(pd.concat(cellframes))==8174*2
assert stats.loc[stats.N.lt(2),'P'].isna().all()
assert len(m[m.subgroup_id.eq('F2_C3')])==38
summary={'status':'completed','at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cells':8174,'F2C3_cells':38,'signature_panels':len(panels),'statistic_rows':len(stats),'N_lt2_no_P':True,'no_labels_changed':True,'target_primary_Fibro_cells':int(((m.subgroup_id=='F2_C3')&(m.lineage_stratum=='Fibro')&(m.dataset=='SCPCP000017')).sum()),'inference':'Exploratory signature projection; not independent reference mapping'}
(O/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
