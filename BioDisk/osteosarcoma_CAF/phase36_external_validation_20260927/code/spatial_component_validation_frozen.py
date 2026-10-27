from pathlib import Path
import json,datetime,re,hashlib,itertools
import pandas as pd,numpy as np,anndata as ad
from scipy.stats import spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase36_external_validation_20260927';E=P/'phase35_external_validation_prep_20260927'
protocol=json.loads((O/'resources/VALIDATION_PROTOCOL_FROZEN.json').read_text());markers=protocol['markers'];meta=pd.read_csv(O/'resources/sample_patient_manifest_corrected.tsv',sep='\t',keep_default_na=False).set_index('sample')
qa=json.loads((O/'delivery/IC_DOWNLOAD_QA.json').read_text());assert len(qa)==2 and all(x['status']=='complete' for x in qa)
emb=pd.read_csv(O/'downloads/GSE293065_ica.cell.embeddings.all.csv.gz',index_col=0);assert emb.index.is_unique and np.isfinite(emb.to_numpy()).all()
ann=pd.read_excel(E/'sources/GSE293065_supplement_dataset.xlsx',header=1);ann=ann[ann.Use.astype(str).str.lower().eq('true')].copy();assert len(ann)==144
groups={
 'fibroblast':ann.loc[ann.Type.eq('Stroma')&ann.broad.isin(['Fibroblast','Myofibroblast']),'IC'].tolist(),
 'endothelial':ann.loc[ann.Type.eq('Stroma')&ann.broad.eq('Endothelial'),'IC'].tolist(),
 'mural':ann.loc[ann.Type.eq('Stroma')&ann.broad.isin(['Smooth muscle','Pericyte']),'IC'].tolist(),
 'cancer':ann.loc[ann.Type.eq('Cancer'),'IC'].tolist(),
}
assert all(groups.values()) and all(set(v)<=set(emb.columns) for v in groups.values())
# Continuous signed component sums: no 2.5% high-pass or Cancer>=0.4 gate imported from the source paper.
scores=pd.DataFrame({name:emb[ids].sum(axis=1) for name,ids in groups.items()})
annotations=[]
for name,ics in groups.items():
 for ic in ics:
  rec=ann[ann.IC.eq(ic)].iloc[0];genes=set(re.split(r'[,;|\s]+',str(rec.contrib_genes)))
  annotations.append(dict(group=name,IC=ic,Type=rec.Type,broad=rec.broad,annotation=rec.Annotation,candidate_contributors=';'.join(g for g in markers if g in genes)))
pd.DataFrame(annotations).to_csv(O/'tables/spatial_component_groups_and_candidate_circularity.tsv',sep='\t',index=False)
rows=[];coverage=[];spots=[]
for sm,r in meta.iterrows():
 a=ad.read_h5ad(E/'objects/GSE293065'/f'{sm}_counts_spatial.h5ad');keys=pd.Index(sm+'_'+a.obs.source_barcode.astype(str));hit=keys.isin(scores.index)
 coverage.append(dict(sample=sm,input_spots=a.n_obs,IC_matched=int(hit.sum()),IC_unmatched=int((~hit).sum()),author_IC_spots_for_sample=int(emb.index.str.startswith(sm+'_').sum())))
 if not hit.any():continue
 assert keys[hit].is_unique
 lib=np.asarray(a.X.sum(axis=1)).ravel();symbols=a.var.gene_symbol.astype(str).to_numpy();z=scores.loc[keys[hit]].copy();z.index=a.obs_names[hit]
 for marker in markers:
  ids=np.flatnonzero(symbols==marker);v=np.asarray(a.X[:,ids].sum(axis=1)).ravel().astype(float)
  expr=np.log1p(np.divide(v,lib,out=np.full(len(v),np.nan),where=lib>0)*10000)[hit];z['expr_'+marker]=expr
  for group in groups:
   target=z[group].to_numpy();ok=np.isfinite(expr)&np.isfinite(target);rho=float(spearmanr(expr[ok],target[ok]).statistic) if ok.sum()>1 and np.ptp(expr[ok])>0 and np.ptp(target[ok])>0 else np.nan
   rows.append(dict(sample=sm,patient_id=r.patient_id,event_month=r.months_from_diagnosis,site_group=r.site_group,marker=marker,component_group=group,rho=rho,matched_spots=int(hit.sum()),spot_P_not_used=True))
 z['sample']=sm;z['patient_id']=r.patient_id;z['site_group']=r.site_group;z['x']=a.obsm['spatial'][hit,0];z['y']=a.obsm['spatial'][hit,1];spots.append(z.reset_index(names='spot_id'))
rows=pd.DataFrame(rows);pd.DataFrame(coverage).to_csv(O/'tables/spatial_IC_alignment_coverage.tsv',sep='\t',index=False);rows.to_csv(O/'tables/sample_marker_component_correlations.tsv',sep='\t',index=False);pd.concat(spots,ignore_index=True).to_csv(O/'tables/spot_candidate_component_evidence.tsv.gz',sep='\t',index=False)
events=rows.groupby(['patient_id','site_group','event_month','marker','component_group']).rho.mean().reset_index();v=events.groupby(['patient_id','marker','component_group']).rho.mean().reset_index();v.to_csv(O/'tables/patient_marker_component_correlations.tsv',sep='\t',index=False)
rng=np.random.default_rng(20260927)
def test(d):
 d=np.asarray(d,float);d=d[np.isfinite(d)];n=len(d)
 if n<2 or np.all(d==0):return np.nan,'not_tested_insufficient_or_zero'
 obs=abs(d.mean())
 if n<=15:
  matrix=np.array(list(itertools.product([-1.,1.],repeat=n)));return float(np.mean(np.abs(matrix@d/n)>=obs-1e-12)),'exact_sign_flip'
 hits=0
 for _ in range(16):hits+=int(np.count_nonzero(np.abs(rng.choice([-1.,1.],size=(4096,n))@d/n)>=obs-1e-12))
 return (hits+1)/65537,'Monte_Carlo_sign_flip_65536'
def bh(x):
 a=np.asarray(x,float);out=np.full(len(a),np.nan);ix=np.flatnonzero(np.isfinite(a));ix=ix[np.argsort(a[ix])]
 if len(ix):out[ix]=np.minimum(1,np.minimum.accumulate((a[ix]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1])
 return out
stats=[]
for (marker,group),t in v.groupby(['marker','component_group']):
 d=t.rho.dropna().to_numpy();p,method=test(d);stats.append(dict(marker=marker,component_group=group,n_patients=len(d),mean_rho=float(np.mean(d)),median_rho=float(np.median(d)),positive_patients=int(np.sum(d>0)),P=p,test=method,family='marker_component_'+group))
s=pd.DataFrame(stats);s['FDR']=s.groupby('family').P.transform(lambda x:bh(x));s['nominal_P_lt_0_05']=s.P.lt(.05);s.to_csv(O/'tables/marker_component_patient_level_tests.tsv',sep='\t',index=False)
# Predeclared specificity comparisons assess relative association, not a protein gating rule.
wide=v.pivot(index=['patient_id','marker'],columns='component_group',values='rho').reset_index();contrasts=[]
for marker,t in wide.groupby('marker'):
 for other in ['cancer','endothelial','mural']:
  d=(t.fibroblast-t[other]).dropna().to_numpy();p,method=test(d);contrasts.append(dict(marker=marker,contrast='fibroblast_minus_'+other,n_patients=len(d),mean_rho_difference=float(np.mean(d)),median_rho_difference=float(np.median(d)),positive_patients=int(np.sum(d>0)),P=p,test=method,family='specificity_fibroblast_minus_'+other))
c=pd.DataFrame(contrasts);c['FDR']=c.groupby('family').P.transform(lambda x:bh(x));c['nominal_P_lt_0_05']=c.P.lt(.05);c.to_csv(O/'tables/marker_component_specificity_contrasts.tsv',sep='\t',index=False)
qa={'at':datetime.datetime.now().astimezone().isoformat(),'source_IC_rows':len(emb),'source_retained_components':len(ann),'groups':groups,'samples_with_alignment':len(rows['sample'].unique()),'matched_spots':sum(x['IC_matched'] for x in coverage),'unmatched_input_spots':sum(x['IC_unmatched'] for x in coverage),'patient_association_tests':len(s),'specificity_contrasts':len(c),'IC_scores_are_signed_sums_not_cell_proportions':True,'no_author_highpass_or_cancer_gate_applied':True,'spot_level_P_not_used':True,'contributor_gene_overlap_disclosed':True,'CAF_identity_or_flow_protein_validation_claimed':False}
(O/'delivery/SPATIAL_COMPONENT_VALIDATION_QA.json').write_text(json.dumps(qa,indent=2));(O/'SPATIAL_COMPONENT_VALIDATION_FINISHED').write_text(qa['at']);print(json.dumps(qa));print(s.to_string(index=False));print(c[c.marker.eq('GPC3')].to_string(index=False))
