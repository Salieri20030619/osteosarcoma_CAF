from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,json,itertools,time,hashlib,traceback
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,confusion_matrix
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase77_C3_external_function_20261010'
for k in ['resources','tables','code','logs']:(O/k).mkdir(exist_ok=True)
PANELS={'published_COL5_4':['COL5A1','COL16A1','EMILIN1','ADAMTS2'],'collagen_core5':['COL3A1','COL5A1','COL5A2','PCOLCE','LOXL2'],'processing2':['PCOLCE','LOXL2']};endpoints=set(sum(PANELS.values(),[]))
protocol={'created_before_model_fit':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'task':'C3 functional generalization with outcome genes held out from mapping','reference':'SCPCP000017 current Fibro, exclude model positive doublets; raw counts; C3=positive; other Fibro=negative','external':'GSE237070 six native-called matrices: 2 primary patients and 4 adjacent patients, analyzed separately','secondary':'Other sources in frozen8174 are source-held-out mapping only, not a wholly new independent discovery cohort','feature_exclusion':'All Matrisome gene sets, all three endpoint panels, all lineage panel genes, mitochondrial/ribosomal and listed sex genes excluded; gene availability intersection fixed across evaluated sources and external query','feature_selection':'Top2000 positive-variance genes by equal-donor variance in training cells only; redo within every leave-one-donor-out fold; no outcome-based tuning','model':'L2 logistic C=1, liblinear, max_iter3000; equal aggregate class weights, and equal donor weights within each class; donor-weighted centering/scaling with SD floor .1; call at >=.5; score not calibrated probability','internal':'Leave each of 18 reference Fibro donors out, repeat full feature selection and model; report donor AUC/recall and full heldout confusion matrix descriptively','mapping_quality_rule':'median AUC among donors with both labels >=.75 AND overall heldout recall >=.5 AND C3 recovered in >=5/7 positive donors; failure limits interpretation rather than changing cutoffs','outcomes':PANELS,'residual':'Fit processing2 ~ published_COL5_4 + log1p(full library) on reference nonC3 Fibro with donors equally weighted; freeze coefficients and apply to queries; compare within donor/site','binary_statistics':'C3-like vs other Fibro-like; both arms >=1 or5; site equal then patient equal; exact two-sided donor signflip N>=2 nonzero; 4 endpoint BH family per source/tissue/minimum; external query lineage eligibility uses original competitive panels after holding out all endpoint genes','continuous_sensitivity':'Within donor/site Spearman frozen logit versus processing residual, average sites within donor, exact patient signflip; separate BH across all source/tissue correlations; at least5 cells per site and nonconstant vectors','no_source_pooling':True,'no_master_labels_changed':True,'no_new_functional_or_CAFl_identity_claim':True,'reference_archive_access':'Figshare api and public file returned403 from both local and server; no author reference label mapping completed'}
(O/'resources/PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2));sys.path.insert(0,str(P/'phase66_doublet_repair_DSC1_20261005/code'));from current_metadata import load_current_cells
m=load_current_cells().set_index('cell_id');a=ad.read_h5ad(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad');assert a.obs_names.tolist()==m.index.tolist();mp=pd.read_csv(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t');presence=sparse.load_npz(P/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz');pi=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample').presence_row
lineage=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text());libs=json.loads((P/'phase62_F1F5_internal_20261005/resources/libraries_canonical.json').read_text());exclude=set(sum(libs['matrisome'].values(),[]))|endpoints|set(sum(lineage.values(),[]))|{'XIST','TSIX','RPS4Y1','RPS4Y2','DDX3Y','UTY','EIF1AY','KDM5D','ZFY','USP9Y','TMSB4Y'}
queries=sorted((P/'phase37_requested_source_search_20260927/objects/GSE237070_reconstructed').glob('*/full_native_called_counts.h5ad'));assert len(queries)==6
mask=m.lineage_stratum.eq('Fibro')&~m.doublet_positive_any;datasets=sorted(m.loc[mask,'dataset'].unique());maps={};common=None;coverage=[]
for ds in datasets:
 dm=m[mask&m.dataset.eq(ds)];pa=np.asarray(presence[pi.loc[dm['sample'].unique()].to_numpy()].sum(0)).ravel()==dm['sample'].nunique();mpp=mp[mp.dataset==ds].drop_duplicates('symbol').set_index('symbol').feature
 valid=mpp[(a.var_names.get_indexer(mpp)>=0)&pa[a.var_names.get_indexer(mpp)]];maps[ds]=valid
 common=set(valid.index) if common is None else common&set(valid.index)
 coverage.append({'dataset':ds,'fibro_cells':len(dm),'donors':dm.donor_record.nunique(),'eligible_measured_symbols':len(valid),'endpoint_coverage':len(endpoints&set(valid.index))})
query_maps={}
for path in queries:
 q=ad.read_h5ad(path,backed='r');v=q.var.copy();valid=v[v.symbol.notna()&~v.symbol.duplicated(keep=False)].symbol;query_maps[path.parent.name]=pd.Series(valid.index,index=valid.values);common&=set(valid.values);q.file.close()
genes=sorted(g for g in common if g not in exclude and not g.startswith(('MT-','RPL','RPS')))
assert len(genes)>=2000;pd.DataFrame({'gene':genes}).to_csv(O/'resources/mapping_eligible_genes.tsv',sep='\t',index=False);pd.DataFrame({'gene':sorted(exclude)}).to_csv(O/'resources/mapping_excluded_genes.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/source_coverage.tsv',sep='\t',index=False)
print('COMMON_FEATURES',len(genes),flush=True)
ix=np.flatnonzero(mask&m.dataset.eq('SCPCP000017'));ref=m.iloc[ix].copy();lib=np.asarray(a.X[ix].sum(1)).ravel();ji=a.var_names.get_indexer(maps['SCPCP000017'].loc[genes]);assert (ji>=0).all();x=a.X[ix][:,ji].toarray().astype(float);x=np.log1p(x/lib[:,None]*1e4);y=ref.subgroup_id.eq('F2_C3').to_numpy();donor=ref.donor_record.to_numpy();assert y.sum()==34

def train(xx,yy,dd):
 dw=np.array([1/sum(dd==d) for d in dd]);dw/=dw.sum();mu=dw@xx;var=dw@(xx-mu)**2;sel=np.argsort(-var,kind='stable')[:2000];sd=np.maximum(np.sqrt(var[sel]),.1);zz=(xx[:,sel]-mu[sel])/sd
 w=np.zeros(len(yy))
 for c in [False,True]:
  ds=np.unique(dd[yy==c])
  for d in ds:w[(dd==d)&(yy==c)]=len(yy)/(2*len(ds)*sum((dd==d)&(yy==c)))
 model=LogisticRegression(C=1,solver='liblinear',max_iter=3000,random_state=20261010).fit(zz,yy,sample_weight=w);assert model.n_iter_[0]<3000
 return sel,mu[sel],sd,model

def predict(xx,t):
 sel,mu,sd,model=t;return model.decision_function((xx[:,sel]-mu)/sd)
folds=[];held=np.zeros(len(ref));recovered=0
for d in sorted(set(donor)):
 ho=donor==d;t=train(x[~ho],y[~ho],donor[~ho]);v=predict(x[ho],t);held[ho]=v;p=v>=0;yy=y[ho]
 auc=roc_auc_score(yy,v) if len(set(yy))==2 else np.nan
 recall=float(p[yy].mean()) if yy.any() else np.nan
 if yy.any() and p[yy].any():recovered+=1
 folds.append({'donor':d,'cells':int(ho.sum()),'C3_cells':int(yy.sum()),'predicted_C3like':int(p.sum()),'AUROC':auc,'recall':recall,'precision':float(yy[p].mean()) if p.any() else np.nan})
 print('FOLD',d,'C3',int(yy.sum()),'AUC',auc,'recall',recall,flush=True)
folds=pd.DataFrame(folds);folds.to_csv(O/'tables/LODO_mapping_performance.tsv',sep='\t',index=False);conf=confusion_matrix(y,held>=0,labels=[False,True]);quality={'median_donor_AUROC':float(folds.AUROC.median()),'pooled_C3_recall':float((held[y]>=0).mean()),'pooled_precision':float(y[held>=0].mean()) if (held>=0).any() else None,'pooled_AUROC_descriptive':float(roc_auc_score(y,held)),'positive_donors_recovered':recovered,'positive_donors':int(folds.C3_cells.gt(0).sum()),'confusion_TN_FP_FN_TP':conf.ravel().tolist()};quality['mapping_gate_passed']=quality['median_donor_AUROC']>=.75 and quality['pooled_C3_recall']>=.5 and recovered>=5
(O/'resources/MAPPING_QUALITY.json').write_text(json.dumps(quality,indent=2));ref[['donor_record','subgroup_id']].assign(heldout_logit=held,true_C3=y,predicted_C3like=held>=0).to_csv(O/'tables/reference_heldout_predictions.tsv.gz',sep='\t');final=train(x,y,donor);sel,mu,sd,clf=final;fg=[genes[i] for i in sel];pd.DataFrame({'gene':fg,'mean':mu,'sd':sd,'coefficient':clf.coef_[0]}).to_csv(O/'resources/frozen_mapping_model.tsv',sep='\t',index=False);(O/'resources/frozen_mapping_intercept.json').write_text(json.dumps({'intercept':float(clf.intercept_[0]),'threshold_logit':0.0}));del x
# Endpoint extraction keeps absent measurements as NA, never zero.
def scores(raw,den,genemap,varnames):
 out={}
 for name,gs in PANELS.items():
  if not set(gs)<=set(genemap.index):out[name]=np.full(raw.shape[0],np.nan);continue
  js=varnames.get_indexer(genemap.loc[gs]);assert (js>=0).all();out[name]=np.log1p(raw[:,js].toarray()/den[:,None]*1e4).mean(1)
 return out
rp=pd.DataFrame(scores(a.X[ix],lib,maps['SCPCP000017'],a.var_names),index=ref.index);nw=np.array([1/sum(donor[~y]==d) for d in donor[~y]]);Z=np.column_stack([np.ones((~y).sum()),rp.loc[~y,'published_COL5_4'],np.log1p(lib[~y])]);coef=np.linalg.lstsq(Z*np.sqrt(nw)[:,None],rp.loc[~y,'processing2']*np.sqrt(nw),rcond=None)[0];(O/'resources/frozen_residual_model.json').write_text(json.dumps({'terms':['intercept','published_COL5_4','log1p_library'],'coefficients':coef.tolist(),'training':'SCPCP000017 nonC3 Fibro, patient equal'},indent=2))
allcells=[]
def append(q,xx,den,gm,vnames,scope):
 idx=vnames.get_indexer(gm.loc[fg]);assert (idx>=0).all();v=np.log1p(xx[:,idx].toarray()/den[:,None]*1e4);logit=((v-mu)/sd)@clf.coef_[0]+clf.intercept_[0];q=q.copy();q['mapping_logit']=logit;q['C3like']=logit>=0;q['library']=den;q['validation_scope']=scope
 for name,vals in scores(xx,den,gm,vnames).items():q[name]=vals
 q['processing_residual']=q.processing2-coef[0]-coef[1]*q.published_COL5_4-coef[2]*np.log1p(den);allcells.append(q)
for ds in datasets:
 if ds=='SCPCP000017':continue
 ii=np.flatnonzero(mask&m.dataset.eq(ds));dm=m.iloc[ii][['dataset','sample','donor_record','site','subgroup_id']].copy();dm['tissue']='frozen_candidate_pool';dm['eligibility']='current_Fibro_doublet_excluded';append(dm,a.X[ii],np.asarray(a.X[ii].sum(1)).ravel(),maps[ds],a.var_names,'source_holdout_internal_pool')
lineageheld={k:[g for g in gs if g not in endpoints] for k,gs in lineage.items()};(O/'resources/query_eligibility_panels.json').write_text(json.dumps(lineageheld,indent=2));qsum=[]
for path in queries:
 q=ad.read_h5ad(path);raw=q.X.tocsr();den=np.asarray(raw.sum(1)).ravel();gm=query_maps[path.parent.name];allgs=sorted(set(sum(lineageheld.values(),[])));ls=[]
 for name,gs in lineageheld.items():
  have=[g for g in gs if g in gm.index];js=q.var_names.get_indexer(gm.loc[have]);ls.append(np.log1p(raw[:,js].toarray()/np.where(den>0,den,np.nan)[:,None]*1e4).mean(1))
 vv=np.column_stack(ls);mx=np.nanmax(vv,axis=1);label=np.array(list(lineageheld),object)[np.nan_to_num(vv,nan=-np.inf).argmax(1)];label[(mx<=0)|(np.sum(vv==mx[:,None],axis=1)!=1)]='Unresolved';keep=label=='Fibroblast';ob=q.obs.copy();ob['lineage_endpoint_held_out']=label;ob.to_csv(O/'tables'/f'{path.parent.name}_eligibility.tsv.gz',sep='\t');dm=pd.DataFrame(index=q.obs_names[keep]);dm['dataset']='GSE237070';dm['sample']=path.parent.name;dm['donor_record']=q.obs.loc[dm.index,'patient'];dm['site']=q.obs.loc[dm.index,'tissue'];dm['tissue']=q.obs.loc[dm.index,'tissue'];dm['subgroup_id']='external_unassigned';dm['eligibility']='endpoint_held_out_Fibroblast_like';append(dm,raw[keep],den[keep],gm,q.var_names,'external_new_cohort');qsum.append({'sample':path.parent.name,'patient':q.obs.patient.iloc[0],'tissue':q.obs.tissue.iloc[0],'all_cells':q.n_obs,'eligible_Fibro_like':int(keep.sum())});print('QUERY',path.parent.name,int(keep.sum()),flush=True)
cells=pd.concat(allcells);cells.index.name='cell_id';cells.reset_index().to_csv(O/'tables/query_cell_scores.tsv.gz',sep='\t',index=False);pd.DataFrame(qsum).to_csv(O/'tables/query_sample_inventory.tsv',sep='\t',index=False)

def test(v):
 v=np.asarray(v);n=len(v)
 if n<2 or not np.any(v):return np.nan
 assert n<=20
 return float(np.mean(np.abs(np.array(list(itertools.product([-1,1],repeat=n)))@v/n)>=abs(v.mean())-1e-12))
effects=[];stats=[];corrs=[];corrstats=[];roster=[]
for (ds,tissue),dg in cells.groupby(['dataset','tissue']):
 for minimum in [1,5]:
  for panel in list(PANELS)+['processing_residual']:
   vals=[]
   for donor,dd in dg.groupby('donor_record'):
    sitevals=[]
    for site,z in dd.groupby('site'):
     l=z[z.C3like];r=z[~z.C3like]
     if panel=='processing2':roster.append({'dataset':ds,'tissue':tissue,'min_cells':minimum,'donor':donor,'site':site,'target_cells':len(l),'reference_cells':len(r),'evaluable':min(len(l),len(r))>=minimum})
     if min(len(l),len(r))>=minimum and l[panel].notna().all() and r[panel].notna().all():sitevals.append(l[panel].mean()-r[panel].mean())
    if sitevals:
     v=float(np.mean(sitevals));vals.append(v);effects.append({'dataset':ds,'tissue':tissue,'min_cells':minimum,'panel':panel,'donor':donor,'effect':v})
   stats.append({'dataset':ds,'tissue':tissue,'min_cells':minimum,'panel':panel,'N':len(vals),'mean_effect':np.mean(vals) if vals else np.nan,'positive':sum(v>0 for v in vals),'P':test(vals),'family':f'{ds}|{tissue}|min{minimum}|four_functional_endpoints'})
 vals=[]
 for donor,dd in dg.groupby('donor_record'):
  vv=[]
  for site,z in dd.groupby('site'):
   if len(z)>=5 and z.mapping_logit.nunique()>1 and z.processing_residual.nunique()>1:vv.append(spearmanr(z.mapping_logit,z.processing_residual).statistic)
  if vv:
   v=np.mean(vv);vals.append(v);corrs.append({'dataset':ds,'tissue':tissue,'donor':donor,'rho':v})
 corrstats.append({'dataset':ds,'tissue':tissue,'N':len(vals),'mean_rho':np.mean(vals) if vals else np.nan,'positive':sum(v>0 for v in vals),'P':test(vals),'family':'all_sources_tissues_continuous_processing_sensitivity'})
for name,rows in [('functional_statistics',stats),('continuous_statistics',corrstats)]:
 d=pd.DataFrame(rows);d['FDR']=np.nan
 for fam,z in d[d.P.notna()].groupby('family'):d.loc[z.index,'FDR']=multipletests(z.P,method='fdr_bh')[1]
 d.to_csv(O/'tables'/f'{name}.tsv',sep='\t',index=False)
for name,rows in [('patient_effects',effects),('continuous_patient_correlations',corrs),('patient_site_roster',roster)]:pd.DataFrame(rows).to_csv(O/'tables'/f'{name}.tsv',sep='\t',index=False)
summary=cells.groupby(['dataset','tissue','donor_record']).agg(cells=('C3like','size'),C3like=('C3like','sum')).reset_index();summary['fraction']=summary.C3like/summary.cells;summary.to_csv(O/'tables/donor_mapping_counts.tsv',sep='\t',index=False)
result={'status':'completed','mapping_quality':quality,'eligible_genes':len(genes),'model_genes':len(fg),'query_cells':len(cells),'query_patient_groups':len(summary),'functional_records':len(stats),'continuous_records':len(corrstats),'external_dataset':'GSE237070','missing_author_reference':'Figshare public access403; published four-gene program comparison only'};(O/'COMPLETE.json').write_text(json.dumps(result,indent=2));print('COMPLETE',json.dumps(result),flush=True)
