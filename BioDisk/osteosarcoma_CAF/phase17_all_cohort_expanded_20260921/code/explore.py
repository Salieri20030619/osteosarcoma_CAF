from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']: os.environ[k]='4'
import json,datetime,hashlib,traceback,itertools
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase17_all_cohort_expanded_20260921'
for d in ['code','tables','objects','figures','logs','delivery']: (O/d).mkdir(parents=True,exist_ok=True)
def writej(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2))
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def truth(s):return s.fillna(False).astype(str).str.lower().isin(['true','1','1.0'])
def status(s,**kw):writej(O/'STATUS.json',dict(status=s,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw));print(s,kw,flush=True)
protocol={'scope':'40 previous sample records plus 12 ScPCA records; descriptive exploration only',
 'molecular_expansion':'n_genes>=300; raw Fib>=2 OR (Fib>=1 AND ECM>=3); exclude strong immune/endo/epi/myogenic/OC evidence at >1 percent with >=3 markers (OC>=2); known reviewed doublets remain excluded',
 'broad_pool':'molecular expansion and not both CNV high, union unchanged strict set; no frozen FibLike scope, no mandatory corrected count gate, no whole-cluster CNV exclusion, no osteo/mural exclusion',
 'supported_layer':'broad pool AND at least one CNV low AND no CNV high AND previous reference QC passed',
 'risk_labels':'single low with opposing high, unassessed, reference failure, mural program, osteo program, ambient absent remain explicit',
 'myCAF':'same customer 6 genes; >=2 genes raw counts>=1 and >=80% measured; operational expression flag, not identity proof',
 'integration':'full source count matrices; shared measured Ensembl genes; logCP10k with full library denominator; HVG2000 PC30 Harmony study_technology theta2; Leiden 0.5; same model refit on strict and supported layers',
 'statistics':'sample/donor-record descriptive effects only; no cell-level disease P values; cross-study independence unresolved; thresholds fixed for this attempt before examining results',
 'source_preservation':True}
writej(O/'PROTOCOL.json',protocol)
customer=json.loads((P/'phase16_scpca_supplement_20260920/resources/customer_panels.json').read_text())
surfaces=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
old=read(P/'phase10_batch_gates_20260919/tables/all_cell_unified_evidence.tsv.gz').set_index('cell_id')
old['strict']=truth(old.unified_strict); old['doublet']=truth(old.independent_doublet_excluded)
manifest=read(P/'phase10_batch_gates_20260919/resources/source_manifest.tsv').set_index('sample')
clinical=read(P/'phase14_client_objectives_20260920/tables/01_clinical_sample_audit_updated.tsv').set_index('sample')
meta=read(P/'phase16_scpca_supplement_20260920/manifest.tsv').set_index('scpca_sample_id')
sources=[]
for sm,g in old.groupby('sample',sort=True):
 g=g.copy();g['clinical_group']=str(clinical.loc[sm,'updated_comparison_role']);g['site']=str(clinical.loc[sm,'site_verified'])
 sources.append((sm,g,Path(manifest.loc[sm,'source_file'])))
for sm,r in meta.iterrows():
 b=P/'phase16_scpca_supplement_20260920'
 raw=read(b/'tables'/f'{sm}_cell_evidence.tsv.gz').set_index('cell_id');ev=read(b/'tables'/f'{sm}_identity_evidence.tsv.gz').set_index('cell_id').loc[raw.index]
 g=pd.DataFrame(index=raw.index)
 for c in ['sample','dataset','patient_id','site','n_genes']:g[c]=raw[c]
 for c in raw:
  if c.startswith(('det_','pct_')):g['raw_'+c]=raw[c]
 g['infer_call']=ev.infer_call;g['cnv2_call']=ev.copykat_call;g['reference_QC_pass']=truth(ev.reference_QC_pass)
 g['ambient_checked']=truth(ev.ambient_checked);g['strict']=truth(ev.unified_strict);g['doublet']=False
 g['original_lineage']=raw.lineage_candidate;g['technology']='ScPCA_'+raw.suspension_type.astype(str)
 g['clinical_group']='B2' if 'lung' in str(r.tissue_location).lower() and r.primary_or_metastasis=='Metastasis' else 'Unknown'
 sources.append((sm,g,b/'objects'/f'{sm}_qc.h5ad'))
assert len(sources)==52 and sum(len(g) for _,g,_ in sources)==299483
status('BUILDING_ALL_COHORT_CANDIDATES',samples=52)
parts=[];summ=[];coverage=[];expr=[];surface=[];all_evidence=[];source_stats=[]
target_genes=sorted(set(sum(customer.values(),[])+surfaces+['PTPRC','PECAM1','EPCAM']))
for sm,g,src in sources:
 g=g.copy();g['sample_key']=g.dataset.astype(str)+'|'+sm;g['donor_record']=g.dataset.astype(str)+'|'+g.patient_id.astype(str)
 g['study_technology']=g.dataset.astype(str)+'|'+g.technology.astype(str)
 fib=(g.raw_det_Fib>=2)|((g.raw_det_Fib>=1)&(g.raw_det_ECM>=3))
 exclusion=np.zeros(len(g),bool)
 for name,n in [('T_NK',3),('B_Plasma',3),('Myeloid',3),('DC',3),('Endo',3),('Epi',3),('Myogenic',3),('OC_specific',2)]:
  exclusion|=((g['raw_det_'+name]>=n)&(g['raw_pct_'+name]>1)).to_numpy()
 low=g.infer_call.eq('low_signal')|g.cnv2_call.eq('diploid')
 high=g.infer_call.eq('elevated_signal')|g.cnv2_call.eq('aneuploid')
 both=g.infer_call.eq('elevated_signal')&g.cnv2_call.eq('aneuploid')
 g['CNV_state']=np.select([both,low&high,low&~high,high],['both_high','low_high_conflict','low_no_high','high_without_low'],default='no_low_evidence')
 g['mural_flag']=(g.raw_det_Mural>=3)&(g.raw_pct_Mural>.2)
 g['osteo_flag']=(g.raw_det_Osteo>=3)&(g.raw_pct_Osteo>.2)
 g['expanded']=((g.n_genes>=300)&fib&~exclusion&~g.doublet&~both)|g.strict
 g['supported']=g.expanded&low&~high&truth(g.reference_QC_pass)
 # single-high without normal support stays a separate risk-only stratum, outside clustering pool.
 g['expanded']=g.expanded&~(high&~low)
 assert (~g.strict|g.expanded).all()
 keep=g.expanded.to_numpy()
 summ.append({'sample':sm,'dataset':g.dataset.iloc[0],'QC_cells':len(g),'strict':int(g.strict.sum()),'supported':int(g.supported.sum()),'expanded':int(keep.sum()),'expanded_CNV_conflict':int((g.expanded&low&high).sum()),'expanded_no_low':int((g.expanded&~low).sum()),'mural_flag':int((g.expanded&g.mural_flag).sum())})
 cols=['sample','dataset','patient_id','site','original_lineage','strict','supported','expanded','CNV_state','reference_QC_pass','ambient_checked','mural_flag','osteo_flag','doublet']
 all_evidence.append(g[cols].reset_index(names='cell_id'))
 if not keep.any():print('SAMPLE',summ[-1],flush=True);continue
 a=ad.read_h5ad(src);assert a.obs_names.is_unique and set(a.obs_names)==set(g.index);a=a[g.index].copy()
 x=a.layers['counts'].tocsr();sy=a.var.symbol.astype(str).to_numpy();den=np.asarray(x.sum(1)).ravel()
 vals=np.full((len(g),len(target_genes)),np.nan)
 for j,gene in enumerate(target_genes):
  ix=np.flatnonzero(sy==gene)
  if len(ix): vals[:,j]=np.asarray(x[:,ix].sum(1)).ravel()
 for name,gs in customer.items():
  ix=[target_genes.index(z) for z in gs];measured=np.isfinite(vals[:,ix]).all(0);cov=float(measured.mean())
  g[name+'_positive']=(np.sum(vals[:,ix]>=1,axis=1)>=2)&(cov>=.8)
  g[name+'_score']=np.nanmean(np.log1p(vals[:,ix]/np.maximum(den,1)[:,None]*1e4),axis=1)
  coverage.append({'sample':sm,'panel':name,'measured':int(measured.sum()),'requested':len(gs),'coverage':cov})
 for layer in ['strict','supported','expanded']:
  mask=g[layer].to_numpy();summ[-1][layer+'_myCAF_raw']=int((mask&g.myCAF_customer_positive).sum())
  for label,sel in [('all',mask),('myCAF',mask&g.myCAF_customer_positive.to_numpy())]:
   for j,gene in enumerate(target_genes):
    expr.append({'sample':sm,'dataset':g.dataset.iloc[0],'donor_record':g.donor_record.iloc[0],'site':g.site.iloc[0],'clinical_group':g.clinical_group.iloc[0],'layer':layer,'set':label,'gene':gene,'cells':int(sel.sum()),'library_counts':float(den[sel].sum()),'gene_counts':float(vals[sel,j].sum()) if np.isfinite(vals[:,j]).all() else np.nan})
 # Surface denominators use all QC cells from each contributing sample, not just candidate cells.
 target=g.expanded&g.myCAF_customer_positive
 dump=(vals[:,[target_genes.index(z) for z in ['PTPRC','PECAM1','EPCAM']]]<1).all(1)
 for pair in [(z,) for z in surfaces]+list(itertools.combinations(surfaces,2)):
  ix=[target_genes.index(z) for z in pair]
  if not np.isfinite(vals[:,ix]).all():continue
  gate=dump&(vals[:,ix]>=1).all(1)
  surface.append({'sample':sm,'donor_record':g.donor_record.iloc[0],'site':g.site.iloc[0],'markers':'+'.join(pair),'QC_cells':len(g),'target_cells':int(target.sum()),'gate_cells':int(gate.sum()),'target_in_gate':int((gate&target).sum()),'mural_in_gate':int((gate&g.mural_flag).sum()),'CNV_high_in_gate':int((gate&high).sum())})
 ids=np.array([z.split('.')[0] for z in a.var_names])
 if g.dataset.iloc[0]=='GSE250015':
  gm=read(P/'geo_cnv_20260918/inputs'/sm/'genes.tsv').set_index('source_gene_id');ids=gm.ensembl_id.reindex(a.var_names).astype(str).to_numpy()
 valid=np.array([z.startswith('ENSG') for z in ids]);assert len(set(ids[valid]))==valid.sum()
 ob=g.loc[g.expanded,['sample','sample_key','dataset','patient_id','donor_record','site','clinical_group','study_technology','original_lineage','strict','supported','CNV_state','reference_QC_pass','ambient_checked','mural_flag','osteo_flag']+[z for z in g if z.endswith(('_positive','_score'))]].copy()
 ob['full_library_counts']=den[keep]
 b=ad.AnnData(x[keep][:,valid],obs=ob,var=pd.DataFrame({'symbol':sy[valid]},index=ids[valid]));parts.append(b)
 source_stats.append({'sample':sm,'path':str(src),'bytes':src.stat().st_size,'mtime_ns':src.stat().st_mtime_ns})
 print('SAMPLE',summ[-1],flush=True)
save(pd.DataFrame(summ).fillna(0),'sample_layer_counts.tsv');save(pd.concat(all_evidence),'all_cell_layer_evidence.tsv.gz');save(pd.DataFrame(coverage),'panel_coverage.tsv')
save(pd.DataFrame(expr),'sample_gene_expression.tsv.gz');save(pd.DataFrame(surface),'surface_screen_all_QC.tsv.gz')
writej(O/'source_stat_snapshot.json',source_stats)
common=sorted(set.intersection(*[set(b.var_names) for b in parts]));assert len(common)>1000
x=sparse.vstack([b[:,common].X for b in parts]).tocsr();ob=pd.concat([b.obs for b in parts]);assert ob.index.is_unique
var=parts[0].var.loc[common].copy();a=ad.AnnData(x.copy(),obs=ob,var=var);a.layers['counts']=x.copy()
a.X=x.astype(float).multiply((1e4/np.maximum(a.obs.full_library_counts.to_numpy(),1))[:,None]).tocsr();a.X.data=np.log1p(a.X.data)
def clean(b):
 for frame in [b.obs,b.var]:
  for c in frame:
   if frame[c].dtype==object or isinstance(frame[c].dtype,pd.StringDtype):frame[c]=frame[c].fillna('').astype(str).astype(object)
  frame.index=pd.Index(frame.index.to_numpy(dtype=object),dtype=object)
clean(a)
status('JOINT_MODELS_RUNNING',expanded_cells=a.n_obs,shared_genes=a.n_vars)
def model(b,label):
 w=b[:,~b.var.symbol.str.startswith(('MT-','RPL','RPS')).to_numpy()].copy();sc.pp.filter_genes(w,min_cells=3)
 sc.pp.highly_variable_genes(w,n_top_genes=min(2000,w.n_vars),flavor='seurat',batch_key='study_technology')
 w=w[:,w.var.highly_variable].copy();sc.pp.scale(w,max_value=10);sc.tl.pca(w,n_comps=min(30,w.n_obs-1,w.n_vars-1),svd_solver='arpack',random_state=20260921)
 z=w.obsm['X_pca'];h=hm.run_harmony(z,b.obs,'study_technology',theta=2,random_state=20260921,ncores=4,max_iter_harmony=30,verbose=False)
 zz=np.asarray(h.Z_corr);zz=zz.T if zz.shape!=z.shape else zz;assert zz.shape==z.shape and np.isfinite(zz).all()
 b.obsm['X_pca']=z;b.obsm['X_harmony']=zz
 sc.pp.neighbors(b,use_rep='X_harmony',n_neighbors=min(15,b.n_obs-1),random_state=20260921)
 sc.tl.leiden(b,resolution=.5,key_added='cluster',flavor='igraph',n_iterations=2,directed=False,random_state=20260921)
 sc.tl.umap(b,random_state=20260921)
 clean(b);b.write_h5ad(O/'objects'/f'{label}_integrated.h5ad',compression='gzip')
 b.obs.to_csv(O/'tables'/f'{label}_cell_labels.tsv.gz',sep='\t',index_label='cell_id')
 print('MODEL_COMPLETE',label,b.n_obs,b.obs.cluster.nunique(),flush=True)
 return b
a=model(a,'expanded');comparisons=[]
for layer in ['strict','supported']:
 b=model(a[a.obs[layer]].copy(),layer);shared=a.obs.loc[b.obs_names,'cluster'].astype(str)
 comparisons.append({'layer':layer,'cells':b.n_obs,'clusters':b.obs.cluster.nunique(),'ARI_vs_expanded_on_shared_cells':adjusted_rand_score(shared,b.obs.cluster)})
save(pd.DataFrame(comparisons),'layer_refit_comparison.tsv')
support=[];profile=[];directions=[]
scores=[c for c in a.obs if c.endswith('_score')]
for cl,g in a.obs.groupby('cluster',observed=True):
 dn=g.groupby('donor_record',observed=True).size()
 support.append({'cluster':cl,'cells':len(g),'samples':g.sample_key.nunique(),'datasets':g.dataset.nunique(),'donor_records_ge10':int((dn>=10).sum()),'donor_records_ge20':int((dn>=20).sum()),'max_donor_record_fraction':float(dn.max()/len(g)),'strict_cells':int(g.strict.sum()),'supported_cells':int(g.supported.sum()),'raw_myCAF_cells':int(g.myCAF_customer_positive.sum()),'mural_fraction':float(g.mural_flag.mean()),'osteo_fraction':float(g.osteo_flag.mean()),'CNV_conflict_fraction':float(g.CNV_state.eq('low_high_conflict').mean())})
 for donor,gg in g.groupby('donor_record',observed=True):
  for s in scores:profile.append({'cluster':cl,'donor_record':donor,'cells':len(gg),'program':s,'mean_logCP10k':float(gg[s].mean())})
 for layer in ['expanded','supported','strict']:
  base=a.obs if layer=='expanded' else a.obs[a.obs[layer]]
  for donor,gg in base.groupby('donor_record',observed=True):
   inside=gg[gg.cluster==cl];outside=gg[gg.cluster!=cl]
   if len(inside)<10 or len(outside)<10:continue
   for s in scores:directions.append({'cluster':cl,'layer':layer,'donor_record':donor,'inside_cells':len(inside),'outside_cells':len(outside),'program':s,'delta_mean_logCP10k':float(inside[s].mean()-outside[s].mean())})
save(pd.DataFrame(support),'expanded_cluster_support.tsv');save(pd.DataFrame(profile),'cluster_donor_programs.tsv');save(pd.DataFrame(directions),'within_donor_layer_effects.tsv')
# Describe gene ranks without cell-level inferential tests.
markers=[]
for cl in a.obs.cluster.unique():
 m=a.obs.cluster.eq(cl).to_numpy();diff=np.asarray(a.X[m].mean(0)-a.X[~m].mean(0)).ravel()
 for ix in np.argsort(diff)[-30:][::-1]:markers.append({'cluster':cl,'gene_id':a.var_names[ix],'symbol':a.var.symbol.iloc[ix],'delta_mean_logCP10k':float(diff[ix])})
save(pd.DataFrame(markers),'descriptive_top30_genes.tsv')
eff=pd.DataFrame(directions);agree=[]
for (cl,pr),g in eff.groupby(['cluster','program']):
 piv=g.pivot(index='donor_record',columns='layer',values='delta_mean_logCP10k')
 for layer in ['strict','supported']:
  if layer not in piv or 'expanded' not in piv:continue
  q=piv[['expanded',layer]].dropna()
  agree.append({'cluster':cl,'program':pr,'comparison':layer,'evaluable_donor_records':len(q),'same_direction_records':int((np.sign(q.expanded)==np.sign(q[layer])).sum()),'expanded_median_effect':q.expanded.median(),'comparison_median_effect':q[layer].median()})
save(pd.DataFrame(agree),'strict_expanded_direction_agreement.tsv')
# Paired site descriptions aggregate repeated samples from a donor within each source.
ex=pd.DataFrame(expr);ex['site_group']=np.where(ex.site.str.contains('lung',case=False,na=False),'lung','other')
agg=ex.groupby(['dataset','donor_record','site_group','layer','set','gene'])[['cells','library_counts','gene_counts']].sum(min_count=1).reset_index()
agg['CP10k']=agg.gene_counts/agg.library_counts.replace(0,np.nan)*1e4;save(agg,'donor_site_gene_expression.tsv.gz')
# 2D display and donor-balanced program heatmap.
fig,axs=plt.subplots(1,3,figsize=(17,5))
xy=a.obsm['X_umap']
for ax,c in zip(axs,['cluster','dataset','CNV_state']):
 for name in sorted(a.obs[c].astype(str).unique()):
  m=a.obs[c].astype(str).eq(name);ax.scatter(xy[m,0],xy[m,1],s=1,alpha=.6,label=name,rasterized=True)
 ax.set_title(c);ax.set_xticks([]);ax.set_yticks([]);ax.legend(markerscale=5,fontsize=6,loc='best')
fig.tight_layout();fig.savefig(O/'figures/expanded_joint_umap.png',dpi=180);plt.close(fig)
for src in source_stats:
 st=Path(src['path']).stat();assert st.st_size==src['bytes'] and st.st_mtime_ns==src['mtime_ns']
qa={'QC_sample_records':52,'QC_cells':299483,'strict_cells':int(a.obs.strict.sum()),'supported_cells':int(a.obs.supported.sum()),'expanded_cells':a.n_obs,'expanded_samples':a.obs.sample_key.nunique(),'shared_genes':a.n_vars,'clusters':a.obs.cluster.nunique(),'expanded_raw_myCAF':int(a.obs.myCAF_customer_positive.sum()),'strict_raw_myCAF':int((a.obs.strict&a.obs.myCAF_customer_positive).sum()),'CNV_conflict_cells':int(a.obs.CNV_state.eq('low_high_conflict').sum()),'CNV_no_low_cells':int(a.obs.CNV_state.eq('no_low_evidence').sum()),'mural_flag_cells':int(a.obs.mural_flag.sum()),'source_files_unchanged':True,'independent_patient_N_established':False,'inferential_tests':False}
assert qa['strict_cells']==4099 and a.n_obs==sum(r['expanded'] for r in summ)
writej(O/'delivery/QA.json',{k:int(v) if isinstance(v,np.integer) else v for k,v in qa.items()})
status('ANALYSIS_COMPLETE',**{k:int(v) if isinstance(v,np.integer) else v for k,v in qa.items()})
(O/'ANALYSIS_COMPLETE').write_text('ALL_COHORT_EXPANDED_EXPLORATION_COMPLETE\n')
