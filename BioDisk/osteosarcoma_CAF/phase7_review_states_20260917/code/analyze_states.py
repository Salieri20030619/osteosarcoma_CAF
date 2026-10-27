"""SYNTAX_CHECK: python -m py_compile code/analyze_states.py
RUN_COMMAND: python -u code/analyze_states.py
EXPECTED_SUCCESS_MARKER: CAF_STATE_ANALYSIS_COMPLETE
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from itertools import product
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import KNeighborsClassifier
from threadpoolctl import threadpool_limits
threadpool_limits(limits=4);sc.settings.n_jobs=4
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=r/'phase7_review_states_20260917';out=p/'analysis';figs=p/'figures';out.mkdir(exist_ok=True);figs.mkdir(exist_ok=True)
SEED=20260917
PANELS={
'Fibroblast_core':['DCN','LUM','PDGFRA','DPT','COL3A1','COL6A1','COL14A1','PI16','TCF21'],
'Contractile':['ACTA2','TAGLN','MYL9','TPM2','CNN1'],
'Inflammatory':['IL6','CXCL12','CXCL14','CCL2','PTGS2','LIF','CXCL1','IL33'],
'Antigen_presentation':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','CIITA'],
'Matrix_remodeling':['COL1A1','COL1A2','COL3A1','COL5A1','COL6A3','FN1','POSTN','LOX','MMP2','CTHRC1','THBS2','COL11A1'],
'Osteogenic':['RUNX2','SP7','ALPL','IBSP','BGLAP','DMP1','PHEX','SPP1','SOST'],
'Mineral_calcium':['MGP','ENPP1','STC1','ANXA2','ANXA5','CASR','TRPM7','TNFRSF11B','TNFSF11'],
'Stress_core':['FOS','FOSB','JUN','JUNB','ATF3','NR4A1','HSPA1A','HSPA1B','DDIT3'],
'Cell_cycle':['MKI67','TOP2A','CDK1','UBE2C','CENPF','ASPM','AURKB','PCNA','MCM2','MCM5','TYMS','RRM2'],
'Mural_competing':['RGS5','CSPG4','NOTCH3','MCAM','MYH11','ABCC9','KCNJ8','COX4I2'],
'Immune_competing':['PTPRC','CD3D','CD3E','TYROBP','LST1','FCER1G','CD79A'],
'Myogenic_competing':['MYOD1','MYOG','MYF5','MYF6','MYH3','MYH8','ACTA1']}
PANELS['Lung_alveolar_reference']=['SPINT2','FGFR4','GPC3','NPNT','CES1']
PANELS['Lung_adventitial_reference']=['SFRP2','PI16','SERPINF1']
CORE=['Contractile','Inflammatory','Antigen_presentation','Matrix_remodeling','Osteogenic','Mineral_calcium']
PARAMS={'seed':SEED,'primary_dataset':'GSE152048','primary_cells':'reference-reviewed strict CAF candidates','main_Leiden_resolution':.5,'sensitivity_resolutions':[.3,.8],'PCA_dimensions_max':30,'neighbors':15,'HVG_requested':2000,'HVG_batch_key':'patient_id','patient_supported_state_min_cells_per_patient':10,'patient_supported_state_min_patients':3,'pseudobulk_state_comparison_min_cells_per_side':20,'gene_effect':'within-patient log2(CPM+1) state minus other CAF; median across paired patients','clinical_comparison':'descriptive donor-level effects; no bone-versus-lung inference with one lung donor','state_gene_tests':'exact patient-level paired sign-flip where at least 3 paired patients; BH separately per state-versus-other-CAF gene family; exploratory data-derived states','program_scores':'Scanpy expression-matched control-gene scores plus raw multi-gene detection; predefined candidate panels, not validated clinical signatures','cross_cohort':'separate projection/expression check, not pooled clinical replication'}
(out/'analysis_parameters.json').write_text(json.dumps(PARAMS,indent=2));(out/'gene_panels.json').write_text(json.dumps(PANELS,indent=2))
coverage=[]
def normalize(a):
 a.X=a.layers['counts'].astype(np.float32).copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
 return a

def genevec(a,g,counts=False):
 ix=np.flatnonzero(a.var.symbol.astype(str).values==g);x=a.layers['counts'] if counts else a.X
 return np.asarray(x[:,ix].sum(1)).ravel() if len(ix) else np.zeros(a.n_obs)

def score(a,label):
 pool=a.var_names[~a.var.symbol.astype(str).str.match(r'^(MT-|RPL|RPS)')].tolist()
 for name,genes in PANELS.items():
  ids=a.var_names[a.var.symbol.isin(genes)].tolist();present=[g for g in genes if (a.var.symbol==g).any()]
  coverage.append({'cohort_set':label,'program':name,'requested_genes':len(genes),'covered_genes':len(present),'missing':';'.join(sorted(set(genes)-set(present)))})
  if len(ids)>=2:sc.tl.score_genes(a,ids,score_name='program_'+name,gene_pool=pool,ctrl_size=50,random_state=SEED,use_raw=False)
  else:a.obs['program_'+name]=np.nan
  a.obs['detected_'+name]=sum((genevec(a,g,True)>0).astype(int) for g in genes)
  a.obs['meanlog_'+name]=sum(genevec(a,g) for g in genes)/max(len(present),1)
 return a

def embedding(a,key):
 sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=2000,batch_key='patient_id')
 technical=set(PANELS['Stress_core']+PANELS['Cell_cycle']+['XIST','RPS4Y1','DDX3Y','EIF1AY'])
 excluded=a.var.symbol.astype(str).str.match(r'^(MT-|RPL|RPS)')|a.var.symbol.isin(technical)
 a.var.loc[excluded,'highly_variable']=False
 b=a[:,a.var.highly_variable].copy();assert b.n_vars>=100
 sc.pp.scale(b,max_value=10);npc=min(30,b.n_obs-1,b.n_vars-1)
 sc.tl.pca(b,n_comps=npc,random_state=SEED,svd_solver='arpack');a.obsm['X_pca']=b.obsm['X_pca']
 sc.pp.neighbors(a,n_neighbors=min(15,a.n_obs-1),n_pcs=npc,random_state=SEED)
 for res in [.3,.5,.8]:sc.tl.leiden(a,resolution=res,key_added=key+'_'+str(res),random_state=SEED,flavor='igraph',n_iterations=2,directed=False)
 a.obs[key]=('C'+a.obs[key+'_0.5'].astype(str)).astype('category');sc.tl.umap(a,random_state=SEED)
 return b

def save_scores(a,label,cluster='caf_cluster'):
 cols=['sample','patient_id','site']+[c for c in [cluster] if c in a.obs]+['program_'+x for x in PANELS]+['meanlog_'+x for x in PANELS]+['detected_'+x for x in PANELS]
 cols += [c for c in ['projection_status','projection_mean_neighbor_distance','projection_distance_cutoff'] if c in a.obs]
 a.obs[cols].to_csv(out/(label+'_cell_programs.tsv.gz'),sep='\t',index_label='cell_id')
 num=['program_'+x for x in PANELS]+['meanlog_'+x for x in PANELS]
 a.obs.groupby(['sample','patient_id','site'],observed=True)[num].median().to_csv(out/(label+'_sample_program_medians.tsv'),sep='\t')

def project(query,train,b):
 ids=b.var_names;means=b.var['mean'].to_numpy();std=b.var['std'].to_numpy();std=np.where(std>0,std,1)
 loc=query.var_names.get_indexer(ids);covered=loc>=0
 if covered.mean()<.9:return None,float(covered.mean())
 mat=np.tile(means,(query.n_obs,1));mat[:,covered]=query.X[:,loc[covered]].toarray();mat=(mat-means)/std;mat=np.minimum(mat,10)
 pc=(mat-np.asarray(b.X.mean(axis=0)).ravel())@b.varm['PCs'];knn=KNeighborsClassifier(n_neighbors=15,weights='distance',n_jobs=4).fit(train.obsm['X_pca'],train.obs.caf_cluster.astype(str))
 train_dist=knn.kneighbors(train.obsm['X_pca'],n_neighbors=16)[0][:,1:].mean(axis=1)
 cutoff=float(np.quantile(train_dist,.99));query_dist=knn.kneighbors(pc)[0].mean(axis=1)
 query.obs['projection_mean_neighbor_distance']=query_dist;query.obs['projection_distance_cutoff']=cutoff
 query.obs['projection_status']=np.where(query_dist<=cutoff,'within_discovery_distance','out_of_distribution')
 prediction=knn.predict(pc).astype(object);prediction[query_dist>cutoff]='Unassigned_projection'
 return prediction,float(covered.mean())

assert (p/'identity/COMPLETE').exists()
a=normalize(ad.read_h5ad(p/'identity/objects/GSE152048_caf_strict.h5ad'));assert a.obs.caf_strict.all() and a.obs_names.is_unique
score(a,'GSE152048_strict');b=embedding(a,'caf_cluster');save_scores(a,'GSE152048_strict')
counts=a.layers['counts'];clusters=sorted(a.obs.caf_cluster.unique(),key=lambda x:int(x[1:]));patients=sorted(a.obs.patient_id.unique())
a.obs.to_csv(out/'strict_cell_annotations.tsv.gz',sep='\t',index_label='cell_id')
ct=pd.crosstab(a.obs.patient_id,a.obs.caf_cluster).reindex(columns=clusters,fill_value=0);ct.to_csv(out/'patient_by_state_cells.tsv',sep='\t');prop=ct.div(ct.sum(1),axis=0);prop.to_csv(out/'patient_by_state_fractions.tsv',sep='\t')
summary=[]
for c in clusters:
 ix=(a.obs.caf_cluster==c).values;pc=ct[c];shares=pc[pc>0]/pc.sum()
 z={'state':c,'cells':int(ix.sum()),'patients_any':int((pc>0).sum()),'patients_at_least_10_cells':int((pc>=10).sum()),'max_patient_fraction':float(shares.max()),'dominant_patient':shares.idxmax(),'patient_entropy':float(-(shares*np.log(shares)).sum()),'cross_patient_support':bool((pc>=10).sum()>=3),'median_genes':float(a.obs.loc[ix,'n_genes'].median()),'median_UMI':float(a.obs.loc[ix,'n_counts'].median())}
 for panel in PANELS:
  z[panel+'_score']=float(a.obs.loc[ix,'program_'+panel].mean());z[panel+'_pct_at_least_3_genes']=float((a.obs.loc[ix,'detected_'+panel]>=3).mean()*100)
 summary.append(z)
summary=pd.DataFrame(summary);summary.to_csv(out/'state_support_and_programs.tsv',sep='\t',index=False)
# Patient-balanced marker effect tables, without treating cells as independent replicates.
marker_tables=[];paired_effects=[]
for c in clusters:
 diffs=[];eligible=[]
 for patient in patients:
  within=(a.obs.patient_id==patient).values;target=within&(a.obs.caf_cluster==c).values;other=within&(a.obs.caf_cluster!=c).values
  if target.sum()<20 or other.sum()<20:continue
  t=np.asarray(counts[target].sum(0)).ravel();o=np.asarray(counts[other].sum(0)).ravel();d=np.log2(t/t.sum()*1e6+1)-np.log2(o/o.sum()*1e6+1)
  diffs.append(d);eligible.append(patient)
  paired_effects.append(pd.DataFrame({'state':c,'patient_id':patient,'ensembl_id':a.var_names,'log2CPM_difference':d}))
 ix=(a.obs.caf_cluster==c).values;ctarget=np.asarray(counts[ix].sum(0)).ravel();crest=np.asarray(counts[~ix].sum(0)).ravel()
 overall=np.log2(ctarget/ctarget.sum()*1e6+1)-np.log2(crest/crest.sum()*1e6+1)
 arr=np.asarray(diffs,dtype=np.float64);effect=np.median(arr,axis=0) if len(diffs) else overall
 m=pd.DataFrame({'state':c,'ensembl_id':a.var_names,'symbol':a.var.symbol.values,'median_paired_log2CPM_difference':effect,'eligible_paired_patients':len(eligible),'patients_positive':(arr>0).sum(0) if len(diffs) else np.zeros(a.n_vars,dtype=int),'positive_patient_fraction':(arr>0).mean(0) if len(diffs) else np.full(a.n_vars,np.nan),'target_pct_detected':np.asarray((counts[ix]>0).mean(0)).ravel()*100,'rest_pct_detected':np.asarray((counts[~ix]>0).mean(0)).ravel()*100,'effect_basis':'within_patient_median' if len(diffs) else 'pooled_cells_descriptive_only'})
 if len(diffs)>=3:
  assert len(diffs)<=12
  signs=np.array(list(product([-1.,1.],repeat=len(diffs))))
  observed=np.abs(arr.sum(axis=0));null=np.abs(signs@arr)
  pvals=(null>=observed[None,:]-1e-10).mean(axis=0)
  m['P_value']=pvals;m['FDR_BH']=multipletests(pvals,method='fdr_bh')[1]
  m['theoretical_min_two_sided_P_at_N']=2/(2**len(diffs));m['test_status']='exploratory_exact_paired_patient_sign_flip'
 else:
  m['P_value']=np.nan;m['FDR_BH']=np.nan;m['theoretical_min_two_sided_P_at_N']=np.nan;m['test_status']='not_tested_fewer_than_3_paired_patients'
 m['FDR_family']='state_'+c+'_vs_other_CAF_all_genes'
 marker_tables.append(m)
markers=pd.concat(marker_tables,ignore_index=True);markers.to_csv(out/'state_gene_effects_all.tsv.gz',sep='\t',index=False)
if paired_effects:pd.concat(paired_effects,ignore_index=True).to_csv(out/'patient_paired_state_gene_effects.tsv.gz',sep='\t',index=False)
top=markers[markers.target_pct_detected>=20].sort_values(['state','median_paired_log2CPM_difference'],ascending=[True,False]).groupby('state',sort=False).head(40);top.to_csv(out/'state_top40_effect_markers.tsv',sep='\t',index=False)
# Exact raw-count pseudobulk by independent patient.
pb=np.vstack([np.asarray(counts[(a.obs.patient_id==pat).values].sum(0)).ravel() for pat in patients]);assert np.all(pb==np.floor(pb))
pd.DataFrame(pb.T,index=a.var_names,columns=patients).to_csv(out/'patient_pseudobulk_counts.tsv.gz',sep='\t',index_label='ensembl_id')
meta=a.obs[['patient_id','sample','site','group','treatment_raw','clinical_context_note']].drop_duplicates().copy();meta['CAF_cells']=meta.patient_id.map(a.obs.patient_id.value_counts());meta=meta.merge(pd.read_csv(p/'review/clinical_enriched_20260917.tsv',sep='\t')[['patient_id','pathological_type_raw','tumor_size_cm_raw','necrosis_rate_raw','Ki67_percent_raw','disease_stage_raw']],on='patient_id',how='left',validate='one_to_one');meta.to_csv(out/'patient_metadata.tsv',sep='\t',index=False)
site={pat:('Lung' if (a.obs.loc[a.obs.patient_id==pat,'site']=='Lung').all() else 'Bone') for pat in patients};lung=np.array([site[pat]=='Lung' for pat in patients]);bone=~lung
program_effects=[];sample_program=a.obs.groupby('patient_id',observed=True)[['program_'+x for x in CORE]].median()
for panel in CORE:
 lv=sample_program.loc[[pat for pat in patients if site[pat]=='Lung'],'program_'+panel];bv=sample_program.loc[[pat for pat in patients if site[pat]=='Bone'],'program_'+panel]
 program_effects.append({'program':panel,'lung_patients':len(lv),'bone_patients':len(bv),'lung_median_patient_score':lv.median(),'bone_median_patient_score':bv.median(),'median_score_difference':lv.median()-bv.median(),'inference':'descriptive_only; no valid replicated lung group' if min(len(lv),len(bv))<2 else 'exploratory descriptive; no inferential test performed'})
pd.DataFrame(program_effects).to_csv(out/'bone_lung_program_effects.tsv',sep='\t',index=False)
if lung.any() and bone.any():
 logcpm=np.log2(pb/pb.sum(1)[:,None]*1e6+1);delta=np.median(logcpm[lung],axis=0)-np.median(logcpm[bone],axis=0)
 pd.DataFrame({'ensembl_id':a.var_names,'symbol':a.var.symbol.values,'lung_minus_bone_median_log2CPM':delta,'lung_patients':int(lung.sum()),'bone_patients':int(bone.sum()),'inference':'descriptive_only'}).to_csv(out/'bone_lung_gene_effects.tsv.gz',sep='\t',index=False)
# Candidate gene detection and mean expression by independent patient.
genes=sorted(set(sum(PANELS.values(),[])+['FAP','PDPN','THY1','CD34','ITGA11','PDGFRB','CD36','ENPP1','STC1','MGP']))
ge=[]
for pat in patients:
 ix=(a.obs.patient_id==pat).values
 for g in genes:
  ge.append({'patient_id':pat,'site':site[pat],'gene':g,'CAF_cells':int(ix.sum()),'pct_detected':float((genevec(a,g,True)[ix]>0).mean()*100),'mean_log1p':float(genevec(a,g)[ix].mean()),'gene_present':bool((a.var.symbol==g).any())})
pd.DataFrame(ge).to_csv(out/'patient_candidate_gene_expression.tsv',sep='\t',index=False)
# Quality sensitivity is applied to the same frozen candidate set.
quality=(a.obs.n_genes>=750)&~a.obs.high_complexity_flag
q=a.obs.loc[quality].groupby('patient_id',observed=True)[['program_'+x for x in CORE]].median();q.to_csv(out/'stricter_quality_patient_programs.tsv',sep='\t');a.obs['quality_sensitivity_keep']=quality
quality_counts=a.obs.groupby('patient_id',observed=True).quality_sensitivity_keep.agg(['sum','count']);quality_counts.to_csv(out/'quality_sensitivity_counts.tsv',sep='\t')
ari={'resolution_03_vs_05':adjusted_rand_score(a.obs['caf_cluster_0.3'],a.obs['caf_cluster_0.5']),'resolution_05_vs_08':adjusted_rand_score(a.obs['caf_cluster_0.5'],a.obs['caf_cluster_0.8'])}
# Freeze discovery before looking at the independent cohort.
(out/'DISCOVERY_FROZEN.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'patients':patients,'cells':a.n_obs,'states':clusters,'parameters':PARAMS,'marker_file_sha256':hashlib.sha256((out/'state_top40_effect_markers.tsv').read_bytes()).hexdigest()},indent=2))
ext=normalize(ad.read_h5ad(p/'identity/objects/GSE152048_caf_extended.h5ad'));score(ext,'GSE152048_extended');mapped,cov=project(ext,a,b);assert mapped is not None
ext.obs['caf_cluster']=mapped;shared=ext.obs_names.intersection(a.obs_names);ext.obs.loc[shared,'caf_cluster']=a.obs.loc[shared,'caf_cluster'].astype(str)
save_scores(ext,'GSE152048_extended');pd.crosstab(ext.obs.patient_id,ext.obs.caf_cluster).reindex(columns=clusters+['Unassigned_projection'],fill_value=0).to_csv(out/'extended_patient_by_state_cells.tsv',sep='\t')
strict_med=a.obs.groupby('patient_id',observed=True)[['program_'+x for x in CORE]].median();ext_med=ext.obs.groupby('patient_id',observed=True)[['program_'+x for x in CORE]].median();common=strict_med.index.intersection(ext_med.index)
sens=[]
for name in CORE:
 x=strict_med.loc[common,'program_'+name];y=ext_med.loc[common,'program_'+name];sens.append({'program':name,'paired_patients':len(common),'strict_extended_spearman':float(spearmanr(x,y).statistic),'median_extended_minus_strict_score':float((y-x).median())})
pd.DataFrame(sens).to_csv(out/'extended_set_sensitivity.tsv',sep='\t',index=False)
old=normalize(ad.read_h5ad(r/'phase3_identity_20260914/identity/objects/caf_strict.h5ad'));score(old,'GSE270231_strict');pred,old_cov=project(old,a,b)
old.obs['caf_cluster']=pred if pred is not None else 'not_projected';old.obs['patient_replication_eligible']=old.obs.patient_id!='UNCONFIRMED';save_scores(old,'GSE270231_strict')
pd.crosstab(old.obs['sample'],old.obs.caf_cluster).to_csv(out/'GSE270231_projected_state_counts.tsv',sep='\t')
# Conditional myCAF subdivision, only when expression and patient coverage support it.
contract=(a.obs.detected_Contractile>=3)&(a.obs.program_Contractile>0);my_counts=a.obs.loc[contract].groupby('patient_id',observed=True).size();my_info={'candidate_cells':int(contract.sum()),'patients_with_20_cells':int((my_counts>=20).sum()),'executed':False}
if contract.sum()>=150 and (my_counts>=20).sum()>=3:
 my=a[contract].copy();embedding(my,'mycaf_cluster');my.obs[['patient_id','sample','mycaf_cluster','program_Contractile','program_Osteogenic','program_Mineral_calcium']].to_csv(out/'mycaf_subcluster_cells.tsv.gz',sep='\t',index_label='cell_id');pd.crosstab(my.obs.patient_id,my.obs.mycaf_cluster).to_csv(out/'mycaf_subcluster_patient_counts.tsv',sep='\t');my.write_h5ad(out/'mycaf_subclusters.h5ad',compression='gzip');my_info['executed']=True
else:my_info['reason']='fewer than 150 contractile candidate cells or fewer than 3 patients with 20 cells; report continuous program instead'
(out/'mycaf_subdivision_status.json').write_text(json.dumps(my_info,indent=2))
pd.DataFrame(coverage).to_csv(out/'program_gene_coverage.tsv',sep='\t',index=False)
# Server-side figures; all associated numerical tables are preserved.
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
fig,axs=plt.subplots(1,3,figsize=(16,4.5))
for ax,key,title in zip(axs,['caf_cluster','patient_id','site'],['CAF state clusters','Independent patient','Sampling site']):sc.pl.umap(a,color=key,ax=ax,show=False,title=title,size=10)
fig.tight_layout();fig.savefig(figs/'CAF_states_patients_sites.png',dpi=180,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));prop.plot.bar(stacked=True,ax=ax,colormap='tab20');ax.set_ylabel('Fraction within strict CAF candidates');ax.set_xlabel('Patient');ax.legend(title='State',bbox_to_anchor=(1.02,1));fig.tight_layout();fig.savefig(figs/'Patient_state_composition.png',dpi=180,bbox_inches='tight');plt.close(fig)
vals=summary.set_index('state')[[x+'_score' for x in CORE]];z=(vals-vals.mean())/vals.std(ddof=0).replace(0,1)
fig,ax=plt.subplots(figsize=(9,max(3,len(clusters)*.5)));im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto');ax.set_xticks(range(len(CORE)),CORE,rotation=35,ha='right');ax.set_yticks(range(len(clusters)),clusters);fig.colorbar(im,ax=ax,label='Across-state z-score');fig.tight_layout();fig.savefig(figs/'CAF_state_programs.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(2,3,figsize=(13,7))
for ax,name in zip(axs.ravel(),CORE):
 x=sample_program['program_'+name].reindex(patients);ax.scatter(np.arange(len(patients)),x,c=['#B94731' if site[pat]=='Lung' else '#3067A8' for pat in patients]);ax.set_xticks(range(len(patients)),[pat.split(':')[-1] for pat in patients],rotation=45);ax.set_title(name);ax.set_ylabel('Median cell score per patient')
fig.tight_layout();fig.savefig(figs/'Patient_program_scores.png',dpi=180);plt.close(fig)
a.uns['phase7_analysis_parameters']=json.dumps(PARAMS);a.write_h5ad(out/'GSE152048_strict_CAF_states.h5ad',compression='gzip')
result={'state':'CAF_STATE_ANALYSIS_COMPLETE','completed_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_strict_cells':a.n_obs,'primary_patients':len(patients),'state_clusters':len(clusters),'cross_patient_supported_clusters':int(summary.cross_patient_support.sum()),'lung_patients_with_strict_CAF':int(lung.sum()),'bone_patients_with_strict_CAF':int(bone.sum()),'extended_cells':ext.n_obs,'external_GSE270231_cells':old.n_obs,'external_confirmed_patients':int(old.obs.loc[old.obs.patient_id!='UNCONFIRMED','patient_id'].nunique()),'external_unknown_mapping_cells':int((old.obs.patient_id=='UNCONFIRMED').sum()),'external_projection_HVG_coverage':old_cov,'external_projection_unassigned_cells':int((old.obs.caf_cluster=='Unassigned_projection').sum()),'projection_interpretation':'nearest-state exploratory projection with distance screen; not independent validation','resolution_ARI':ari,'mycaf_subdivision':my_info,'inferential_clinical_tests_performed':False}
(out/'analysis_summary.json').write_text(json.dumps(result,indent=2));(out/'COMPLETE').write_text('CAF_STATE_ANALYSIS_COMPLETE\n');print(json.dumps(result,indent=2),flush=True)
