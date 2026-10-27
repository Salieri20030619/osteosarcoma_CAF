from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,itertools,hashlib,inspect
import pandas as pd,numpy as np,anndata as ad,scanpy as sc,harmonypy
from sklearn.metrics import adjusted_rand_score
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
threadpool_limits(4)
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');J=P/'phase9_joint_20260919';O=P/'phase9_method_audit_20260919'
for d in ['tables','resources','code','logs']:(O/d).mkdir(parents=True,exist_ok=True)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
a=ad.read_h5ad(J/'objects/joint_CAF_integrated.h5ad');ob=a.obs.copy();ob['joint_cluster']=ob.joint_cluster.astype(str);cl=sorted(ob.joint_cluster.unique(),key=int)
# Record versions and real implementation defaults, no reliance on R docs for Python defaults.
(O/'resources/software.json').write_text(json.dumps({'scanpy':sc.__version__,'harmonypy':harmonypy.__version__,'run_harmony_signature':str(inspect.signature(harmonypy.run_harmony))},indent=2))
for f in [P/'phase7_review_states_20260917/identity/integration_parameters.json',P/'geo_preprocess_20260918/review/lineage_rules.json']:
 if f.exists():(O/'resources'/f.name).write_text(f.read_text())
# Threshold effects on within-study, within-donor evaluability.
rows=[]
for cut in [5,10,20,30,50]:
 for ds in sorted(ob.dataset.unique()):
  if ds=='GSE250015':continue
  for c in cl:
   n=0;total=0
   for d,g in ob[(ob.dataset==ds)&ob.donor_known].groupby('donor_record',observed=True):
    k=int(g.joint_cluster.eq(c).sum());r=len(g)-k
    if k>=cut and r>=cut:n+=1;total+=k
   rows.append({'min_cells_each_side':cut,'dataset':ds,'cluster':c,'paired_donor_records':n,'target_cells_supported':total,'min_two_sided_signflip_P':2**(1-n) if n else np.nan})
save(pd.DataFrame(rows),'pseudobulk_threshold_audit.tsv')
s=pd.read_csv(J/'tables/joint_cluster_review.tsv',sep='\t');rows=[]
for cut in [.6,.7,.8,.9]:
 for _,r in s.iterrows():rows.append({'donor_fraction_cutoff':cut,'cluster':r.cluster,'max_known_donor_fraction':r.max_known_donor_fraction,'flagged':r.max_known_donor_fraction>cut})
save(pd.DataFrame(rows),'donor_dominance_threshold_audit.tsv')
# Recompute eligibility gates from stored evidence without changing labels or source files.
e=pd.read_csv(P/'geo_identity_20260919/tables/all_cell_identity_evidence.tsv.gz',sep='\t',low_memory=False)
(O/'resources/identity_columns.json').write_text(json.dumps(e.columns.tolist(),indent=2))
for col in list(e.columns):
 if col.startswith('raw_') and col[4:].startswith(('det_','pos_','pct_')):e[col[4:]]=e[col]
fib=e.lineage_phase3.eq('FibLike');base_label=e.caf_supported_candidate.astype(bool)
def comp(e,prefix='',frac=.2,det='det_',osteo=True):
 out=np.zeros(len(e),bool)
 for name,n,p in [('Mural',3,frac),('OC',2,frac/2),('T_NK',3,frac),('Myeloid',4,frac),('Endo',3,frac),('Epi',3,frac),('Myogenic',2,frac)]:
  pc='OC_specific' if name=='OC' else name
  out|=(e[prefix+det+name]>=n)&(e[prefix+'pct_'+pc]>p)
 if osteo:out|=(e[prefix+det+'Osteo']>=3)&(e[prefix+'pct_Osteo']>frac)
 return out

def mixed(e,prefix='',frac=.2):
 immune=sum(e[prefix+'pct_'+k] for k in ['T_NK','B_Plasma','Myeloid','DC']);mes=sum(e[prefix+'pct_'+k] for k in ['Fib','Osteo','Mural'])
 return ((immune>frac)&(mes>frac))|((e[prefix+'pct_Endo']>frac)&(mes>frac))|((e[prefix+'pct_Epi']>frac)&(mes>frac))

def gate(ng=500,cnv=.2,frac=.2,fibgenes=2):
 raw=fib&e.joint.eq('both_low')&(e.det_Fib>=fibgenes)&~mixed(e,frac=frac)&~e.CTA_context.astype(bool)&(e.det_CTA<2)&~comp(e,frac=frac)&(e.cluster_both_high_fraction.fillna(1)<cnv)&(e.cluster_any_high_fraction.fillna(1)<cnv)&(e.n_genes>=ng)&e.reference_QC_pass.astype(bool)
 cor=(e.corrected_det_Fib>=fibgenes)&~comp(e,'corrected_',frac=frac)&~mixed(e,'corrected_',frac=frac)&(e.corrected_det_CTA<2)&e.ambient_checked.astype(bool)
 return raw&cor
baseline=gate();assert (baseline==base_label).all(),int((baseline!=base_label).sum())
rows=[]
conditions=[('baseline','500|0.2|0.2|2',{})]+[('n_genes',v,{'ng':v}) for v in [300,800,1000]]+[('cluster_CNV_risk_fraction',v,{'cnv':v}) for v in [.1,.3,.5]]+[('competing_program_percent',v,{'frac':v}) for v in [.1,.5,1.0]]+[('Fib_genes',v,{'fibgenes':v}) for v in [1,3,4]]
for name,val,kw in conditions:
 mask=gate(**kw)
 for ds in ['ALL']+sorted(e.dataset.unique()):
  m=np.ones(len(e),bool) if ds=='ALL' else e.dataset.eq(ds).to_numpy();rows.append({'varied_rule':name,'value':val,'dataset':ds,'candidates':int((mask&m).sum()),'added_vs_primary':int((mask&~baseline&m).sum()),'removed_vs_primary':int((~mask&baseline&m).sum()),'scope':'diagnostic_only_no_identity_change'})
save(pd.DataFrame(rows),'identity_one_factor_threshold_audit.tsv')
# Characterize cells excluded by lineage rules even with low CNV, without promoting them to CAFs.
rows=[]
for ds,g in e[fib&e.joint.eq('both_low')].groupby('dataset'):
 rows.append({'dataset':ds,'FibLike_dual_low':len(g),'raw_osteo_conflict':int(g.raw_osteo_conflict.sum()),'mural_competing':int(((g.det_Mural>=3)&(g.pct_Mural>.2)).sum()),'myogenic_competing':int(((g.det_Myogenic>=2)&(g.pct_Myogenic>.2)).sum()),'mixed_lineage':int(mixed(g).sum())})
save(pd.DataFrame(rows),'excluded_program_scope.tsv')
# Program scores are affected by coverage and library complexity; quantify associations, no hypothesis test.
rows=[]
for c,g in ob.groupby('joint_cluster',observed=True):
 rows.append({'cluster':c,'cells':len(g),'median_n_genes':g.n_genes.median(),'median_full_library_counts':g.full_library_counts.median(),'n_genes_q10':g.n_genes.quantile(.1),'n_genes_q90':g.n_genes.quantile(.9)})
save(pd.DataFrame(rows),'cluster_complexity_audit.tsv')
rows=[]
for key in [c for c in ob if c.startswith('program_')]:
 for ds,g in ob.groupby('dataset',observed=True):rows.append({'program':key,'dataset':ds,'rho_n_genes':spearmanr(g[key],g.n_genes).statistic,'rho_library':spearmanr(g[key],g.full_library_counts).statistic})
save(pd.DataFrame(rows),'program_complexity_correlations.tsv')
# Assess reference-program detectability on each full source, not arbitrary subtype assignment.
panels={'CD36':['CD36','FABP4','APOE','STEAP4','APOLD1'],'MYH11':['MYH11','RGS5','MCAM','MYH9','TGFBI'],'C7':['C7','PTGDS','TNXB','FGL2','IGF2','CXCL12'],'FNDC1':['FNDC1','COMP','MMP2','COL11A1','COL10A1','P4HA3']}
(O/'resources/reference_panels.json').write_text(json.dumps(panels,indent=2));rows=[];coverage=[]
manifest=pd.read_csv(J/'resources/source_manifest.tsv',sep='\t')
for r in manifest.to_dict('records'):
 b=ad.read_h5ad(r['path']);symbols=b.var.symbol.astype(str).to_numpy();x=b.layers['counts'];bo=ob.loc[b.obs_names];den=bo.full_library_counts.to_numpy()
 for name,genes in panels.items():
  present=[g for g in genes if g in symbols];coverage.append({'dataset':r['dataset'],'path':r['path'],'program':name,'measured':len(present),'requested':len(genes),'missing':','.join(g for g in genes if g not in symbols)})
  counts=np.column_stack([np.asarray(x[:,np.flatnonzero(symbols==g)].sum(1)).ravel() for g in present]);anchor=np.asarray(x[:,np.flatnonzero(symbols==name)].sum(1)).ravel() if name in symbols else None
  for cut in [1,2]:
   ng=(counts>=cut).sum(1)
   for (sm,c),inds in bo.groupby(['sample_key','joint_cluster'],observed=True).indices.items():
    rows.append({'dataset':r['dataset'],'sample_key':sm,'joint_cluster':c,'program':name,'count_threshold':cut,'cells':len(inds),'genes_measured':len(present),'anchor_detected_cells':int((anchor[inds]>=cut).sum()) if anchor is not None else np.nan,'ge2_program_genes_cells':int((ng[inds]>=2).sum()),'anchor_and_ge2_genes_cells':int(((anchor[inds]>=cut)&(ng[inds]>=2)).sum()) if anchor is not None else np.nan})
save(pd.DataFrame(rows),'reference_program_detectability.tsv');save(pd.DataFrame(coverage),'reference_program_coverage.tsv')
print('EVIDENCE_AND_THRESHOLDS_COMPLETE',flush=True)
# Limited diagnostic neighborhood/Leiden grid on frozen primary embedding; no main labels overwritten.
rows=[];overlap=[]
for k in [10,15,30]:
 b=ad.AnnData(np.zeros((a.n_obs,1)),obs=ob.copy());b.obsm['audit_embedding']=a.obsm['joint_main'].copy();sc.pp.neighbors(b,use_rep='audit_embedding',n_neighbors=k,random_state=20260919)
 for resolution in [.1,.2,.3,.5,.8]:
  for seed in [20260919,20260920]:
   name=f'k{k}_r{resolution}_s{seed}';sc.tl.leiden(b,resolution=resolution,key_added='audit',random_state=seed,flavor='igraph',n_iterations=2,directed=False);labs=b.obs.audit.astype(str)
   rows.append({'analysis':name,'neighbors':k,'resolution':resolution,'seed':seed,'clusters':labs.nunique(),'ARI_vs_main':adjusted_rand_score(ob.joint_cluster,labs)})
   for c in cl:
    u=ob.joint_cluster==c;best=max((int((u&(labs==v)).sum())/int((u|(labs==v)).sum()),v) for v in labs.unique());overlap.append({'analysis':name,'main_cluster':c,'best_cluster':best[1],'Jaccard':best[0]})
save(pd.DataFrame(rows),'neighborhood_resolution_seed_audit.tsv');save(pd.DataFrame(overlap),'neighborhood_cluster_overlap.tsv')
assert all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'] for r in manifest.to_dict('records'))
(O/'AUDIT_COMPLETE').write_text('READ_ONLY_ANALYSIS_AUDIT_COMPLETE\n');print('AUDIT_COMPLETE',flush=True);print(pd.DataFrame(rows).to_string(index=False))