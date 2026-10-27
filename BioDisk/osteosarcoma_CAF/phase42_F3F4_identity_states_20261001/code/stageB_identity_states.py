# SYNTAX_CHECK: python -m py_compile stageB_identity_states.py
# RUN_COMMAND: batch_bench/bin/python -u stageB_identity_states.py
# EXPECTED_SUCCESS_MARKER: STAGE_B_F3F4_IDENTITY_STATES_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,hashlib,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';O=P/'phase42_F3F4_identity_states_20261001';T=P/'phase30_two_layer_classification_20260926'
for sub in ['code','logs','tables','objects','figures','delivery','resources']:(O/sub).mkdir(parents=True,exist_ok=True)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
def writej(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
rc=1
try:
 assert (A/'setup.exit_code').read_text().strip()=='0' and (A/'COMPLETE').read_text().strip()=='SUCCESS'
 full=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=full.obs.copy();full.file.close();ob.index.name='cell_id'
 save(ob.reset_index(),'all_candidate_identity_registry.tsv.gz')
 joins=ob[ob.in_supported.astype(bool)].groupby(['expanded_class_id','supported_class_id'],observed=True).size().reset_index(name='cells');save(joins,'supported_expanded_class_correspondence.tsv')
 identity=[];profile=[];state_rows=[];modelqa=[];composition=[]
 fields=[c for c in ['copykat_call','infercnv_call','CNV_evidence','candidate_identity_status','provider_consensus_celltype_annotation','historical_broad_lineage','provisional_lineage','previous_lineage_review','known_doublet'] if c in ob]
 for layer in ['supported','expanded']:
  label='supported_class_id' if layer=='supported' else 'expanded_class_id';frame=ob[ob.in_supported.astype(bool)].copy() if layer=='supported' else ob.copy()
  for cls in ['F1','F2','F3','F4','F5']:
   g=frame[frame[label].eq(cls)]
   for field in fields:
    for value,z in g.groupby(field,observed=True,dropna=False):identity.append({'layer':layer,'class_id':cls,'field':field,'value':str(value),'cells':len(z),'donors':z.donor_record.nunique(),'datasets':z.dataset.nunique()})
   for c in [x for x in frame if x.startswith(('lineage_','program_'))]:
    if len(g):
     vals=pd.to_numeric(g[c]).groupby(g.donor_record).mean();profile.append({'layer':layer,'class_id':cls,'feature':c,'donors':len(vals),'donor_equal_mean':vals.mean(),'donor_median':vals.median(),'cells':len(g)})
  integrated=ad.read_h5ad(T/'objects'/f'{layer}_integrated.h5ad');assert set(integrated.obs_names)==set(frame.index)
  for cls in ['F3','F4']:
   g=frame[frame[label].eq(cls)]
   if not len(g):modelqa.append({'layer':layer,'class_id':cls,'cells':0,'state':'class_absent'});continue
   status('FITTING_INTERNAL_STATE',layer=layer,class_id=cls,cells=len(g))
   a=integrated[g.index].copy();a.obs=g.copy();a.layers['counts']=a.layers['counts'].tocsr()
   a.X=a.layers['counts'].astype(np.float32).multiply((1e4/pd.to_numeric(a.obs.full_library_counts).to_numpy())[:,None]).tocsr();a.X.data=np.log1p(a.X.data)
   sc.pp.highly_variable_genes(a,n_top_genes=2000,flavor='seurat',batch_key='study_technology')
   w=a[:,a.var.highly_variable].copy();sc.pp.scale(w,max_value=10);sc.tl.pca(w,n_comps=min(50,w.n_obs-1,w.n_vars-1),svd_solver='arpack',random_state=20261001)
   z=w.obsm['X_pca'];a.obsm['X_pca']=z
   if a.obs.study_technology.nunique()>1:
    result=hm.run_harmony(z,a.obs,'study_technology',random_state=20261001,ncores=4,verbose=False);h=np.asarray(result.Z_corr);h=h.T if h.shape!=z.shape else h;assert h.shape==z.shape and np.isfinite(h).all()
   else:h=z
   a.obsm['X_harmony']=h;sc.pp.neighbors(a,use_rep='X_harmony',n_neighbors=min(15,a.n_obs-1),random_state=20261001)
   for resolution,key in [(1,'internal_state'),(.5,'internal_state_res05')]:sc.tl.leiden(a,resolution=resolution,flavor='igraph',n_iterations=2,directed=False,key_added=key,random_state=20261001)
   sc.tl.umap(a,random_state=20261001)
   a.obs['state_id']=[f'{layer}_{cls}_C{x}' for x in a.obs.internal_state.astype(str)]
   for col in a.obs.select_dtypes(include='object'):a.obs[col]=a.obs[col].astype(str)
   a.obs['UMAP1']=a.obsm['X_umap'][:,0];a.obs['UMAP2']=a.obsm['X_umap'][:,1];a.obs.index.name='cell_id'
   sf=a.obs.reset_index();sf['layer']=layer;sf['parent_class']=cls;state_rows.append(sf)
   save(sf,f'{layer}_{cls}_internal_cell_states.tsv.gz')
   for state,sg in a.obs.groupby('state_id',observed=True):
    donors=sg.groupby('donor_record',observed=True).size();composition.append({'layer':layer,'class_id':cls,'state_id':state,'cells':len(sg),'donors':len(donors),'datasets':sg.dataset.nunique(),'largest_donor_fraction':donors.max()/len(sg),'largest_donor':donors.idxmax(),'supported_cells':int(sg.in_supported.astype(bool).sum())})
    for field in fields:
     for value,zg in sg.groupby(field,observed=True,dropna=False):identity.append({'layer':layer,'class_id':state,'field':field,'value':str(value),'cells':len(zg),'donors':zg.donor_record.nunique(),'datasets':zg.dataset.nunique()})
   fig,axes=plt.subplots(1,3,figsize=(14,4.2))
   for ax,col,title in zip(axes,['internal_state','dataset','candidate_identity_status'],['Internal states','Source','Identity evidence']):
    vals=a.obs[col].astype(str);levels=sorted(vals.unique());colors=plt.get_cmap('tab20')
    for k,level in enumerate(levels):
     mask=vals.eq(level).to_numpy();coords=a.obsm['X_umap'][mask];ax.scatter(coords[:,0],coords[:,1],s=5,alpha=.85,color=colors(k%20),label=level if len(level)<15 else f'{col[:3]}{k+1}')
    ax.set_title(title);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(fontsize=6,markerscale=2,loc='best');ax.spines[['top','right']].set_visible(False)
   fig.suptitle(f'{layer.title()} {cls}');fig.tight_layout();fig.savefig(O/'figures'/f'{layer}_{cls}_internal_UMAP.png',dpi=180);fig.savefig(O/'figures'/f'{layer}_{cls}_internal_UMAP.pdf');plt.close(fig)
   legend=[]
   for col in ['internal_state','dataset','candidate_identity_status']:
    for k,level in enumerate(sorted(a.obs[col].astype(str).unique())):legend.append({'panel':col,'plot_label':level if len(level)<15 else f'{col[:3]}{k+1}','full_label':level})
   save(pd.DataFrame(legend),f'{layer}_{cls}_figure_legend.tsv')
   a.write_h5ad(O/'objects'/f'{layer}_{cls}_internal_shared_embedding.h5ad',compression='gzip')
   modelqa.append({'layer':layer,'class_id':cls,'cells':a.n_obs,'embedding_genes':a.n_vars,'HVG':int(a.var.highly_variable.sum()),'states_res1':a.obs.internal_state.nunique(),'states_res05':a.obs.internal_state_res05.nunique(),'resolution_ARI':adjusted_rand_score(a.obs.internal_state,a.obs.internal_state_res05),'representation':'existing measured shared gene space; full-source expression retained separately for all hypothesis tests','gene_tests_performed':False,'state':'complete'})
   del a,w;gc.collect()
  del integrated;gc.collect()
 save(pd.DataFrame(identity),'class_and_state_identity_distribution.tsv');save(pd.DataFrame(profile),'class_patient_equal_lineage_program_profiles.tsv');save(pd.DataFrame(composition),'internal_state_patient_source_composition.tsv');save(pd.concat(state_rows,ignore_index=True),'all_F3F4_internal_cell_states.tsv.gz')
 mapping=O/'tables/all_F3F4_internal_cell_states.tsv.gz'
 protocol={'states_frozen_before_new_gene_P':True,'mapping_sha256':hashlib.sha256(mapping.read_bytes()).hexdigest(),'random_seed':20261001,'primary_resolution':1,'comparison_resolution':.5,'matches_prior_registered_graph_parameters':True,'HVG':2000,'PC_max':50,'neighbors_max':15,'Harmony_key':'study_technology','no_resolution_P_tuning':True,'full_gene_tests_use_stage41_full_counts_not_shared_embedding_features':True,'source_membership_or_CNV_calls_not_changed':True,'internal_states_operational_not_new_validated_subtypes':True}
 writej(O/'resources/INTERNAL_STATE_PROTOCOL_FROZEN.json',protocol);writej(O/'delivery/IDENTITY_STATE_QA.json',{'models':modelqa,'mapping_rows':sum(len(x) for x in state_rows),'frozen_candidate_set_unchanged':True,'CNV_results_reused':True,'CNV_new_runs':0,'class_correspondence_cells':int(joins.cells.sum())})
 (O/'COMPLETE').write_text('SUCCESS\n');status('complete',models=modelqa);print('STAGE_B_F3F4_IDENTITY_STATES_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();status('failed');(O/'FAILED').write_text(traceback.format_exc())
finally:(O/'identity_states.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)