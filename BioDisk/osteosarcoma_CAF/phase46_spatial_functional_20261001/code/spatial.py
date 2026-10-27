from common import *
import anndata as ad
from scipy.spatial import cKDTree
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
V=P/'phase35_external_validation_prep_20260927';D=P/'phase36_external_validation_20260927'
protocol=json.loads((O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json').read_text());programs=protocol['programs'];meta=pd.read_csv(D/'resources/sample_patient_manifest_corrected.tsv',sep='\t',keep_default_na=False).set_index('sample')
emb=pd.read_csv(D/'downloads/GSE293065_ica.cell.embeddings.all.csv.gz',index_col=0);ann=pd.read_excel(V/'sources/GSE293065_supplement_dataset.xlsx',header=1);ann=ann[ann.Use.astype(str).str.lower().eq('true')]
groups={'cancer':ann.loc[ann.Type.eq('Cancer'),'IC'].tolist(),'fibroblast':ann.loc[ann.Type.eq('Stroma')&ann.broad.isin(['Fibroblast','Myofibroblast']),'IC'].tolist(),'Tcell_IC':ann.loc[ann.Type.eq('Lymphoid')&ann.broad.eq('T cells'),'IC'].tolist()}
assert all(groups.values());dump(O/'resources/AUTHOR_COMPONENT_GROUPS.json',groups)
annotation=[]
for group,ids in groups.items():
 for ic in ids:
  r=ann[ann.IC.eq(ic)].iloc[0];contrib={canon(g.strip()) for g in str(r.contrib_genes).split(',') if canon(g.strip())}
  for name,genes in programs.items():annotation.append(dict(group=group,IC=ic,program=name,overlap=';'.join(sorted(set(genes)&contrib)),overlap_N=len(set(genes)&contrib)))
pd.DataFrame(annotation).to_csv(O/'tables/component_program_contributor_overlap.tsv',sep='\t',index=False)
comp=pd.DataFrame({g:emb[ids].sum(1) for g,ids in groups.items()});rows=[];coverage=[];spots=[];graphs=[]
for sm,r in meta.iterrows():
 dump(O/'SPATIAL_STATUS.json',{'state':'processing','sample':sm,'finished_samples':len(graphs)})
 a=ad.read_h5ad(V/'objects/GSE293065'/f'{sm}_counts_spatial.h5ad');x=a.X.tocsr();lib=np.asarray(x.sum(1)).ravel();assert np.isfinite(x.data).all() and (x.data>=0).all();gm=feature_map(a)
 keys=pd.Index(sm+'_'+a.obs.source_barcode.astype(str));match=keys.isin(comp.index);xy=np.asarray(a.obsm['spatial'],float);assert np.isfinite(xy).all();tree=cKDTree(xy);ds=tree.query(xy,k=2)[0][:,1];pitch=float(np.median(ds[ds>0]));pairs=tree.query_pairs(pitch*1.05,output_type='ndarray');ii=np.r_[pairs[:,0],pairs[:,1]];jj=np.r_[pairs[:,1],pairs[:,0]];graph=sparse.csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(a.n_obs,a.n_obs));degree=np.asarray(graph.sum(1)).ravel()
 scores={}
 for name,genes in {**programs,'CD8_RNA':['CD8A','CD8B'],'Tcell_RNA':['CD3D','CD3E','TRAC']}.items():
  have=[g for g in genes if g in gm];idx=[gm[g] for g in have];raw=x[:,idx].toarray().astype(float) if idx else np.empty((a.n_obs,0));log=np.log1p(raw/np.where(lib>0,lib,np.nan)[:,None]*1e4);scores[name]=np.mean(log,1) if have else np.full(a.n_obs,np.nan)
  coverage.append(dict(sample=sm,program=name,requested=len(genes),present=len(have),missing=';'.join(sorted(set(genes)-set(have))),full_library_zero_spots=int((lib==0).sum()),object_SHA256=sha(V/'objects/GSE293065'/f'{sm}_counts_spatial.h5ad') if name=='CD8_RNA' else ''))
 z=pd.DataFrame(scores,index=a.obs_names);z['sample']=sm;z['patient_id']=r.patient_id;z['site_group']=r.site_group;z['event_month']=r.months_from_diagnosis or 'primary_or_unknown';z['library_UMI']=lib;z['x']=xy[:,0];z['y']=xy[:,1];z['degree']=degree
 for group in groups:z[group]=np.nan;z.loc[match,group]=comp.loc[keys[match],group].to_numpy()
 ng=np.divide(graph@z.CD8_RNA.to_numpy(),degree,out=np.full(a.n_obs,np.nan),where=degree>0);z['neighbor_CD8_RNA']=ng
 controls=z[['cancer','fibroblast']].to_numpy();controls=np.column_stack([controls,np.log1p(lib),xy])
 for name in programs:
  s=z[name].to_numpy();values={'CD8_raw':corr(s,z.CD8_RNA),'CD8_adjusted':corr(s,z.CD8_RNA,controls),'neighbor_CD8_adjusted':corr(s,ng,controls),'Tcell_adjusted':corr(s,z.Tcell_RNA,controls),'fibroblast_component':corr(s,z.fibroblast),'cancer_component':corr(s,z.cancer)}
  for ep,v in values.items():rows.append(dict(sample=sm,patient_id=r.patient_id,site_group=r.site_group,event_month=r.months_from_diagnosis or 'primary_or_unknown',program=name,endpoint=ep,effect=v,spots=a.n_obs,IC_matched=int(match.sum()),no_spot_P=True))
 graphs.append(dict(sample=sm,spots=a.n_obs,matched_IC=int(match.sum()),unmatched_IC=int((~match).sum()),pitch=pitch,radius=pitch*1.05,undirected_edges=len(pairs),isolated_spots=int((degree==0).sum()),accepted_image=False,accepted_tumor_mask=False,CD8_cells_not_measured=True));spots.append(z.reset_index(names='spot_id'));np.savez_compressed(O/'objects'/f'{sm}_adjacency.npz',row=ii,col=jj,shape=graph.shape)
 print('PROCESSED',sm,a.n_obs,int(match.sum()),flush=True)
frame=pd.DataFrame(rows);frame.to_csv(O/'tables/sample_spatial_associations.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/program_gene_coverage.tsv',sep='\t',index=False);pd.DataFrame(graphs).to_csv(O/'tables/spatial_geometry_acceptance.tsv',sep='\t',index=False)
events=frame.groupby(['patient_id','site_group','event_month','program','endpoint'],dropna=False).effect.mean().reset_index();sites=events.groupby(['patient_id','site_group','program','endpoint']).effect.mean().reset_index();patients=sites.groupby(['patient_id','program','endpoint']).effect.mean().reset_index();patients.to_csv(O/'tables/patient_spatial_effects.tsv',sep='\t',index=False);sites.to_csv(O/'tables/patient_site_spatial_effects.tsv',sep='\t',index=False)
stats=test_table(patients,['program','endpoint'],'effect',lambda r:'spatial_'+r.endpoint);stats.to_csv(O/'tables/spatial_patient_statistics.tsv',sep='\t',index=False)
site_stats=test_table(sites,['site_group','program','endpoint'],'effect',lambda r:'spatial_site_'+r.site_group+'_'+r.endpoint);site_stats.to_csv(O/'tables/spatial_site_sensitivity_statistics.tsv',sep='\t',index=False)
allspots=pd.concat(spots,ignore_index=True);allspots.to_csv(O/'tables/all_spot_RNA_and_components.tsv.gz',sep='\t',index=False)
main=patients[patients.program.eq('NABA_CORE_MATRISOME')&patients.endpoint.eq('CD8_adjusted')].sort_values('patient_id');fig,ax=plt.subplots(figsize=(8,5));ax.axvline(0,color='gray',lw=1);ax.scatter(main.effect,np.arange(len(main)),color='#257b7f');ax.set_yticks(np.arange(len(main)),main.patient_id.str.replace('GSE293065_',''));ax.set_xlabel('Patient mean partial rank correlation');ax.set_title('Matrix RNA and CD8 RNA');fig.tight_layout();fig.savefig(O/'figures/patient_matrix_CD8_effects.png',dpi=180);fig.savefig(O/'figures/patient_matrix_CD8_effects.pdf');plt.close(fig)
wide=stats.pivot(index='program',columns='endpoint',values='mean_effect');fig,ax=plt.subplots(figsize=(10,8));im=ax.imshow(wide,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_yticks(np.arange(len(wide)),[z.replace('HALLMARK_','').replace('NABA_','').replace('_',' ') for z in wide.index],fontsize=8);ax.set_xticks(np.arange(len(wide.columns)),[z.replace('_',' ') for z in wide.columns],rotation=45,ha='right',fontsize=8);fig.colorbar(im,ax=ax,label='Patient mean rank correlation');fig.tight_layout();fig.savefig(O/'figures/spatial_functional_heatmap.png',dpi=180);fig.savefig(O/'figures/spatial_functional_heatmap.pdf');plt.close(fig)
with PdfPages(O/'figures/all_patient_spatial_RNA_maps.pdf') as pdf:
 for sm,z in allspots.groupby('sample',sort=True):
  fig,axs=plt.subplots(1,3,figsize=(12,4))
  for ax,name,label in zip(axs,['NABA_CORE_MATRISOME','CD8_RNA','cancer'],['Matrix RNA','CD8 RNA','Cancer signed IC']):
   v=ax.scatter(z.x,z.y,c=z[name],s=5,cmap='viridis',rasterized=True);fig.colorbar(v,ax=ax,shrink=.65);ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(label);ax.set_axis_off()
  fig.suptitle(sm+' '+z.patient_id.iloc[0]+'  RNA coordinates without histology');fig.tight_layout();pdf.savefig(fig)
  if sm in ['PT1','MR1']:fig.savefig(O/'figures'/f'{sm}_spatial_RNA.png',dpi=180)
  plt.close(fig)
qa={'at':datetime.datetime.now().astimezone().isoformat(),'samples':len(graphs),'patients':patients.patient_id.nunique(),'spots':len(allspots),'matched_IC':sum(g['matched_IC'] for g in graphs),'statistics_rows':len(stats),'nominal':int(stats.nominal.sum()),'FDR':int(stats.FDR.lt(.05).sum()),'geometry_endpoint':'unassessable','immune_exclusion_not_directly_measured':True,'actual_patient_statistics':True,'primary':stats[stats.program.eq('NABA_CORE_MATRISOME')&stats.endpoint.eq('CD8_adjusted')].to_dict('records'),'protocol_SHA256':sha(O/'resources/FUNCTIONAL_SIGNATURES_FROZEN.json'),'no_proxy':True}
assert qa['samples']==26 and qa['patients']==19 and qa['spots']==53309
dump(O/'delivery/SPATIAL_QA.json',qa);(O/'SPATIAL_COMPLETE').write_text('SUCCESS');dump(O/'SPATIAL_STATUS.json',{'state':'complete',**qa});print(json.dumps(qa),flush=True)
