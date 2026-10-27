from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,gzip,gc,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase53_de_novo_reannotation_20261002';B=S/'blind_review';V=P/'review_intake_20261004/received_review';O=P/'phase54_raw_review_20261004'
for d in ['code','tables','logs','delivery','resources','cells']: (O/d).mkdir(parents=True,exist_ok=True)
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8388608),b''):h.update(c)
 return h.hexdigest()
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
panels={'Fibro':['COL3A1','DCN','LUM','PDGFRA','COL6A1','COL6A2','DPT','CFD','PI16'],'Osteo':['RUNX2','SP7','ALPL','BGLAP','IBSP','IFITM5'],'Mural':['RGS5','CSPG4','NOTCH3','MCAM','MYH11','CNN1'],'Myeloid':['LST1','TYROBP','FCER1G','C1QA','CD163','CSF1R'],'Endo':['PECAM1','VWF','KDR','EMCN','CDH5','RAMP2'],'Epi':['SFTPB','SFTPC','SLC34A2','EPCAM','KRT8','KRT18'],'Platelet':['PF4','PPBP','GP9','ITGA2B'],'Chondro':['SOX9','ACAN','COL2A1','COMP','COL9A1']}
protocol={'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Raw-count and reused DecontX/available native doublet review; all 107 samples; 15 Fibro, 14 Mesench and all mixed-program targets. No cell exclusion or formal label replacement.','panels':panels,'held_out_identity_genes':['COL1A1','LOX','MMP14'],'feature_policy':'Lexicographically first source-common feature per symbol for fixed panels; frozen pair feature chosen from received evidence by matching marginal detection. All mapping exported.','codetection':'Observed raw joint detection; corrected >0 and >=1/2 counts are descriptive sensitivity, not cell-entry or biological classification gates.','statistics':'Descriptive identity review only; no inferential P values, cells are not biological replicates. Patient summaries within source; no assumed cross-source independence.','ambient':'Reuse completed DecontX output from phase31, automatic broad clustering without empty droplets; sensitivity evidence, not proof of ambient origin.','doublets':'Reuse completed native model or available upstream input obs evidence; missing remains unassessed.','CNV':'Same-sample continuous RNA CNV reference and osteogenic-expression comparison; exploratory similarity, not DNA-confirmed clone.','versions':{'numpy':np.__version__,'pandas':pd.__version__,'anndata':ad.__version__}}
dump('resources/PROTOCOL_FROZEN.json',protocol)
comp=read(V/'unblinded_comparison/cluster_annotation_comparison.tsv');details=read(V/'blind_submission/cluster_review_details.tsv');pairs=read(V/'blind_submission/evidence/mixed_program_evidence.tsv');ev=read(V/'blind_submission/evidence/selected_gene_evidence.tsv.gz')
targets=set(details.loc[details.lineage_code.isin(['Fibro','Mesench'])|details.lineage_code.str.startswith('Mix'),'cluster_id'])|set(pairs.cluster_id)
write(comp[comp.cluster_id.isin(targets)],'review_target_registry.tsv')
roster=read(B/'anonymous_cell_cluster_roster.tsv.gz').set_index('cell_id');key=read(S/'unblinding/cell_identifier_key.tsv.gz').set_index('cell_id');tech=read(B/'sample_technical_metadata.tsv');manifest=read(B/'RAW_COUNT_MANIFEST.tsv').set_index('relative_path')
selected={};frozen_maps=[]
for pair in pairs.itertuples(index=False):
 for gene,det in [(pair.gene_A,pair.detection_A),(pair.gene_B,pair.detection_B)]:
  v=ev[(ev.cluster_id==pair.cluster_id)&(ev.gene==gene)].copy();assert len(v)
  v['distance']=abs(v.detection_fraction-det);v=v.sort_values(['distance','feature_id']);assert v.iloc[0].distance<1e-10
  selected[(pair.cluster_id,gene)]=v.iloc[0].feature_id
  frozen_maps.append({'cluster_id':pair.cluster_id,'gene':gene,'feature_id':v.iloc[0].feature_id,'frozen_detection':det,'features_with_same_detection':int((v.distance<1e-10).sum())})
write(pd.DataFrame(frozen_maps),'frozen_pair_feature_mapping.tsv')
common={};panelmap={}
for source in tech.source_id.unique():
 g=read(B/source/'genes.tsv').fillna('');common[source]=set(g.feature)
 panelmap[source]={gene:sorted(g.loc[g.symbol.eq(gene),'feature'].astype(str))[0] for gene in sorted(set(sum(panels.values(),[]))) if g.symbol.eq(gene).any()}
write(pd.DataFrame([{'source_id':s,'gene':g,'feature_id':f} for s,m in panelmap.items() for g,f in m.items()]),'fixed_panel_feature_mapping.tsv')
allgene=[];allpairs=[];allqc=[];integrity=[];missingfeatures=[];d08diff=[];d08features=[];d08agg=None
for ti,t in enumerate(tech.sort_values('sample_id').itertuples(index=False)):
 sid=t.sample_id;source=t.source_id;f=B/'raw_counts'/f'{sid}.h5ad';mr=manifest.loc[f'raw_counts/{sid}.h5ad'];assert sha(f)==mr.SHA256
 a=ad.read_h5ad(f);x=a.X.tocsr().astype(np.float64);obs=roster.loc[a.obs_names].copy();assert obs.sample_id.eq(sid).all() and len(obs)==t.cells
 k=key.loc[a.obs_names];sample=k.original_sample.iloc[0];obs['original_cell_id']=k.original_cell_id.to_numpy();obs['tissue']=a.obs.tissue.to_numpy();obs['technology']=a.obs.technology.to_numpy()
 rawlib=np.asarray(x.sum(1)).ravel();assert np.array_equal(rawlib,a.obs.library_counts.to_numpy())
 obs['raw_library']=rawlib;obs['raw_detected_features']=x.getnnz(axis=1)
 inputdir=P/'phase26_client_reanalysis_20260922/inputs'/sample;ic=read(inputdir/'cells.tsv');ig=read(inputdir/'genes.tsv');assert ic.cell_id.tolist()==obs.original_cell_id.tolist() and ig.ensembl.tolist()==a.var_names.tolist()
 # Input obs can contain native provider doublet metadata; preserve supplied fields only.
 upstream=ad.read_h5ad(inputdir/'full_called_counts.h5ad',backed='r');dc=[c for c in upstream.obs.columns if 'doublet' in c.lower() or 'scrublet' in c.lower()]
 for c in dc:obs['upstream_'+c]=upstream.obs[c].to_numpy()
 upstream.file.close()
 de=P/'phase31_customer_completion_20260926/decontx'/sample;ok=(de/'COMPLETE').is_file() and (de/'exit_code').read_text().strip()=='0';assert ok
 cont=read(de/'cell_contamination.tsv');assert cont.cell_id.tolist()==obs.original_cell_id.tolist();assert np.array_equal(cont.original_UMI.to_numpy(),rawlib)
 with gzip.open(de/'corrected_counts.mtx.gz','rb') as h:
  with threadpool_limits(limits=4):y=mmread(h).tocsr().T.tocsr()
 assert y.shape==x.shape and np.isfinite(y.data).all() and (y.data>=0).all()
 corlib=np.asarray(y.sum(1)).ravel();assert np.allclose(corlib,cont.corrected_UMI,rtol=1e-5,atol=1e-4)
 obs['corrected_library']=corlib;obs['decontx_contamination']=cont.contamination.to_numpy()
 dbl=P/'phase31_customer_completion_20260926/doublet'/sample
 obs['native_doublet_class']='unassessed';obs['native_doublet_score']=np.nan
 if (dbl/'COMPLETE').is_file() and (dbl/'exit_code').read_text().strip()=='0':
  dr=read(dbl/'doublet_results.tsv').set_index('cell_id');assert dr.index.is_unique
  obs['native_doublet_class']=dr.scDblFinder_class.reindex(obs.original_cell_id).to_numpy();obs['native_doublet_score']=dr.scDblFinder_score.reindex(obs.original_cell_id).to_numpy();assert obs.native_doublet_class.notna().all()
 fmap={v:i for i,v in enumerate(a.var_names)};genes=list(panelmap[source]);idx=[fmap[panelmap[source][g]] for g in genes];rr=x[:,idx].toarray();cc=y[:,idx].toarray()
 rn=np.log1p(rr/np.maximum(rawlib,1)[:,None]*1e4);cn=np.log1p(cc/np.maximum(corlib,1)[:,None]*1e4)
 for name,gs in panels.items():
  cols=[genes.index(g) for g in gs if g in genes]
  obs['raw_'+name+'_mean_logCP10k']=rn[:,cols].mean(1) if cols else np.nan
  obs['corrected_'+name+'_mean_logCP10k']=cn[:,cols].mean(1) if cols else np.nan
  obs['raw_'+name+'_detected_genes']=(rr[:,cols]>0).sum(1) if cols else np.nan
  obs['corrected_'+name+'_genes_ge1']=(cc[:,cols]>=1).sum(1) if cols else np.nan
 sel=obs.r10.isin(targets).to_numpy();targetobs=obs.loc[sel].copy()
 for j,g in enumerate(genes):targetobs['raw_'+g]=rr[sel,j];targetobs['corrected_'+g]=cc[sel,j]
 targetobs.to_csv(O/'cells'/f'{sid}_target_cells.tsv.gz',sep='\t',index_label='cell_id')
 for cl in obs.r10.unique():
  mask=obs.r10.eq(cl).to_numpy();rest=~mask;n=int(mask.sum());is_target=cl in targets
  if is_target:
   for j,g in enumerate(genes):
    rcp=rr[mask,j].sum()/max(rawlib[mask].sum(),1)*1e4;rcpr=rr[rest,j].sum()/max(rawlib[rest].sum(),1)*1e4
    ccp=cc[mask,j].sum()/max(corlib[mask].sum(),1)*1e4;ccpr=cc[rest,j].sum()/max(corlib[rest].sum(),1)*1e4
    allgene.append({'source_id':source,'sample_id':sid,'patient_id':t.patient_id,'tissue':t.tissue,'cluster_id':cl,'cells':n,'gene':g,'feature_id':panelmap[source][g],'raw_total':rr[mask,j].sum(),'corrected_total':cc[mask,j].sum(),'raw_detection':(rr[mask,j]>0).mean(),'corrected_detection_gt0':(cc[mask,j]>0).mean(),'corrected_detection_ge1':(cc[mask,j]>=1).mean(),'raw_CP10k':rcp,'corrected_CP10k':ccp,'raw_vs_rest_log2_effect':np.log2((rcp+.1)/(rcpr+.1)) if rest.any() else np.nan,'corrected_vs_rest_log2_effect':np.log2((ccp+.1)/(ccpr+.1)) if rest.any() else np.nan})
  q={'source_id':source,'sample_id':sid,'patient_id':t.patient_id,'cluster_id':cl,'cells':n,'target':is_target,'native_doublet_assessed':int(obs.loc[mask,'native_doublet_class'].ne('unassessed').sum()),'native_doublet_calls':int(obs.loc[mask,'native_doublet_class'].eq('doublet').sum())}
  for column in ['raw_library','corrected_library','raw_detected_features','decontx_contamination','native_doublet_score']+[c for c in obs if c.endswith('_mean_logCP10k')]:
   q[column+'_median']=obs.loc[mask,column].median()
  allqc.append(q)
 for pair in pairs[pairs.cluster_id.isin(obs.r10.unique())].itertuples(index=False):
  mask=obs.r10.eq(pair.cluster_id).to_numpy();i=fmap[selected[(pair.cluster_id,pair.gene_A)]];j=fmap[selected[(pair.cluster_id,pair.gene_B)]]
  rv=x[mask,:][:,[i,j]].toarray();cv=y[mask,:][:,[i,j]].toarray();z={'source_id':source,'sample_id':sid,'patient_id':t.patient_id,'cluster_id':pair.cluster_id,'gene_A':pair.gene_A,'gene_B':pair.gene_B,'feature_A':a.var_names[i],'feature_B':a.var_names[j],'cells':int(mask.sum()),'raw_A_detected':int((rv[:,0]>0).sum()),'raw_B_detected':int((rv[:,1]>0).sum()),'raw_both':int((rv>0).all(1).sum()),'raw_neither':int((rv==0).all(1).sum()),'corrected_both_gt0':int((cv>0).all(1).sum()),'corrected_both_ge1':int((cv>=1).all(1).sum()),'corrected_both_ge2':int((cv>=2).all(1).sum()),'raw_A_total':rv[:,0].sum(),'raw_B_total':rv[:,1].sum(),'corrected_A_total':cv[:,0].sum(),'corrected_B_total':cv[:,1].sum()}
  allpairs.append(z)
 if source=='D08':
  cg=read(B/source/'genes.tsv');cls=read(B/source/'clusters.tsv');clcol='cluster_id' if 'cluster_id' in cls else 'cluster';cllist=cls[clcol].tolist()
  commonidx=[fmap[v] for v in cg.feature];non=np.flatnonzero(~a.var_names.isin(cg.feature));commonlib=np.asarray(x[:,commonidx].sum(1)).ravel();diff=rawlib-commonlib
  assert np.array_equal(diff,np.asarray(x[:,non].sum(1)).ravel())
  d08diff.append({'sample_id':sid,'patient_id':t.patient_id,'cells':len(obs),'full_counts':rawlib.sum(),'common_counts':commonlib.sum(),'excluded_feature_counts':diff.sum(),'excluded_measured_features':len(non)})
  sums=np.asarray(x[:,non].sum(0)).ravel()
  for j,total in zip(non,sums):
   if total>0:d08features.append({'sample_id':sid,'feature_id':a.var_names[j],'gene_symbol':a.var.iloc[j].gene_symbol,'counts':total})
  if d08agg is None:d08agg=np.zeros((len(cg),len(cllist)),dtype=np.float64)
  for ci,cl in enumerate(cllist):d08agg[:,ci]+=np.asarray(x[obs.r10.eq(cl).to_numpy()][:,commonidx].sum(0)).ravel()
 integrity.append({'sample_id':sid,'original_sample':sample,'source_id':source,'cells':len(obs),'raw_SHA256_verified':True,'raw_cell_gene_order_matches_correction_input':True,'corrected_library_matches':True,'target_cells':int(sel.sum()),'native_doublet_status':'available' if obs.native_doublet_class.ne('unassessed').all() else 'unassessed','upstream_doublet_fields':';'.join(dc)})
 print(json.dumps({'sample':sid,'done':ti+1,'of':len(tech),'target_cells':int(sel.sum())}),flush=True)
 dump('PROGRESS.json',{'state':'raw_and_corrected_review_running','completed_samples':ti+1,'total_samples':len(tech)})
 del a,x,y,rr,cc,rn,cn,obs,targetobs;gc.collect()
write(pd.DataFrame(allgene),'sample_gene_before_after.tsv.gz');write(pd.DataFrame(allpairs),'sample_pair_codetection.tsv');write(pd.DataFrame(allqc),'all_cluster_sample_QC_panels.tsv.gz');write(pd.DataFrame(integrity),'sample_integrity.tsv');write(pd.DataFrame(d08diff),'D08_count_reconciliation.tsv');write(pd.DataFrame(d08features),'D08_excluded_feature_counts.tsv')
with gzip.open(B/'D08/primary_cluster_counts.mtx.gz','rb') as h:
 with threadpool_limits(limits=4):original=mmread(h).toarray()
assert np.array_equal(original,d08agg);assert pd.DataFrame(d08diff).excluded_feature_counts.sum()==122334
pairsdf=pd.DataFrame(allpairs);num=['cells','raw_A_detected','raw_B_detected','raw_both','raw_neither','corrected_both_gt0','corrected_both_ge1','corrected_both_ge2','raw_A_total','raw_B_total','corrected_A_total','corrected_B_total'];ps=pairsdf.groupby(['cluster_id','gene_A','gene_B'])[num].sum().reset_index()
ps=ps.merge(pairs[['cluster_id','gene_A','gene_B','detection_A','detection_B','minimum_codetection_fraction']],on=['cluster_id','gene_A','gene_B'],validate='one_to_one')
assert np.allclose(ps.raw_A_detected/ps.cells,ps.detection_A,atol=1e-10) and np.allclose(ps.raw_B_detected/ps.cells,ps.detection_B,atol=1e-10)
for n in ['raw_both','corrected_both_gt0','corrected_both_ge1','corrected_both_ge2']:ps[n+'_fraction']=ps[n]/ps.cells
ps['A_count_retention']=ps.corrected_A_total/ps.raw_A_total;ps['B_count_retention']=ps.corrected_B_total/ps.raw_B_total;write(ps,'cluster_pair_codetection.tsv')
# One summary row per patient; multiple technical samples are averaged without pretending to add independent replication.
gd=pd.DataFrame(allgene);pdg=gd.groupby(['source_id','cluster_id','patient_id','gene'])[['raw_vs_rest_log2_effect','corrected_vs_rest_log2_effect','raw_detection','corrected_detection_ge1']].mean().reset_index();write(pdg,'patient_gene_direction.tsv.gz')
summary=pdg.groupby(['source_id','cluster_id','gene']).agg(patients=('patient_id','nunique'),raw_mean_effect=('raw_vs_rest_log2_effect','mean'),corrected_mean_effect=('corrected_vs_rest_log2_effect','mean'),raw_positive_patients=('raw_vs_rest_log2_effect',lambda x:int((x>0).sum())),corrected_positive_patients=('corrected_vs_rest_log2_effect',lambda x:int((x>0).sum()))).reset_index();write(summary,'cluster_gene_patient_summary.tsv')
dump('RAW_REVIEW_QA.json',{'samples':len(integrity),'cells':sum(x['cells'] for x in integrity),'target_clusters':len(targets),'target_cells':sum(x['target_cells'] for x in integrity),'raw_hashes_verified':len(integrity),'pair_marginal_detection_matches_frozen':True,'D08_all_common_feature_cluster_counts_match':True,'D08_excluded_counts':122334,'no_cell_filters':True,'formal_labels_unchanged':True})
(O/'RAW_REVIEW_COMPLETE').write_text('SUCCESS\n')
print('RAW_REVIEW_COMPLETE',flush=True)