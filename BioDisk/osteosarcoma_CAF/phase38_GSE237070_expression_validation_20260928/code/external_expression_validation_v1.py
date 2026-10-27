from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,time,datetime,hashlib,traceback,fcntl,gc,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase37_requested_source_search_20260927';O=P/'phase38_GSE237070_expression_validation_20260928';B=P/'phase30_two_layer_classification_20260926'
for d in ['code','tables','resources','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
lock=(O/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'EXPRESSION_RESULTS_READY.json').exists()
def dump(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(p)
def status(state,**kw):dump(O/'STATUS.json',{'state':state,'at':datetime.datetime.now().astimezone().isoformat(),**kw})
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False)
lineage=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text())
scope=json.loads((S/'tables/GSE237070_VALIDATION_SCOPE.json').read_text());markers=scope['markers'];programs=scope['programs'];target=sorted(set(markers+sum(programs.values(),[])))
blind={k:[g for g in gs if g not in markers] for k,gs in lineage.items()}
protocol={'created':datetime.datetime.now().astimezone().isoformat(),'source':'GSE237070 reconstructed raw/native-called counts;6samples4patients;2primary+4adjacent','scope':'external expression evidence for candidate classification and analytical flow; not protein/functional validation','annotation':'reuse phase27 competitive mean log1p(CP10k) lineage panels, preserve exact tie/allzero as unresolved','validation_annotation':'also compute all-eight-surface-markers-held-out lineage scores to avoid using tested surface markers in their own grouping; no new expression or confidence cutoff','class_mapping':'separately map query fibroblast-like cells to donor-averaged phase30 supported and expanded frozen class prototypes using Pearson expression correlation on each layer frozen HVGs intersect query; exclude all59predeclared target genes from mapping features; no confidence cutoff; classes are reference-like phenotypes, not proven identities','unassessed_identity':'CNV and doublet status in this new source remain unassessed; do not assign it to the frozen supported1224 layer or call it confirmed CAF; adjacent fibroblasts are not called CAF','markers':markers,'programs':programs,'stats':'within each tissue, donor-paired group-vs-group effects; genes log2(raw group-pseudobulk CP10k+0.1); programs mean cellular log1p CP10k; exact two-sided sign flips; n<2 or all-zero endpoint -> P NA','families':'BH separately by annotation scheme, tissue and group-contrast, and by marker/program endpoint type; complete prespecified endpoint family, untestable rows retained','nominal_gate':'raw P<0.05 only; FDR reported not entry gate','sample_size':'primary N<=2, adjacent N<=4; exploratory directions/effects and consistency first; no cell-level biological P','coexpression':'raw RNA>0 descriptive co-detection only; not protein positivity or CAF inclusion gate','no_added_QC_threshold':True,'no_new_functional_or_CNV_prerequisite':True}
dump(O/'resources/PROTOCOL_FROZEN.json',protocol)
queries=sorted((S/'objects/GSE237070_reconstructed').glob('*/full_native_called_counts.h5ad'));assert len(queries)==6
first=ad.read_h5ad(queries[0],backed='r');qgenes=pd.Index(first.var_names.astype(str));first.file.close()
refs={};prototype_rows=[];coverage=[];account=[];allobs=[];aggregates=[];coexpr=[];source=[]
rc=1
try:
 status('preparing_frozen_reference_prototypes')
 for layer in ['supported','expanded']:
  a=ad.read_h5ad(B/'objects'/f'{layer}_integrated.h5ad');lab=pd.read_csv(B/'tables'/f'{layer}_classified_cells.tsv.gz',sep='\t').set_index('cell_id').loc[a.obs_names]
  features=[g for g in a.var_names[a.var.highly_variable] if g in qgenes and not any(z in target for z in str(a.var.loc[g,'symbol']).split(';'))]
  idx=a.var_names.get_indexer(features);x=a.X[:,idx].tocsr();classes=sorted(lab.class_id.unique());protos=[]
  for cl in classes:
   dmeans=[]
   for donor in lab.loc[lab.class_id.eq(cl),'donor_record'].unique():
    keep=lab.class_id.eq(cl)&lab.donor_record.eq(donor);dmeans.append(np.asarray(x[keep.to_numpy()].mean(0)).ravel())
   mean=np.mean(dmeans,axis=0);protos.append(mean);prototype_rows.append({'layer':layer,'class_id':cl,'reference_cells':int(lab.class_id.eq(cl).sum()),'donor_records':len(dmeans),'mapping_features':len(features)})
  pro=np.asarray(protos);pro-=pro.mean(1,keepdims=True);pro/=np.linalg.norm(pro,axis=1,keepdims=True)
  refs[layer]=(features,classes,pro);pd.DataFrame({'ensembl':features}).to_csv(O/'resources'/f'{layer}_mapping_features.tsv',sep='\t',index=False)
  del a,x;gc.collect()
 save(pd.DataFrame(prototype_rows),'reference_prototypes.tsv')
 for path in queries:
  a=ad.read_h5ad(path);sm=str(a.obs['sample'].iloc[0]);donor=str(a.obs.patient.iloc[0]);tissue=str(a.obs.tissue.iloc[0]);status('annotating_and_summarizing',sample=sm,finished_samples=len(source))
  x=a.X.tocsr();den=np.asarray(x.sum(1)).ravel();assert np.allclose(den,a.obs.n_counts) and np.isfinite(x.data).all() and (x.data>=0).all()
  sy=a.var.symbol.astype(str);genes=sorted(set(target+sum(lineage.values(),[])));gm={g:np.flatnonzero(sy.to_numpy()==g) for g in genes}
  raw=np.column_stack([np.asarray(x[:,gm[g]].sum(1)).ravel() if len(gm[g]) else np.full(a.n_obs,np.nan) for g in genes]);norm=raw/np.where(den>0,den,np.nan)[:,None]*1e4;lg=np.log1p(norm);gidx={g:i for i,g in enumerate(genes)}
  obs=a.obs.copy();obs.index.name='cell_id';obs['CNV_evidence']='unassessed_external_source';obs['doublet_evidence']='unassessed_external_source'
  for scheme,panels in [('frozen',lineage),('surface_markers_held_out',blind)]:
   scores=[]
   for name,gs in panels.items():
    have=[gidx[g] for g in gs if len(gm[g])];v=np.mean(lg[:,have],axis=1) if have else np.full(a.n_obs,np.nan);obs[f'{scheme}_score_{name}']=v;scores.append(v)
    coverage.append({'sample':sm,'scheme':scheme,'panel':name,'present_genes':len(have),'requested_genes':len(gs),'missing':';'.join(g for g in gs if not len(gm[g]))})
   values=np.column_stack(scores);v=np.where(np.isfinite(values),values,-np.inf);maximum=v.max(1);labels=np.asarray(list(panels),dtype=object)[v.argmax(1)];labels[(maximum<=0)|(np.sum(v==maximum[:,None],axis=1)!=1)]='Unresolved';obs[f'lineage_{scheme}']=labels
  pscores={}
  for name,gs in programs.items():
   have=[gidx[g] for g in gs if len(gm[g])];pscores[name]=np.mean(lg[:,have],axis=1) if have else np.full(a.n_obs,np.nan);obs['program_'+name]=pscores[name]
  fib=obs.lineage_surface_markers_held_out.eq('Fibroblast').to_numpy()
  for layer,(features,classes,pro) in refs.items():
   qi=a.var_names.get_indexer(features);v=x[fib][:,qi].astype(float).multiply((1e4/den[fib])[:,None]).tocsr();v.data=np.log1p(v.data)
   if v.shape[0]:
    variance=np.asarray(v.multiply(v).sum(1)).ravel()-np.asarray(v.sum(1)).ravel()**2/v.shape[1];norms=np.sqrt(np.maximum(variance,0));corr=np.asarray(v@pro.T)/np.where(norms>0,norms,np.nan)[:,None]
    label=np.asarray(classes,dtype=object)[np.nan_to_num(corr,nan=-np.inf).argmax(1)];mx=np.nanmax(corr,axis=1);label[(norms==0)|(np.sum(corr==mx[:,None],axis=1)!=1)]='Unresolved'
   else:corr=np.empty((0,len(classes)));label=np.array([],dtype=object)
   obs[f'{layer}_reference_like']='not_fibroblast_like';obs.loc[fib,f'{layer}_reference_like']=label
   for j,cl in enumerate(classes):obs[f'{layer}_correlation_{cl}']=np.nan;obs.loc[fib,f'{layer}_correlation_{cl}']=corr[:,j]
  obs.reset_index().to_csv(O/'tables'/f'{sm}_cell_annotations.tsv.gz',sep='\t',index=False)
  types={f'lineage_{s}':obs[f'lineage_{s}'].astype(str).to_numpy() for s in ['frozen','surface_markers_held_out']}
  for layer in refs:types[layer+'_reference_like']=obs[layer+'_reference_like'].astype(str).to_numpy()
  for grouping,labels in types.items():
   masks={g:labels==g for g in sorted(set(labels)) if g!='not_fibroblast_like'}
   if grouping.endswith('_reference_like'):
    classes=refs[grouping.replace('_reference_like','')][1]
    for cl in classes:masks['rest_of_'+cl]=np.isin(labels,classes)&(labels!=cl)
   for group,mask in masks.items():
    n=int(mask.sum())
    if not n:continue
    lib=den[mask].sum();row={'sample':sm,'patient':donor,'tissue':tissue,'grouping':grouping,'group':group,'cells':n,'library_UMI':lib}
    for gene in markers:
     rawsum=float(raw[mask,gidx[gene]].sum());row['count_'+gene]=rawsum;row['gene_'+gene]=float(np.log2(rawsum/lib*1e4+.1)) if lib>0 else np.nan;row['fraction_'+gene]=float(np.mean(raw[mask,gidx[gene]]>0))
    for name,pv in pscores.items():row['program_'+name]=float(pv[mask].mean())
    aggregates.append(row)
    for g1,g2 in [('PDGFRB','MCAM'),('GPC3','THY1')]:
     b1=raw[mask,gidx[g1]]>0;b2=raw[mask,gidx[g2]]>0;coexpr.append({'sample':sm,'patient':donor,'tissue':tissue,'grouping':grouping,'group':group,'cells':n,'pair':g1+'__'+g2,'both_detected':int((b1&b2).sum()),'fraction_both':float((b1&b2).mean()),'RNA_not_protein':True})
  for scheme in ['frozen','surface_markers_held_out']:
   for k,n in obs[f'lineage_{scheme}'].value_counts().items():account.append({'sample':sm,'patient':donor,'tissue':tissue,'scheme':scheme,'lineage':k,'cells':int(n)})
  source.append({'sample':sm,'path':str(path),'cells':a.n_obs,'features':a.n_vars,'all_native_cells_retained':True});del a,x,raw,norm,lg,obs;gc.collect()
 save(pd.DataFrame(account),'sample_lineage_counts.tsv');save(pd.DataFrame(coverage),'panel_gene_coverage.tsv');agg=pd.DataFrame(aggregates);save(agg,'sample_group_marker_program_values.tsv');save(pd.DataFrame(coexpr),'RNA_pair_codetection.tsv')
 status('patient_level_effects')
 effects=[];tests=[]
 for (grouping,tissue),frame in agg.groupby(['grouping','tissue']):
  comparisons=[('Fibroblast',b) for b in lineage if b!='Fibroblast'] if grouping.startswith('lineage_') else [(cl,'rest_of_'+cl) for cl in refs[grouping.replace('_reference_like','')][1]]
  for left,right in comparisons:
   joined=frame.loc[frame.group.eq(left)].merge(frame.loc[frame.group.eq(right)],on=['sample','patient','tissue'],suffixes=('_left','_right'));assert joined.patient.is_unique
   for kind,ends in [('gene',markers),('program',list(programs))]:
    for name in ends:
     col=kind+'_'+name;vals=joined[col+'_left'].to_numpy()-joined[col+'_right'].to_numpy();valid=np.isfinite(vals);vals=vals[valid];n=len(vals)
     zero=(joined['count_'+name+'_left'].to_numpy()[valid].sum()+joined['count_'+name+'_right'].to_numpy()[valid].sum()==0) if kind=='gene' else bool(n and np.all(joined[col+'_left'].to_numpy()[valid]==0)&np.all(joined[col+'_right'].to_numpy()[valid]==0))
     p=np.nan
     if n>=2 and not zero:
      signs=np.asarray(list(itertools.product([-1,1],repeat=n)));null=np.abs(signs@vals/n);stat=abs(vals.mean());p=float(np.mean(null>=stat-1e-12));assert p>=2/2**n-1e-12
     row={'grouping':grouping,'tissue':tissue,'left':left,'right':right,'kind':kind,'endpoint':name,'N_patients':n,'mean_effect':float(vals.mean()) if n else np.nan,'median_effect':float(np.median(vals)) if n else np.nan,'positive_patients':int((vals>0).sum()),'negative_patients':int((vals<0).sum()),'P':p,'untestable_reason':'all_zero' if zero else 'less_than_two_patients' if n<2 else '','family':f'{grouping}|{tissue}|{left}_vs_{right}|{kind}'};tests.append(row)
     for patient,value in zip(joined.patient.to_numpy()[valid],vals):effects.append({k:row[k] for k in ['grouping','tissue','left','right','kind','endpoint']}|{'patient':patient,'effect':float(value)})
 test=pd.DataFrame(tests);test['FDR']=np.nan
 for fam,indices in test.groupby('family').groups.items():
  v=test.loc[indices,'P'];q=multipletests(v.fillna(1),method='fdr_bh')[1];test.loc[indices,'FDR']=np.where(v.notna(),q,np.nan)
 test['nominal_P_lt_0_05']=test.P.lt(.05);save(test,'patient_level_marker_program_statistics.tsv');save(pd.DataFrame(effects),'individual_patient_effects.tsv')
 complete={'at':datetime.datetime.now().astimezone().isoformat(),'state':'expression_results_ready_for_review','samples':len(source),'patients':4,'primary_patients':2,'called_cells':sum(x['cells'] for x in source),'annotation_and_reference_mapping_complete':True,'statistics_rows':len(test),'nominal_rows':int(test.nominal_P_lt_0_05.sum()),'CNV_identity_not_newly_established':True,'plots_and_interpretation_review_pending':True,'source_cells':source}
 assert complete['called_cells']==47059;dump(O/'EXPRESSION_RESULTS_READY.json',complete);status('expression_results_ready_for_review',**{k:v for k,v in complete.items() if k not in ['state','at','source_cells']});rc=0
except Exception as e:status('failed',error=repr(e));traceback.print_exc();dump(O/'ERROR.json',{'error':repr(e),'traceback':traceback.format_exc()})
finally:(O/'analysis.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)
