from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,itertools,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase30_two_layer_classification_20260926'
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(f,x):
 t=f.with_name(f.name+'.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS.json',{'status':s,'at':now(),**kw});print(s,kw,flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def flips(d):
 n,g=d.shape
 if n<2:return np.full(g,np.nan),'not_tested_no_biological_replication'
 if n<=12:
  signs=np.array(list(itertools.product([-1,1],repeat=n)),dtype=np.float64);method='exact_paired_signflip'
 else:
  signs=np.random.default_rng(20260926).choice([-1.,1.],size=(19999,n));method='paired_signflip_MC20000'
 obs=np.abs(d.mean(0));p=np.empty(g)
 for start in range(0,g,256):
  block=d[:,start:start+256];null=np.abs(signs@block/n);ex=(null>=obs[start:start+256][None,:]-1e-12).sum(0)
  p[start:start+256]=(ex/len(signs)) if n<=12 else (ex+1)/(len(signs)+1)
 return p,method
def family_fdr(d):
 d=d.copy();valid=d.P.notna()
 q=multipletests(d.P.fillna(1),method='fdr_bh')[1];d['FDR']=np.where(valid,q,np.nan)
 d['nominal_exploratory']=valid&(d.P<.05)&(d.median_delta>0)
 return d
SURF=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
annotations=pd.read_csv(O/'tables/class_annotation_frozen.tsv',sep='\t')
prot=json.loads((O/'resources/CLASSIFICATION_AND_STATISTICS_PROTOCOL.json').read_text())
assert hashlib.sha256((O/'tables/class_annotation_frozen.tsv').read_bytes()).hexdigest()==prot['mapping_sha256']
rc=1
try:
 assert (O/'models.exit_code').read_text().strip()=='0'
 programs_all=[];surface_tests=[];top_markers=[];qa_rows=[];labelsets={};classprofiles=[];paired_roster=[];models={}
 for layer in ['supported','expanded']:
  a=ad.read_h5ad(O/'objects'/f'{layer}_integrated.h5ad')
  am=annotations[annotations.layer==layer].set_index('cluster');am.index=am.index.astype(str)
  for col in ['class_id','class_name','state_tags']:a.obs[col]=a.obs.cluster.astype(str).map(am[col])
  assert a.obs.class_id.notna().all()
  a.obs.to_csv(O/'tables'/f'{layer}_classified_cells.tsv.gz',sep='\t',index_label='cell_id');labelsets[layer]=a.obs.copy();models[layer]=a
  save(a.obs.groupby(['class_id','class_name','dataset','sample','donor_record','site','clinical_group'],observed=True).size().reset_index(name='cells'),layer+'_class_sample_counts.tsv')
  programs=[c for c in a.obs if c.startswith('program_')]
  classes=sorted(a.obs.class_id.unique())
  for cls,g in a.obs.groupby('class_id',observed=True):
   pr={'layer':layer,'class_id':cls,'class_name':g.class_name.iloc[0],'cells':len(g),'sample_records':g['sample'].nunique(),'donor_records':g.donor_record.nunique(),'dual_reference_cells':int(g.dual_reference.sum()),'supported_cells':int(g.in_supported.sum())}
   for panel in programs:pr[panel]=float(g.groupby('donor_record',observed=True)[panel].mean().mean())
   classprofiles.append(pr)
  pools={ds:[ds] for ds in a.obs.dataset.unique()}
  pools['verified_independent_GEO3']=['GSE152048','GSE162454','GSE318840']
  poolsummary=[]
  for pool,datasets in pools.items():
   status('PATIENT_LEVEL_MARKER_TESTS',layer=layer,pool=pool)
   ob=a.obs[a.obs.dataset.isin(datasets)];pool_genes=[];pool_programs=[]
   for cls in classes:
    gd=[];pdiff=[];local_roster=[]
    for dn,g in ob.groupby('donor_record',observed=True):
     gene_sites=[];program_sites=[];nt=nr=0;sites=[]
     for site,z in g.groupby('site',observed=True):
      target=z.class_id.eq(cls);other=~target
      if not target.any() or not other.any():continue
      it=a.obs.index.get_indexer(z.index[target]);ir=a.obs.index.get_indexer(z.index[other])
      tc=np.asarray(a.layers['counts'][it].sum(0)).ravel();oc=np.asarray(a.layers['counts'][ir].sum(0)).ravel()
      td=float(a.obs.iloc[it].full_library_counts.sum());od=float(a.obs.iloc[ir].full_library_counts.sum())
      gene_sites.append(np.log2(tc/td*1e4+.1)-np.log2(oc/od*1e4+.1))
      program_sites.append(z.loc[target,programs].mean().to_numpy()-z.loc[other,programs].mean().to_numpy())
      nt+=len(it);nr+=len(ir);sites.append(str(site))
     if gene_sites:
      gd.append(np.mean(gene_sites,axis=0));pdiff.append(np.mean(program_sites,axis=0))
      local_roster.append({'layer':layer,'pool':pool,'class_id':cls,'donor_record':dn,'dataset':g.dataset.iloc[0],'matched_sites':';'.join(sites),'target_cells':nt,'other_cells':nr,'matched_site_count':len(sites)})
    n=len(gd);paired_roster.extend(local_roster)
    for kind,features,values in [('gene',list(a.var_names),gd),('program',programs,pdiff)]:
     if n:
      d=np.vstack(values).astype(float);pv,method=flips(d);med=np.median(d,axis=0);mean=d.mean(0);pos=(d>0).sum(0);neg=(d<0).sum(0)
     else:
      d=np.empty((0,len(features)));pv=np.full(len(features),np.nan);med=mean=np.full(len(features),np.nan);pos=neg=np.zeros(len(features),int);method='not_tested_no_biological_replication'
     frame=pd.DataFrame({'layer':layer,'pool':pool,'class_id':cls,'feature':features,'n_independent_donors':n,'median_delta':med,'mean_delta':mean,'positive_donors':pos,'negative_donors':neg,'P':pv,'test':method,'family':layer+'|'+pool+'|'+kind+'_class_hypotheses'})
     if kind=='gene':
      frame['symbol']=a.var.symbol.to_numpy();frame['predefined_surface_candidate']=frame.symbol.isin(SURF)
      pool_genes.append(frame)
     else:pool_programs.append(frame)
    if n:
     # Complete patient-level program effects, with repeated samples aggregated inside donor/site.
     pr=pd.DataFrame(pdiff,columns=programs);pr=pd.concat([pd.DataFrame(local_roster).reset_index(drop=True),pr],axis=1)
     save(pr,f'{layer}_{pool}_{cls}_patient_program_effects.tsv')
   gf=family_fdr(pd.concat(pool_genes,ignore_index=True));pf=family_fdr(pd.concat(pool_programs,ignore_index=True))
   save(gf,f'{layer}_{pool}_all_gene_statistics.tsv.gz');programs_all.append(pf)
   surface_tests.append(gf[gf.predefined_surface_candidate].copy())
   for cls,g in gf.groupby('class_id'):
    positive=g[g.median_delta>0].sort_values(['nominal_exploratory','P','median_delta'],ascending=[False,True,False])
    top_markers.extend(positive.head(30).to_dict('records'))
   poolsummary.append({'layer':layer,'pool':pool,'gene_hypotheses':len(gf),'gene_tests_with_replication':int(gf.P.notna().sum()),'nominal_positive_gene_markers':int(gf.nominal_exploratory.sum()),'FDR_positive_gene_markers':int(((gf.FDR<.05)&(gf.median_delta>0)).sum()),'program_tests_with_replication':int(pf.P.notna().sum()),'nominal_positive_programs':int(pf.nominal_exploratory.sum()),'min_gene_P':float(gf.P.min()) if gf.P.notna().any() else None,'max_independent_donors':int(gf.n_independent_donors.max())})
   del gf;gc.collect()
  qa_rows.extend(poolsummary)
 save(pd.concat(programs_all,ignore_index=True),'class_program_statistics.tsv')
 save(pd.concat(surface_tests,ignore_index=True),'class_surface_gene_statistics.tsv')
 save(pd.DataFrame(top_markers),'class_top_marker_statistics.tsv')
 save(pd.DataFrame(paired_roster),'independent_patient_comparison_roster.tsv')
 save(pd.DataFrame(classprofiles),'class_profiles.tsv')
 save(pd.DataFrame(qa_rows),'statistics_summary.tsv')
 shared=models['supported'].obs.index
 cross=pd.crosstab(models['supported'].obs.class_id,models['expanded'].obs.loc[shared,'class_id'])
 cross.to_csv(O/'tables/class_correspondence_counts.tsv',sep='\t')
 save(models['expanded'].obs[models['expanded'].obs.dual_reference].groupby(['class_id','dataset','donor_record'],observed=True).size().reset_index(name='dual_reference_cells'),'dual114_reference_distribution.tsv')
 # Fit the newly classified contractile candidate population independently in both layers; never impose the old 230-cell ceiling.
 ns={};src=(O/'code/layer_models.py').read_text().split("lock=(O/")[0];exec(compile(src,str(O/'code/layer_models.py'),'exec'),ns)
 myqa=[]
 for layer,a in models.items():
  b=a[a.obs.class_id.eq('F2')].copy();b.X=b.layers['counts'].copy();b.obs['parent_layer_cluster']=b.obs.cluster.astype(str)
  label=layer+'_myCAF';status('MYCAF_INTERNAL_REFIT',layer=layer,cells=b.n_obs)
  b=ns['model'](b,label);pr=ns['profile'](b,label);ns['plot'](b,label,pr);ns['clean'](b)
  b.write_h5ad(O/'objects'/f'{label}_integrated.h5ad',compression='gzip');b.obs.to_csv(O/'tables'/f'{label}_cell_labels.tsv.gz',sep='\t',index_label='cell_id')
  myqa.append({'layer':layer,'myCAF_like_cells':b.n_obs,'internal_expression_clusters':b.obs.cluster.nunique(),'not_a_validated_myCAF_lineage':True})
 # Screen fixed surface candidates against the full called-cell background. Measurements do not imply protein/FACS purity.
 gate_defs=[(s,) for s in SURF]+list(itertools.combinations(SURF,2));gate_rows=[]
 for sm in sorted(models['expanded'].obs['sample'].unique()):
  bg=ad.read_h5ad(P/'phase27_identity_review_20260926/objects'/f'{sm}_marker_counts.h5ad')
  measured=set(bg.uns['measured_genes']);genes=list(bg.var_names);raw=bg.X.toarray()
  assert set(['PTPRC','PECAM1','EPCAM'])<=measured
  dump=(raw[:,[genes.index(g) for g in ['PTPRC','PECAM1','EPCAM']]]==0).all(1)
  valid_defs=[g for g in gate_defs if set(g)<=measured]
  gates=np.column_stack([dump&(raw[:,[genes.index(z) for z in g]]>0).all(1) for g in valid_defs])
  risk=bg.obs.CNV_evidence.isin(['method_conflict','single_method_alteration','concordant_alteration']).to_numpy()
  unassessed=bg.obs.CNV_evidence.eq('unassessed').to_numpy();mural=bg.obs.provisional_lineage.eq('Mural').to_numpy()
  ng=gates.sum(0);nr=(gates&risk[:,None]).sum(0);nu=(gates&unassessed[:,None]).sum(0);nm=(gates&mural[:,None]).sum(0)
  for layer,ob in labelsets.items():
   ids=ob.index[ob['sample'].eq(sm)];labels=ob.class_id.reindex(bg.obs_names).fillna('outside_research_layer')
   for cls in sorted(ob.loc[ids,'class_id'].unique()):
    target=labels.eq(cls).to_numpy();nt=int(target.sum());hit=(gates&target[:,None]).sum(0)
    for j,g in enumerate(valid_defs):
     gate_rows.append(dict(layer=layer,sample=sm,dataset=bg.obs.dataset.iloc[0],donor_record=bg.obs.donor_record.iloc[0],class_id=cls,markers='+'.join(g),full_called_cells=len(bg),target_cells=nt,RNA_gate_cells=int(ng[j]),target_in_gate=int(hit[j]),CNV_alteration_or_conflict_in_gate=int(nr[j]),CNV_unassessed_in_gate=int(nu[j]),mural_in_gate=int(nm[j])))
 gt=pd.DataFrame(gate_rows);save(gt,'surface_full_background_sample_counts.tsv.gz')
 cols=['full_called_cells','target_cells','RNA_gate_cells','target_in_gate','CNV_alteration_or_conflict_in_gate','CNV_unassessed_in_gate','mural_in_gate']
 dt=gt.groupby(['layer','dataset','donor_record','class_id','markers'],as_index=False)[cols].sum()
 dt['target_coverage']=dt.target_in_gate/dt.target_cells;dt['target_fraction_of_RNA_gate']=dt.target_in_gate/dt.RNA_gate_cells.replace(0,np.nan);dt['CNV_risk_fraction']=dt.CNV_alteration_or_conflict_in_gate/dt.RNA_gate_cells.replace(0,np.nan)
 save(dt,'surface_full_background_donor_screen.tsv')
 rank=dt.groupby(['layer','class_id','markers'],as_index=False).agg(donor_records=('donor_record','size'),median_target_coverage=('target_coverage','median'),median_target_fraction=('target_fraction_of_RNA_gate','median'),median_CNV_risk_fraction=('CNV_risk_fraction','median'))
 rank['descriptive_rank_score']=rank.median_target_coverage*rank.median_target_fraction
 save(rank.sort_values(['layer','class_id','descriptive_rank_score'],ascending=[True,True,False]),'surface_class_candidate_ranking.tsv')
 qa={'at':now(),'classification_frozen_before_tests':True,'supported_cells':1224,'expanded_cells':8174,'reference_only_cells':114,'myCAF_refits':myqa,'statistics_pools':qa_rows,'class_scheme':'F1 inflammatory/adventitial; F2 contractile-active matrix; F3 matrix remodeling; F4 osteogenic-matrix overlap; F5 mixed/unresolved. Technical states and CNV identity grades remain separate.','nominal_significance':'P<0.05 with positive median effect identifies exploratory candidates, full BH FDR retained','clinical_metastasis_specificity_established':False,'no_additional_cell_entry_filter':True}
 writej(O/'delivery/CLASSIFICATION_QA.json',qa);(O/'CLASSIFICATION_COMPLETE').write_text(now()+chr(10));status('TWO_LAYER_CLASSIFICATION_AND_MARKER_EXPLORATION_COMPLETE',supported_cells=1224,expanded_cells=8174,myCAF_refits=myqa);rc=0
except Exception as ex:status('FAILED',error=str(ex));traceback.print_exc()
finally:(O/'classification.exit_code').write_text(str(rc)+chr(10))
raise SystemExit(rc)
