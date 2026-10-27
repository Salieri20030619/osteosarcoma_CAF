# SYNTAX_CHECK: python -m py_compile stageD_clinical.py
# RUN_COMMAND: batch_bench/bin/python -u stageD_clinical.py
# EXPECTED_SUCCESS_MARKER: STAGE_D_CLINICAL_FRACTIONS_COMPLETE
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,datetime,traceback,re,ast,itertools,math,warnings
import numpy as np,pandas as pd,anndata as ad,statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';B=P/'phase42_F3F4_identity_states_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase44_patient_mechanism_20261001'
for sub in ['code','logs','tables','objects','resources','delivery','figures']:(O/sub).mkdir(parents=True,exist_ok=True)
def writej(f,x):
 t=f.with_suffix('.tmp');t.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str));t.replace(f)
def status(s,**kw):writej(O/'STATUS_CLINICAL.json',{'state':s,'pid':os.getpid(),'at':datetime.datetime.now().astimezone().isoformat(),**kw});print(s,kw,flush=True)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
source=C/'code/stageC_full_gene_statistics.py';nodes=[x for x in ast.parse(source.read_text()).body if isinstance(x,ast.FunctionDef) and x.name in ['mean_valid','tests','fdr_by_family']];signs_cache={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'))
rng=np.random.default_rng(20261001)
def independent(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float);a=a[np.isfinite(a)];b=b[np.isfinite(b)]
 if len(a)<2 or len(b)<2:return np.nan,'not_tested_insufficient_patients'
 z=np.r_[a,b]
 if np.all(z==z[0]):return np.nan,'not_tested_constant_outcome'
 obs=abs(b.mean()-a.mean());n=len(a);total=math.comb(len(z),n);hits=0
 if total<=200000:
  for ix in itertools.combinations(range(len(z)),n):
   chosen=z[list(ix)];other=(z.sum()-chosen.sum())/len(b);hits+=abs(other-chosen.mean())>=obs-1e-12
  return hits/total,'exact_patient_label_permutation'
 for _ in range(19999):
  ix=rng.choice(len(z),size=n,replace=False);chosen=z[ix];other=(z.sum()-chosen.sum())/len(b);hits+=abs(other-chosen.mean())>=obs-1e-12
 return (hits+1)/20000,'patient_label_permutation_MC19999_plus_observed'
def site_group(site):
 s=str(site).lower();lung='lung' in s;bone=bool(re.search(r'(?:\bbone\b|femur|tibia|fibula|humerus|radius|ulna|pelvi|ilium|sacr|spin|vertebr|clavicl|scapul|\brib\b)',s));mixed='chest' in s
 if lung and not bone and not mixed:return 'Lung'
 if bone and not lung and not mixed:return 'Bone'
 return 'Unknown_or_mixed'
rc=1
try:
 assert (A/'COMPLETE').read_text().strip()=='SUCCESS' and (B/'COMPLETE').read_text().strip()=='SUCCESS'
 a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=a.obs.copy();a.file.close();m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);assert m['sample'].is_unique
 cov=pd.read_csv(A/'tables/source_sample_gene_coverage.tsv',sep='\t').set_index('sample');states=pd.read_csv(B/'tables/all_F3F4_internal_cell_states.tsv.gz',sep='\t',keep_default_na=False)
 m['strict_site_group']=m.site.map(site_group)
 save(m,'clinical_manifest_with_strict_site_group.tsv')
 writej(O/'resources/CLINICAL_PROTOCOL_FROZEN.json',{'target':'all frozen operational classes and internal F3/F4 states; source-specific patient fractions','class_denominators':['all_frozen_candidates_in_layer','all_source_called_cells'],'internal_state_denominators':['all_frozen_candidates_in_layer','all_source_called_cells','parent_class_cells'],'aggregation':'sum sample counts within source/patient/event; each patient equal weight in contrasts','clinical_contrasts':['strict_Bone_vs_Lung','reporting_A_vs_B1','reporting_B1_vs_B2'],'reporting_groups':'original frozen labels and evidence retained; no inferred A, no unknown-to-nonmetastatic conversion','strict_site':'unambiguous bone anatomical labels or lung without bone/chest mixing; mixed/unknown not forced','paired':'second group minus first, exact patient sign-flip; paired patients removed from independent comparison','independent':'two-sided patient label permutation; exact up to 200000 allocations, MC19999+observed beyond','zero_handling':'zero numerator with valid denominator is real operational fraction; zero denominator NA; all-constant/zero effects not tested','adjusted_model':'optional independent-patient OLS fraction~target_group+treatment+technology, HC3; complete recorded covariates, full rank and residual df>=2; target_group P, separate adjusted family','family':'layer|source|contrast|design|clinical_proportions, all classes/states/denominators in same family; separate adjusted family','causality':'operational fraction association, no transfer/treatment or metastasis causal claim','no_global_cross_source_pooling':True})
 rows=[];defs=[]
 for layer in ['supported','expanded']:
  label='supported_class_id' if layer=='supported' else 'expanded_class_id';v=ob[ob.in_supported.astype(bool)] if layer=='supported' else ob
  st=states[states.layer.eq(layer)];state_map=st.set_index('cell_id').state_id.to_dict();specs=[(x,'class','') for x in ['F1','F2','F3','F4','F5']]+[(s,'internal',p) for s,p in st[['state_id','parent_class']].drop_duplicates().itertuples(index=False,name=None)]
  for feature,kind,parent in specs:
   for den in ['all_frozen_candidates_in_layer','all_source_called_cells']+(['parent_class_cells'] if kind=='internal' else []):defs.append({'layer':layer,'feature':feature,'kind':kind,'parent':parent,'denominator':den})
  for sample,meta in m.set_index('sample').iterrows():
   g=v[v['sample'].astype(str).eq(sample)];cls=g[label].astype(str);state=pd.Series([state_map.get(x,'') for x in g.index],index=g.index);called=int(cov.loc[sample,'called_cells'])
   for feature,kind,parent in specs:
    numerator=int(cls.eq(feature).sum()) if kind=='class' else int(state.eq(feature).sum())
    denominators={'all_frozen_candidates_in_layer':len(g),'all_source_called_cells':called}
    if kind=='internal':denominators['parent_class_cells']=int(cls.eq(parent).sum())
    for den,total in denominators.items():rows.append({'layer':layer,'sample':sample,'dataset':meta.dataset,'donor_record':meta.donor_record,'clinical_group':meta.clinical_group,'group_evidence':meta.group_evidence,'strict_site_group':meta.strict_site_group,'treatment':meta.treatment,'technology':meta.technology,'feature':feature,'kind':kind,'parent':parent,'denominator':den,'numerator':numerator,'denominator_count':total,'fraction':numerator/total if total else np.nan})
 sample=pd.DataFrame(rows);save(sample,'sample_operational_class_state_fractions.tsv.gz')
 allpatients=[];result=[];roster=[];adjusted=[];qa=[]
 planned=[('strict_Bone_vs_Lung','strict_site_group','Bone','Lung'),('reporting_A_vs_B1','clinical_group','A','B1'),('reporting_B1_vs_B2','clinical_group','B1','B2')]
 for layer in ['supported','expanded']:
  for ds in m.dataset.unique():
   g=sample[(sample.layer==layer)&(sample.dataset==ds)]
   for name,col,first,second in planned:
    x=g[g[col].isin([first,second])].copy();pt=[]
    for (donor,group,feature,den),z in x.groupby(['donor_record',col,'feature','denominator'],observed=True):
     total=z.denominator_count.sum();num=z.numerator.sum();treatment=';'.join(sorted(set(z.treatment)));technology=';'.join(sorted(set(z.technology)))
     pt.append({'layer':layer,'dataset':ds,'contrast':name,'donor_record':donor,'group':group,'feature':feature,'denominator':den,'numerator':num,'denominator_count':total,'fraction':num/total if total else np.nan,'samples':';'.join(z['sample']),'treatment':treatment,'technology':technology,'group_evidence':';'.join(sorted(set(z.group_evidence)))})
    if not pt:continue
    patient=pd.DataFrame(pt);allpatients.append(patient)
    for definition in [z for z in defs if z['layer']==layer]:
     feature=definition['feature'];den=definition['denominator'];z=patient[(patient.feature==feature)&(patient.denominator==den)&patient.fraction.notna()].copy()
     a=z[z.group.eq(first)].set_index('donor_record');b=z[z.group.eq(second)].set_index('donor_record');paired=sorted(set(a.index)&set(b.index))
     delta=b.loc[paired,'fraction'].to_numpy()-a.loc[paired,'fraction'].to_numpy() if paired else np.array([])
     p,N,method=tests(delta[:,None]);result.append({'layer':layer,'dataset':ds,'contrast':name,'design':'paired','feature':feature,'kind':definition['kind'],'parent':definition['parent'],'denominator':den,'N_first':len(paired),'N_second':len(paired),'N_pairs':len(paired),'mean_effect_second_minus_first':float(delta.mean()) if len(delta) else np.nan,'positive_patients':int((delta>0).sum()),'negative_patients':int((delta<0).sum()),'P':p[0],'test':method[0],'family':f'{layer}|{ds}|{name}|paired|clinical_proportions'})
     left=a.drop(index=paired);right=b.drop(index=paired);p,method=independent(left.fraction,right.fraction);effect=right.fraction.mean()-left.fraction.mean() if len(left) and len(right) else np.nan
     result.append({'layer':layer,'dataset':ds,'contrast':name,'design':'independent_paired_removed','feature':feature,'kind':definition['kind'],'parent':definition['parent'],'denominator':den,'N_first':len(left),'N_second':len(right),'N_pairs':0,'mean_effect_second_minus_first':effect,'positive_patients':np.nan,'negative_patients':np.nan,'P':p,'test':method,'family':f'{layer}|{ds}|{name}|independent|clinical_proportions'})
     for event,zz in [(first,a),(second,b)]:
      for donor,row in zz.iterrows():roster.append({'layer':layer,'dataset':ds,'contrast':name,'feature':feature,'denominator':den,'donor_record':donor,'group':event,'paired':donor in paired,'fraction':row.fraction,'samples':row.samples})
     # Independent multivariable model, no substitution of the minimum covariate P.
     model=pd.concat([left.assign(target_group=0),right.assign(target_group=1)]);reason='not_fit';coef=Pval=se=np.nan
     missing=model.treatment.astype(str).str.lower().isin(['','field_absent','na','unknown','not_reported','unavailable'])|model.technology.astype(str).str.lower().isin(['','na','unknown','not_reported','unavailable'])|model.treatment.str.contains(';')|model.technology.str.contains(';')
     if len(left)<2 or len(right)<2:reason='insufficient_group_replication'
     elif missing.any():reason='missing_or_mixed_recorded_covariates'
     else:
      covariate=pd.get_dummies(model[['treatment','technology']],drop_first=True,dtype=float);X=pd.concat([model[['target_group']].astype(float),covariate],axis=1);X=sm.add_constant(X,has_constant='add')
      if np.linalg.matrix_rank(X.to_numpy())<X.shape[1] or len(X)-X.shape[1]<2:reason='rank_or_residual_df_not_evaluable'
      elif model.fraction.nunique()<2:reason='constant_outcome'
      else:
       fit=sm.OLS(model.fraction.astype(float),X).fit(cov_type='HC3');coef=float(fit.params['target_group']);Pval=float(fit.pvalues['target_group']);se=float(fit.bse['target_group']);reason='fit_HC3_target_group'
     adjusted.append({'layer':layer,'dataset':ds,'contrast':name,'feature':feature,'kind':definition['kind'],'denominator':den,'N':len(model),'target_effect':coef,'SE_HC3':se,'P':Pval,'status':reason,'family':f'{layer}|{ds}|{name}|adjusted_clinical_proportions'})
     if len(paired):
      for donor,value in zip(paired,delta):qa.append({'comparison_key':f'{layer}|{ds}|{name}|{feature}|{den}','donor_record':donor,'paired_effect':value,'expected_P':result[-2]['P'],'expected_N':len(delta)})
    status('CLINICAL_SOURCE',layer=layer,dataset=ds,contrast=name,patient_event_rows=len(patient))
 if result:clinical=fdr_by_family(pd.DataFrame(result));save(clinical,'all_clinical_fraction_statistics.tsv');save(clinical[clinical.nominal],'nominal_clinical_fraction_statistics.tsv')
 else:clinical=pd.DataFrame()
 if adjusted:adj=fdr_by_family(pd.DataFrame(adjusted));save(adj,'adjusted_clinical_fraction_statistics.tsv')
 else:adj=pd.DataFrame()
 save(pd.concat(allpatients,ignore_index=True),'patient_event_class_state_fractions.tsv.gz');save(pd.DataFrame(roster),'clinical_patient_analysis_roster.tsv.gz');save(pd.DataFrame(qa),'independent_clinical_paired_QA_input.tsv')
 out={'at':datetime.datetime.now().astimezone().isoformat(),'sample_fraction_rows':len(sample),'clinical_test_rows':len(clinical),'valid_tests':int(clinical.P.notna().sum()) if len(clinical) else 0,'nominal_tests':int(clinical.P.lt(.05).sum()) if len(clinical) else 0,'FDR_tests':int(clinical.FDR.lt(.05).sum()) if len(clinical) else 0,'adjusted_fitted_tests':int(adj.P.notna().sum()) if len(adj) else 0,'same_patient_removed_from_independent':True,'zero_denominator_NA':True,'unknown_clinical_labels_not_inferred':True,'proportions_of_operational_classes_not_absolute_CAF_abundance':True,'C2_required':False}
 writej(O/'delivery/CLINICAL_FRACTION_QA.json',out);(O/'CLINICAL_COMPLETE').write_text('SUCCESS\n');status('complete',**out);print('STAGE_D_CLINICAL_FRACTIONS_COMPLETE',flush=True);rc=0
except Exception:traceback.print_exc();status('failed');(O/'CLINICAL_FAILED').write_text(traceback.format_exc())
finally:(O/'clinical.exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)