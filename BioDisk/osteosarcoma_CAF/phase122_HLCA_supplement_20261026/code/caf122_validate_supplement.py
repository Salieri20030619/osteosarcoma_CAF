from pathlib import Path
import json,hashlib,zipfile,numpy as np,pandas as pd,anndata as ad
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
rec=json.loads((O/'ANALYSIS_COMPLETE.json').read_text());check('main_query6423',rec['query_cells']==6423)
q=ad.read_h5ad(O/'objects/query_program_counts.h5ad');r=ad.read_h5ad(O/'objects/reference_http_program_counts.h5ad');pred=pd.read_csv(O/'tables/HLCA_current_cell_predictions.tsv.gz',sep='\t');check('prediction_order',q.obs.cell_id.tolist()==pred.cell_id.tolist());check('unique_query_ids',q.obs.cell_id.is_unique)
for node,n in [('Lung_04',233),('Bone_07',131),('Bone_07_s1',35),('Lung_04_s2',36)]:check('node_count_'+node,int(q.obs['leaf' if '_s' in node else 'parent'].eq(node).sum())==n,n)
for label,a,den in [('query',q,q.obs.library_counts),('reference',r,r.obs.raw_sum)]:
 check(label+'_raw_integer',np.isfinite(a.X.data).all() and (a.X.data>=0).all() and np.equal(a.X.data,np.floor(a.X.data)).all());check(label+'_subset_within_library',np.all(np.asarray(a.X.sum(1)).ravel()<=np.asarray(den)+1e-5))
sc=pd.read_csv(O/'tables/query_cell_program_scores.tsv.gz',sep='\t');rc=pd.read_csv(O/'tables/reference_joint_cell_program_scores.tsv.gz',sep='\t');coverage=json.loads((O/'resources/JOINT_PROGRAM_COVERAGE.json').read_text())
for label,a,table,den in [('query',q,sc,q.obs.library_counts.to_numpy()),('reference',r,rc,r.obs.raw_sum.to_numpy())]:
 for name,entry in coverage.items():
  if name not in table:continue
  gs=entry['genes'];ix=a.var_names.get_indexer(gs);rows=np.unique(np.linspace(0,a.n_obs-1,min(20,a.n_obs),dtype=int));manual=np.log1p(a.X[rows,:][:,ix].toarray()*10000/den[rows,None]).mean(1)
  check(label+'_recomputed_'+name,np.allclose(manual,table[name].to_numpy()[rows],rtol=1e-6,atol=1e-8))
d=pd.read_csv(O/'tables/query_candidate_patient_program_scores.tsv',sep='\t');check('patient_summary_unique',not d.duplicated(['node_id','source','patient_key']).any());b=pd.read_csv(O/'tables/normal_IPF_background_comparison.tsv',sep='\t');check('percentile_bounds',b.median_reference_percentile.between(0,1).all());check('background_no_cross_dataset_p',not any(c in b for c in ['p','pvalue','fdr']))
t=pd.read_csv(O/'tables/IPF_within_study_patient_tests.tsv',sep='\t');valid=t.p.notna();check('valid_tests_have_replicates',((t.loc[valid,'ipf_patients']>=3)&(t.loc[valid,'healthy_patients']>=3)).all());check('valid_p_range',t.loc[valid,'p'].between(0,1).all());check('valid_fdr_range',t.loc[valid,'fdr'].between(0,1).all());check('no_p_for_insufficient_replicates',t.loc[~valid,'fdr'].isna().all())
if (O/'resources/GENE_EVIDENCE_COMPLETE.json').exists():
 marker=pd.read_csv(O/'tables/marker_patient_evidence.tsv',sep='\t');check('marker_detection_fraction_bounds',marker.fraction_detected.between(0,1).all());check('marker_expression_finite',np.isfinite(marker.mean_logCP10k).all())
 ipfgenes=['CCL2','COL1A1','CTHRC1','MMP19','SERPINE1','HIF1A'];mm=marker[(marker.kind=='query')&marker.gene.isin(ipfgenes)].groupby(['group','study','donor']).mean_logCP10k.mean().reset_index();mm['patient_key']=mm.study.astype(str)+'::'+mm.donor.astype(str);merged=mm.merge(d,left_on=['group','patient_key'],right_on=['node_id','patient_key'],validate='one_to_one');check('single_gene_reconstructs_IPF6',len(merged)==len(d) and np.allclose(merged.mean_logCP10k,merged.HLCA_IPF_reference))
if (O/'resources/FROZEN_RESUME_PARAMETER_CHECK.json').exists():
 frozen=json.loads((O/'resources/FROZEN_RESUME_PARAMETER_CHECK.json').read_text());check('sensitivity_reference_parameters_frozen',frozen['status']=='PASS' and frozen['frozen_max_change']==0 and frozen['encoder_gene_weight_max_change']==0)
 sensitivity=json.loads((O/'resources/MAPPING_FROZEN_EXTENDED_REVIEW.json').read_text());check('valid_frozen_sensitivity_complete',sensitivity['status']=='PASS_COMPUTATION')
pd.DataFrame(checks).to_csv(O/'tables/supplement_accuracy_checks.tsv',sep='\t',index=False);result=dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',checks=len(checks),passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']]);(O/'SUPPLEMENT_QA.json').write_text(json.dumps(result,indent=2),encoding='utf-8');assert result['status']=='PASS',result
files=[p for p in (O/'tables').glob('*.tsv') if not p.name.startswith('model_feature_map')]+[O/n for n in ['ANALYSIS_COMPLETE.json','SUPPLEMENT_QA.json','MAPPING_COMPLETE.json','REFERENCE_HTTP_COMPLETE.json']]+[O/'resources'/n for n in ['JOINT_PROGRAM_COVERAGE.json','MAPPING_REVIEW.json','HISTORICAL_MAPPING_REVIEW.json','MAPPING_FROZEN_EXTENDED_REVIEW.json','FROZEN_RESUME_PARAMETER_CHECK.json','GENE_EVIDENCE_COMPLETE.json','REFERENCE_FULL_DOWNLOAD_COMPLETE.json'] if (O/'resources'/n).exists()]
manifest={str(p.relative_to(O)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files};(O/'delivery/supplement_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
with zipfile.ZipFile(O/'delivery/supplement_review_bundle.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
 z.write(O/'delivery/supplement_manifest.json','supplement_manifest.json')
print(json.dumps(result),flush=True)
