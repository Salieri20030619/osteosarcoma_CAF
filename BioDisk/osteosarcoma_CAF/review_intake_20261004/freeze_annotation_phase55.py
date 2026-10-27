from pathlib import Path
import json,hashlib,datetime,shutil
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase53_de_novo_reannotation_20261002';R=P/'phase54_raw_review_20261004';V=P/'review_intake_20261004/received_review';O=P/'phase55_annotation_freeze_20261004'
for sub in ['code','tables','resources','delivery','logs']:(O/sub).mkdir(parents=True,exist_ok=True)
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
assert (R/'REVIEW_COMPLETE').exists() and (R/'finalize.exit_code').read_text().strip()=='0'
comp=read(V/'unblinded_comparison/cluster_annotation_comparison.tsv');dec=read(R/'tables/CLUSTER_REVIEW_DECISIONS.tsv');assert comp.curation_id.is_unique
fib=comp[comp.lineage_code.eq('Fibro')];agreed=fib[fib.proposed_lineage.eq('Fibroblast_like')];added=fib[~fib.proposed_lineage.eq('Fibroblast_like')];flag=fib[fib.malignancy_assessment.str.startswith('CNV支持恶性候选')];patientflag=fib[fib.malignancy_assessment.str.startswith('部分患者')];mixed=comp[comp.lineage_code.str.startswith('Mix')]
assert len(agreed)==11 and agreed.cells.sum()==45721 and len(added)==4 and added.cells.sum()==45291 and len(flag)==3 and len(mixed)==13
scopes={'agreed11':agreed.cluster_id.tolist(),'agreed11_no_flag':agreed.loc[~agreed.cluster_id.isin(flag.cluster_id),'cluster_id'].tolist(),'all15':fib.cluster_id.tolist(),'all15_no_flag':fib.loc[~fib.cluster_id.isin(flag.cluster_id),'cluster_id'].tolist()}
protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope_definitions':scopes,'main_descriptive_scope':'agreed11','CNV_cluster_flags_from_frozen_review':flag.cluster_id.tolist(),'patient_specific_CNV_flags_not_removed':patientflag.cluster_id.tolist(),'fixed_background':'All input cells outside all15 and outside the 13 frozen mixed-program clusters; identical background across scope sensitivities.','secondary_comparator':'Frozen blind Osteo clusters; same patient and same site.','no_new_numeric_cell_gate':True,'no_cells_deleted':True,'no_confirmed_CAF_label_created':True,'original_and_reviewed_labels_parallel':True,'scope_sensitivity_is_not_independent_validation':True,'old_F3_F4_not_reused_as_CAF_identity':True}
dump(O/'resources/ANNOTATION_AND_SCOPE_PROTOCOL_FROZEN.json',protocol)
cols=['cell_id','sample','dataset','patient_id','site','technology','full_library','analysis_library','r05','r10','r20','curation_id','proposed_lineage','state','status','copykat_call','infercnv_call','infercnv_source','CNV_evidence_current','identity_interpretation']
a=read(S/'tables/all_cell_reannotation_v1.tsv.gz',usecols=cols).rename(columns={'cell_id':'original_cell_id','patient_id':'original_donor_record','r05':'original_r05','r10':'original_r10','r20':'original_r20','state':'original_state','status':'original_status'})
k=read(S/'unblinding/cell_identifier_key.tsv.gz')[['original_cell_id','cell_id','sample_id','source_id','patient_id']];a=a.merge(k,on='original_cell_id',validate='one_to_one');assert len(a)==876648 and a.cell_id.is_unique
cr=comp[['curation_id','cluster_id','reviewed_lineage','lineage_code','lineage_confidence','malignancy_assessment','comparison_category']].rename(columns={'reviewed_lineage':'review_expression_lineage','lineage_code':'review_lineage_code'})
a=a.merge(cr,on='curation_id',validate='many_to_one');assert a.cluster_id.notna().all()
a=a.merge(dec[['cluster_id','phase54_review_decision']],on='cluster_id',how='left',validate='many_to_one');a.phase54_review_decision=a.phase54_review_decision.fillna('not_in_phase54_deep_review')
a['candidate_tier']=np.select([a.cluster_id.isin(agreed.cluster_id),a.cluster_id.isin(added.cluster_id)],['agreed_fibroblast_like','additional_fibroblast_like'],default='outside_fibroblast_candidate')
a['mixed_program_flag']=a.cluster_id.isin(mixed.cluster_id);a['CNV_review_cluster_flag']=a.cluster_id.isin(flag.cluster_id);a['CNV_patient_specific_flag']=a.cluster_id.isin(patientflag.cluster_id)
a['fixed_background']=~a.cluster_id.isin(fib.cluster_id)&~a.mixed_program_flag;a['osteogenic_reference']=a.review_lineage_code.eq('Osteo')
for scope,cls in scopes.items():a['scope_'+scope]=a.cluster_id.isin(cls)
a['confirmed_CAF_identity']=False
reviewcells=read(R/'tables/target_cell_review_overlay.tsv.gz',usecols=['cell_id','doublet_assessed_any','doublet_positive_any','decontx_contamination']);a=a.merge(reviewcells,on='cell_id',how='left',validate='one_to_one')
a.to_csv(O/'tables/all_cell_annotation_v2.tsv.gz',sep='\t',index=False)
cluster=comp.merge(dec[['cluster_id','phase54_review_decision','interpretation_ZH']],on='cluster_id',how='left',validate='one_to_one');cluster.to_csv(O/'tables/cluster_annotation_v2.tsv',sep='\t',index=False)
counts=[];pat=[]
for sm,g in a.groupby('sample',sort=True):
 for scope in scopes:
  n=int(g['scope_'+scope].sum());counts.append({'dataset':g.dataset.iloc[0],'source_id':g.source_id.iloc[0],'sample':sm,'sample_id':g.sample_id.iloc[0],'patient_id':g.patient_id.iloc[0],'donor_record':g.original_donor_record.iloc[0],'site':g.site.iloc[0],'scope':scope,'candidate_cells':n,'full_called_cells':len(g),'fraction_of_full_called':n/len(g),'CNV_flagged_candidate_cells':int((g['scope_'+scope]&g.CNV_review_cluster_flag).sum())})
for (ds,donor,site),g in a.groupby(['dataset','original_donor_record','site'],sort=True):
 for scope in scopes:pat.append({'dataset':ds,'donor_record':donor,'patient_id':g.patient_id.iloc[0],'site':site,'scope':scope,'samples':g['sample'].nunique(),'candidate_cells':int(g['scope_'+scope].sum()),'full_called_cells':len(g),'fraction_of_full_called':g['scope_'+scope].mean()})
pd.DataFrame(counts).to_csv(O/'tables/sample_candidate_distribution.tsv',sep='\t',index=False);pd.DataFrame(pat).to_csv(O/'tables/patient_site_candidate_distribution.tsv',sep='\t',index=False)
pats=pd.DataFrame(pat).groupby(['dataset','donor_record','patient_id','scope'],as_index=False).agg(candidate_cells=('candidate_cells','sum'),full_called_cells=('full_called_cells','sum'),sites=('site','nunique'),site_equal_mean_fraction=('fraction_of_full_called','mean'));pats['fraction_pooled_cells']=pats.candidate_cells/pats.full_called_cells;pats.to_csv(O/'tables/patient_candidate_distribution.tsv',sep='\t',index=False)
ss=pd.DataFrame(counts).groupby(['dataset','scope'],as_index=False).agg(candidate_cells=('candidate_cells','sum'),full_called_cells=('full_called_cells','sum'),samples_total=('sample','nunique'),samples_with_candidates=('candidate_cells',lambda x:int((x>0).sum())));pp=pats.groupby(['dataset','scope']).agg(patient_records=('donor_record','nunique'),patients_with_candidates=('candidate_cells',lambda x:int((x>0).sum()))).reset_index();ss=ss.merge(pp,on=['dataset','scope'],validate='one_to_one');ss.to_csv(O/'tables/source_candidate_distribution.tsv',sep='\t',index=False)
sourcefiles=[S/'tables/all_cell_reannotation_v1.tsv.gz',S/'unblinding/cell_identifier_key.tsv.gz',V/'blind_submission/BLIND_FREEZE.json',V/'unblinded_comparison/cluster_annotation_comparison.tsv',R/'tables/CLUSTER_REVIEW_DECISIONS.tsv']
dump(O/'resources/INPUT_SHA256.json',[{'path':str(p),'SHA256':sha(p)} for p in sourcefiles])
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':'annotation_v2_frozen','cells':len(a),'clusters':a.cluster_id.nunique(),'samples':a['sample'].nunique(),'scope_counts':{scope:int(a['scope_'+scope].sum()) for scope in scopes},'CNV_flagged_clusters':flag.cluster_id.tolist(),'all_cell_ids_unique':True,'input_and_output_cell_sets_identical':set(a.cell_id)==set(k.cell_id),'scope_order_nesting_verified':bool((a.scope_agreed11_no_flag<=a.scope_agreed11).all() and (a.scope_agreed11<=a.scope_all15).all()),'annotation_SHA256':sha(O/'tables/all_cell_annotation_v2.tsv.gz'),'original_formal_tables_replaced':False,'new_confirmed_CAF_labels':0}
assert qa['clusters']==307 and qa['samples']==107 and qa['scope_counts']=={'agreed11':45721,'agreed11_no_flag':41193,'all15':91012,'all15_no_flag':86484}
dump(O/'delivery/FINAL_ACCEPTANCE.json',qa);(O/'ANNOTATION_FREEZE_COMPLETE').write_text('SUCCESS\n');print(json.dumps(qa,ensure_ascii=True))