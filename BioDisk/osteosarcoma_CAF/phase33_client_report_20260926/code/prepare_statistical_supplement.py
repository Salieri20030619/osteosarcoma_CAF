from pathlib import Path
import json,gc
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase33_client_report_20260926';M=P/'phase30_two_layer_classification_20260926';Q=P/'phase31_customer_completion_20260926'
OUT=O/'report_statistics';OUT.mkdir(exist_ok=True)
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
cm=read(Q/'resources/clinical_mapping_frozen.tsv').set_index('sample')
contrasts={r[3]:r[:3] for r in json.loads((Q/'clinical_internal/CONTRASTS_FROZEN.json').read_text())['contrasts']}
r30=read(M/'tables/independent_patient_comparison_roster.tsv')
r31=read(Q/'clinical_internal/independent_donor_roster.tsv.gz')
summaries=[]
for layer in ['supported','expanded']:
 a=ad.read_h5ad(M/f'objects/{layer}_integrated.h5ad');ob=read(M/f'tables/{layer}_classified_cells.tsv.gz').set_index('cell_id').loc[a.obs_names]
 for col in ['clinical_group','treatment','site']:ob[col]=ob['sample'].map(cm[col])
 ob['anatomic_site']=ob.site.str.lower().map(lambda s:'Lung' if 'lung' in s else ('Bone' if any(k in s for k in ['bone','femur','tibia','humerus','calcaneus','ulna']) else 'Other_or_unspecified'))
 my=read(M/f'tables/{layer}_myCAF_cell_labels.tsv.gz').set_index('cell_id');ob['my_subgroup']=my.cluster.reindex(ob.index)
 def zero_features(ids,features):
  if not len(ids):return np.zeros(len(features),bool)
  ix=a.obs_names.get_indexer(ids);assert (ix>=0).all()
  sums=np.asarray(a.layers['counts'][ix].sum(0)).ravel();fi=a.var_names.get_indexer(features)
  result=np.zeros(len(features),bool);good=fi>=0;result[good]=sums[fi[good]]==0
  for j,feature in enumerate(features):
   if not good[j] and feature in ob.columns:result[j]=bool(ob.loc[ids,feature].eq(0).all())
  return result
 def apply_mask(frame,indices,ids):
  fs=frame.loc[indices,'feature'].astype(str).to_numpy();mask=zero_features(ids,fs)
  ids_mask=indices[mask]
  if not len(ids_mask):return
  bad=frame.loc[ids_mask,'P'].notna()&~np.isclose(frame.loc[ids_mask,'P'].fillna(1),1)
  assert not bad.any(), 'All-zero feature unexpectedly had non-null/non-unit historical P'
  frame.loc[ids_mask,'historical_P_before_evaluability_annotation']=frame.loc[ids_mask,'P']
  frame.loc[ids_mask,'historical_FDR_before_evaluability_annotation']=frame.loc[ids_mask,'FDR']
  frame.loc[ids_mask,['P','FDR']]=np.nan
  frame.loc[ids_mask,'evaluation_status']='all_zero_in_compared_cells_not_tested'
 def init(d):
  d['evaluation_status']=np.where(d.P.notna(),'tested','no_valid_repeats_or_unavailable')
  d['historical_P_before_evaluability_annotation']=np.nan;d['historical_FDR_before_evaluability_annotation']=np.nan
  return d
 def save(d,name):
  d['nominal_P_lt_0_05']=d.P.lt(.05);d.to_csv(OUT/name,sep='\t',index=False)
  summaries.append(dict(file=name,rows=len(d),all_zero_rows=int(d.evaluation_status.eq('all_zero_in_compared_cells_not_tested').sum()),nominal_rows=int(d.nominal_P_lt_0_05.sum())))
 # Main-class complete gene families.
 for f in sorted((M/'tables').glob(f'{layer}_*_all_gene_statistics.tsv.gz')):
  d=init(read(f))
  for (pool,cls),g in d.groupby(['pool','class_id']):
   rr=r30[(r30.layer==layer)&(r30.pool==pool)&(r30.class_id==cls)];ids=[]
   for r in rr.itertuples():
    ids.extend(ob.index[ob.donor_record.eq(r.donor_record)&ob.site.isin(str(r.matched_sites).split(';'))])
   apply_mask(d,g.index,pd.Index(list(dict.fromkeys(ids))))
  save(d,'S02_'+f.name);del d;gc.collect()
 # Current internal and clinical expression families.
 for analysis in ['myCAF_internal','clinical_expression']:
  f=Q/f'clinical_internal/{layer}_{analysis}_all_feature_statistics.tsv.gz';d=init(read(f))
  for key,g in d.groupby(['study','target','contrast','paired'],dropna=False):
   study,target,contrast,paired=key
   rr=r31[(r31.layer==layer)&(r31.analysis==analysis)&(r31.study==study)&(r31.target==target)&(r31.contrast==contrast)&(r31.paired==paired)]
   if rr.empty:continue
   ids=[]
   if analysis=='myCAF_internal':
    cl=int(target.split('_')[-1]);scope=ob[ob.dataset.eq(study)&ob.my_subgroup.notna()&ob.donor_record.isin(rr.donor_record)]
    for (_,site),s in scope.groupby(['donor_record','site']):
     if s.my_subgroup.eq(cl).any() and s.my_subgroup.ne(cl).any():ids.extend(s.index)
   else:
    if target=='all_CAF_candidates':tm=np.ones(len(ob),bool)
    elif target.startswith('my_subgroup_'):tm=ob.my_subgroup.eq(int(target.split('_')[-1])).to_numpy()
    else:tm=ob.class_id.eq(target).to_numpy()
    field,lev1,lev2=contrasts[contrast]
    for cond,lev in [('first',lev1),('second',lev2)]:
     donors=rr.loc[rr.condition.eq(cond),'donor_record']
     ids.extend(ob.index[tm&ob.dataset.eq(study).to_numpy()&ob[field].eq(lev).to_numpy()&ob.donor_record.isin(donors).to_numpy()])
   apply_mask(d,g.index,pd.Index(list(dict.fromkeys(ids))))
  save(d,'S03_'+f.name);del d;gc.collect()
 del a;gc.collect()
pd.DataFrame(summaries).to_csv(OUT/'statistics_evaluability_manifest.tsv',sep='\t',index=False)
(OUT/'README.txt').write_text('Current report tables annotate features that are all zero in the actually compared cells as not tested (P/FDR NA). Historical numerical values, where present, are retained in explicitly named columns. Original source files are unchanged. Nonzero-feature P/FDR, directions, cohorts, cell populations and predefined hypothesis families are unchanged. The historical family computation used unit P placeholders for untestable features; this reporting annotation preserves the original predefined family and all other FDR values. Nominal screening uses original P < 0.05 in either direction; no FDR entry gate. No expression or cell inclusion threshold is added.\n')
(O/'delivery/STATISTICS_EVALUABILITY_QA.json').write_text(json.dumps({'files':len(summaries),'rows':sum(r['rows'] for r in summaries),'all_zero_rows_annotated':sum(r['all_zero_rows'] for r in summaries),'source_files_unchanged':True,'nonzero_P_and_FDR_unchanged':True,'no_candidate_cells_removed':True},indent=2))
print(json.dumps(summaries))
