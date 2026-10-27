from pathlib import Path
import os,json,re,datetime,hashlib,gzip
import pandas as pd,h5py
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase24_asset_inventory_20260921';T=O/'tables'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(d,n):d.to_csv(T/n,sep='\t',index=False)
files=[]
for root,dirs,names in os.walk(P):
 dirs[:]=[d for d in dirs if d not in ['Rlib','__pycache__','.git','phase24_asset_inventory_20260921']]
 for name in names:
  p=Path(root)/name
  try:s=p.stat();files.append({'path':str(p),'bytes':s.st_size,'mtime_ns':s.st_mtime_ns,'partial':name.endswith(('.part','.tmp'))})
  except OSError:pass
ff=pd.DataFrame(files);save(ff,'server_file_inventory.tsv');paths=[Path(r['path']) for r in files];pathset=set(map(str,paths))
old=read(P/'phase17_all_cohort_expanded_20260921/tables/sample_layer_counts.tsv').set_index('sample')
clinical=read(P/'phase19_external_validation_20260921/tables/52_sample_clinical_overview.tsv').set_index('sample')
meta=json.loads((O/'sources/scpca_source_metadata.json').read_text())
scp=[];headers=[]
def inspect_h5(p):
 try:
  with h5py.File(p,'r') as f:
   x=f['X'];shape=x.attrs.get('shape',x.shape if isinstance(x,h5py.Dataset) else None)
   return {'path':str(p),'readable_header':True,'cells':int(shape[0]),'genes':int(shape[1]),'X_encoding':str(x.attrs.get('encoding-type','dense')),'layers':';'.join(f['layers'].keys()) if 'layers' in f else ''}
 except Exception as e:return {'path':str(p),'readable_header':False,'error':str(e)}
for r in meta:
 sm=r['scpca_sample_id'];hits=[p for p in paths if sm in p.parts and p.name.endswith('.h5ad') and p.name.startswith('SCPCL')]
 hitmap={k:[p for p in hits if p.name.endswith('_'+k+'_rna.h5ad')] for k in ['filtered','processed','unfiltered']}
 row={k:r.get(k,'FIELD_ABSENT') for k in ['scpca_project_id','scpca_sample_id','scpca_library_id','participant_id','tissue_location','primary_or_metastasis','disease_timing','age','sex','sample_type','treatment_status','is_xenograft','is_cell_line','filtered_cell_count','processed_cells','technology','workflow_version']}
 row['source_metadata_date']='2026-09-10';row['non_PDX']=r['is_xenograft']=='False';row['past_52_sample_input']=sm in old.index
 for k,pp in hitmap.items():
  row[k+'_paths']=';'.join(map(str,pp));row[k+'_file_count']=len(pp)
  for p in pp:headers.append(dict(sample=sm,kind=k,**inspect_h5(p)))
 row['asset_state']='matrices_present' if hitmap['filtered'] and hitmap['processed'] else 'metadata_only_in_audited_locations'
 row['past_strict_CAF_cells']=int(old.loc[sm,'strict']) if sm in old.index else None
 row['past_supported_CAF_cells']=int(old.loc[sm,'supported']) if sm in old.index else None
 row['past_expanded_CAF_cells']=int(old.loc[sm,'expanded']) if sm in old.index else None
 row['formal_scope_decision']='PDX_separate_context_not_patient_CAF' if not row['non_PDX'] else 'include_in_source_universe_no_cell_count_gate'
 row['prior_omission_reason']='not_in_previous_12_sample_download_selection' if row['non_PDX'] and not row['past_52_sample_input'] else 'previously_included' if row['past_52_sample_input'] else 'PDX_separate_from_patient_tissue'
 scp.append(row)
sc=pd.DataFrame(scp);save(sc,'ScPCA_90_library_assets.tsv');save(pd.DataFrame(headers),'ScPCA_matrix_header_checks.tsv')
# Customer workbook contains one additional bulk-only record; do not inflate scRNA count.
books=json.loads((O/'sources/customer_workbooks.json').read_text());books_scp=next(v for k,v in books.items() if k.startswith('ScPCA_'))
allrows=next(rr for k,rr in books_scp['sheets'].items() if k.startswith('全部'))
wb=pd.DataFrame(allrows[1:],columns=allrows[0]);extra=wb[~wb['样本ID'].isin(sc.scpca_sample_id)];save(extra,'customer_records_without_scRNA_metadata.tsv')
pairrows=next(rr for k,rr in books_scp['sheets'].items() if k.startswith('肺转移'))
pairs=pd.DataFrame(pairrows[1:],columns=pairrows[0]).merge(sc,left_on='样本ID',right_on='scpca_sample_id',how='left',validate='one_to_one');save(pairs,'bone_lung_pair_assets.tsv')
legacy=[]
for rel in ['OS_CAF_phase2/code/manifest.tsv','phase5_nonpdx_qc_20260916/manifest.tsv','geo_preprocess_20260918/manifest.tsv']:
 for r in read(P/rel).to_dict('records'):
  sm=r['sample'];ds=old.loc[sm,'dataset'];entries=[]
  if 'matrix_dir' in r:
   d=Path(r['matrix_dir']);entries=[p for p in d.rglob('*') if p.is_file() and (p.name.endswith(('.h5ad','.mtx','.mtx.gz')))]
  else:
   name=r['source_url'].rsplit('/',1)[1];entries=[p for p in paths if p.name==name]
  recovered=P/'phase23_loose_first_20260921/inputs'/f'{sm}.h5ad';marker=recovered.with_suffix('.json')
  rec=json.loads(marker.read_text()) if marker.exists() else {}
  one=dict(sample=sm,dataset=ds,source_manifest=str(P/rel),source_matrix_paths=';'.join(map(str,entries)),source_path_count=len(entries),recovered_full_called_ready=recovered.exists() and marker.exists(),recovered_cells=rec.get('entry_cells'),historical_QC_cells=int(old.loc[sm,'QC_cells']),historical_strict_CAF=int(old.loc[sm,'strict']),historical_supported_CAF=int(old.loc[sm,'supported']),historical_expanded_CAF=int(old.loc[sm,'expanded']))
  if sm in clinical.index:one.update(clinical.loc[sm].to_dict())
  legacy.append(one)
le=pd.DataFrame(legacy);save(le,'GEO_40_sample_assets.tsv')
cs=json.loads((P/'phase7_review_states_20260917/candidate_sources_20260917/GEO_candidate_sample_manifest.json').read_text());omitted=pd.DataFrame([r for r in cs if r['dataset']=='GSE237070']);omitted['matrix_found_in_project']=False;omitted['role']=omitted['source'].map(lambda z:'excluded_from_A_B1_B2_by_customer' if 'Lymph' in z else 'primary_tumor_pending_matrix' if 'Primary' in z else 'adjacent_reference_pending_matrix');save(omitted,'GSE237070_source_records.tsv')
summary=[]
for ds,g in le.groupby('dataset'):
 summary.append(dict(dataset=ds,source_type='patient_tumor_scRNA',expected_sample_records=len(g),matrix_sample_records=int(g.source_path_count.gt(0).sum()),past_input_records=len(g),past_strict_nonzero_records=int(g.historical_strict_CAF.gt(0).sum()),past_supported_nonzero_records=int(g.historical_supported_CAF.gt(0).sum()),past_expanded_nonzero_records=int(g.historical_expanded_CAF.gt(0).sum()),gap='none_at_source_manifest_level; prior_CAF_gates_restrict_contribution'))
for ds,g in sc.groupby('scpca_project_id'):
 g=g[g.non_PDX];summary.append(dict(dataset=ds,source_type='non_PDX_patient_scRNA',expected_sample_records=len(g),matrix_sample_records=int(g.filtered_file_count.gt(0).sum()),past_input_records=int(g.past_52_sample_input.sum()),past_strict_nonzero_records=int(g.past_strict_CAF_cells.gt(0).sum()),past_supported_nonzero_records=int(g.past_supported_CAF_cells.gt(0).sum()),past_expanded_nonzero_records=int(g.past_expanded_CAF_cells.gt(0).sum()),gap=f'{int(g.filtered_file_count.eq(0).sum())} non_PDX samples have metadata but no matrix in audited project'))
for ds,typ,n,note in [('GSE237070','tumor_and_reference',2,'2 primary tumor sample records missing matrices; adjacent tissue separate; lymph nodes excluded per customer'),('GSE217792','bone_reference',None,'no matrix found; customer reference not implemented'),('HLCA','lung_reference',None,'no atlas object/classifier found; prior 5/3-gene panel is a substitute'),('Gong_2022','lung_signature',None,'six-gene human PTGS2 coexpression list exists; not full signature'),('Budhathoki_2026','source_identity_pending',None,'must resolve overlap/accession before counting as independent cohort')]:
 summary.append(dict(dataset=ds,source_type=typ,expected_sample_records=n,matrix_sample_records=0,past_input_records=0,gap=note))
su=pd.DataFrame(summary);save(su,'cohort_asset_summary.tsv')
# Capture outstanding inherited work without starting any new analysis.
procs=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  cwd=str((p/'cwd').resolve())
  if str(P/'phase23_loose_first_20260921') in cmd or cwd.startswith(str(P/'phase23_loose_first_20260921')):procs.append({'pid':int(p.name),'cwd':cwd,'command':cmd})
 except (OSError,PermissionError):pass
(O/'inherited_process_snapshot.json').write_text(json.dumps(procs,indent=2))
checks={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'phase':'ASSET_INVENTORY_ONLY','new_formal_analysis_started':False,'project_file_count':len(files),'ScPCA_metadata_libraries':len(sc),'ScPCA_non_PDX':int(sc.non_PDX.sum()),'ScPCA_PDX':int((~sc.non_PDX).sum()),'ScPCA_non_PDX_matrix_samples':int((sc.non_PDX & sc.filtered_file_count.gt(0)).sum()),'ScPCA_non_PDX_metadata_only':int((sc.non_PDX & sc.filtered_file_count.eq(0)).sum()),'legacy_GEO_source_records':len(le),'prior_analysis_records':len(old),'pair_sample_records':len(pairs),'pair_patient_ids':pairs.participant_id.nunique(),'ScPCA_header_readable':all(r['readable_header'] for r in headers),'formal_scope_ready':False,'reason':'previous 52-record scope omits 55 non-PDX ScPCA records and required reference/matrix assets; resolve acquisition and exact source identities before formal analysis','scan_scope':str(P),'patient_N_not_inferred_from_sample_counts':True}
assert len(sc)==90 and sc.scpca_sample_id.is_unique and len(le)==40 and len(pairs)==9 and not extra.empty
(O/'INVENTORY_QA.json').write_text(json.dumps(checks,indent=2));(O/'INVENTORY_COMPLETE').write_text('CENSUS_COMPLETE_NOT_FORMAL_ANALYSIS_COMPLETE\n')
print(json.dumps(checks,indent=2));print(su.to_string(index=False))
