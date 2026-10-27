from pathlib import Path
import json,hashlib,re,gzip,subprocess,shutil
import pandas as pd,numpy as np,h5py
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase53_de_novo_reannotation_20261002');B=R/'blind_review';U=R/'unblinding'
maps=json.loads((U/'identifier_maps.json').read_text())
keys=pd.read_csv(U/'cell_identifier_key.tsv.gz',sep='\t')
new=pd.read_csv(R/'tables/all_cell_reannotation_v1.tsv.gz',sep='\t',usecols=['cell_id','dataset','sample','patient_id','proposed_lineage','curation_id','CNV_evidence_current'])
roster=pd.read_csv(B/'anonymous_cell_cluster_roster.tsv.gz',sep='\t')
assert len(keys)==len(new)==len(roster)==876648
assert keys.cell_id.is_unique and keys.original_cell_id.is_unique and new.cell_id.is_unique and roster.cell_id.is_unique
assert set(keys.original_cell_id)==set(new.cell_id) and set(keys.cell_id)==set(roster.cell_id)
template=pd.read_csv(B/'reviewer_response_template.tsv',sep='\t',keep_default_na=False)
assert len(template)==307 and template.cluster_id.is_unique
assert (template[['reviewed_lineage','malignancy_assessment','supporting_genes','genomic_evidence','uncertainty','reviewer_notes']]=='').all().all()
assert set(roster.r10)==set(template.cluster_id)
rawmanifest=pd.read_csv(B/'RAW_COUNT_MANIFEST.tsv',sep='\t');raw_metadata_checks=[]
for row in rawmanifest.itertuples(index=False):
 f=B/row.relative_path;assert f.stat().st_size==row.bytes and hashlib.sha256(f.read_bytes()).hexdigest()==row.SHA256
 with h5py.File(f,'r') as h:
  assert set(h['obs'].keys())=={'cell_id','source_id','sample_id','patient_id','tissue','technology','library_counts','detected_features'}
  ids=[x.decode() if isinstance(x,bytes) else str(x) for x in h['obs/cell_id'][:]]
  assert len(ids)==row.cells and all(re.fullmatch(r'C_[a-f0-9]{24}',x) for x in ids)
 raw_metadata_checks.append({'file':row.relative_path,'cells':row.cells,'hash_verified':True,'allowed_obs_schema':True})
assert len(raw_metadata_checks)==107
# Verify R-to-HDF5-to-Python orientation and continuous values for every available cell.
expected=pd.read_csv(R/'tables/copykat_cell_continuous_metrics.tsv.gz',sep='\t',float_precision='round_trip').set_index('cell_id')
cnv_checks=[]
def strings(v):return [x.decode() if isinstance(x,bytes) else str(x) for x in np.asarray(v).ravel()]
for f in sorted((R/'objects/cnv').glob('*_copykat_profiles.h5')):
 with h5py.File(f,'r') as h:
  ids=strings(h['cell_id'][:]);bins=strings(h['bin_id'][:]);v=h['profiles'][:]
 assert v.shape==(len(ids),len(bins)) # Native R writer yields cell x bin in h5py.
 err=float(np.max(np.abs(np.abs(v).mean(1)-expected.loc[ids,'mean_abs_native_signal'].to_numpy())))
 assert err<1e-12
 cnv_checks.append({'file':f.name,'cells':len(ids),'bins':len(bins),'square_matrix':len(ids)==len(bins),'max_error':err})
assert len(cnv_checks)==84 and sum(x['cells'] for x in cnv_checks)==568754
assert not any(x['square_matrix'] for x in cnv_checks),'Initial plot orientation requires review for square matrices'
pd.DataFrame(cnv_checks).to_csv(R/'delivery/CNV_HDF5_INTEROPERABILITY_QA.tsv',sep='\t',index=False)
# Verify every primary expression matrix and evidence matrix has the same row/column identity.
for ds,code in maps['sources'].items():
 d=B/code;g=pd.read_csv(d/'genes.tsv',sep='\t',keep_default_na=False);c=pd.read_csv(d/'clusters.tsv',sep='\t')
 assert hashlib.sha256((d/'primary_cluster_counts.mtx.gz').read_bytes()).hexdigest()==hashlib.sha256((R/'objects'/ds/'primary_cluster_counts.mtx.gz').read_bytes()).hexdigest()
 with h5py.File(d/'all_cluster_marker_evidence.h5','r') as h:
  assert strings(h['feature_id'][:])==g.feature.tolist()
  assert strings(h['cluster_id'][:])==c.cluster_id.tolist()
  for field in ['mean_log_expression','rest_mean_log_expression','detection_fraction','patient_mean_log2_difference','positive_patients','negative_patients']:
   assert h[field].shape==(len(g),len(c))
  det=h['detection_fraction'][:];assert np.isfinite(det).all() and (det>=0).all() and (det<=1).all()
 assert c.full_library.notna().all() and (c.full_library>=c.analysis_library).all()
# No original source/sample/patient strings or assignment columns in blind text and PDF evidence.
tokens=list(maps['sources'])+list(maps['samples'])+list(maps['patients'])
pat=re.compile('|'.join(re.escape(x) for x in sorted(tokens,key=len,reverse=True)))
forbidden=['supported_class_id','expanded_class_id','proposed_lineage','curation_id','HPCA_label','BlueprintEncode_label','provider_consensus_celltype_annotation']
checked=[];pdf_tool=shutil.which('pdftotext')
for f in B.rglob('*'):
 if not f.is_file() or 'raw_counts' in f.parts:continue
 name=str(f.relative_to(B));assert not pat.search(name)
 text=None
 if f.suffix in ['.txt','.json','.tsv']:text=f.read_text()
 elif f.name.endswith('.tsv.gz'):
  with gzip.open(f,'rt') as z:text=z.read()
 elif f.suffix=='.pdf':
  assert pdf_tool,'PDF text validation tool unavailable'
  z=subprocess.run([pdf_tool,str(f),'-'],capture_output=True,text=True);assert z.returncode==0;text=z.stdout
 if text is not None:
  assert not pat.search(text),(name,pat.search(text).group() if pat.search(text) else '')
  header=text.splitlines()[0] if text.splitlines() else ''
  assert not any(x in header for x in forbidden),(name,header)
  checked.append(name)
# Genomic group tables and segment tables must remain pairwise connected after anonymization.
for f in (B/'regional_CNV').glob('*_regions.tsv.gz'):
 pair=f.with_name(f.name.replace('_regions.tsv.gz','_groups.tsv.gz'))
 reg=pd.read_csv(f,sep='\t');groups=pd.read_csv(pair,sep='\t')
 assert set(reg.group_id).issubset(set(groups.group_id))
 assert set(groups.cell_id).issubset(set(roster.cell_id))
q={'cells':876648,'primary_clusters':307,'raw_files_verified':107,'raw_count_arrays_prior_verified':True,'annotation_cell_identity_complete':True,'blind_template_empty':True,'blind_text_and_PDF_files_scanned':len(checked),'no_original_source_sample_patient_strings_in_blind_text':True,'no_assignment_columns_in_blind_text_headers':True,'CNV_profile_files':84,'CNV_profile_cells':568754,'CNV_max_R_to_Python_error':max(x['max_error'] for x in cnv_checks),'all_primary_matrix_identity_checks_passed':True,'native_region_group_links_passed':True}
(R/'delivery/FINAL_DATA_QA.json').write_text(json.dumps(q,indent=2));(R/'FINAL_DATA_QA_COMPLETE').write_text('SUCCESS');print(json.dumps(q),flush=True)
