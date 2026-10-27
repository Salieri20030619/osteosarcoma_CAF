from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,datetime,sys,gc
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=R/'phase115_incremental_function_20261023'
for d in ['code','resources','tables','objects','logs','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
tracked=[]
def read(n):
 p=R/n;tracked.append(p);return pd.read_csv(p,sep='\t',low_memory=False)
def jread(n):
 p=R/n;tracked.append(p);return json.loads(p.read_text())
programs=jread('phase66_doublet_repair_DSC1_20261005/resources/PROGRAMS.json')['sets']
programs.update(jread('phase63_matrix_resistance_20261005/resources/POSTHOC_COLLAGEN_STRUCTURE_PROTOCOL.json')['panels'])
programs.update(jread('phase31_customer_completion_20260926/resources/customer_program_panels.json'))
programs.update({'F2C3_collagen_support':['COL3A1','COL5A1','COL5A2','PCOLCE','LOXL2']})
lineages=jread('phase27_identity_review_20260926/resources/lineage_panels.json')
baseline=read('phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz')
assert len(baseline)==8174 and baseline.cell_id.is_unique
query=read('phase93_expanded_admission_20261016/tables/expanded_registry_8980.tsv.gz');query=query[query.registry_origin.eq('NCC_addition_806')].copy()
cross=read('phase99_six_patient_admission_check_20261020/tables/author_corrected_crosswalk_pending_anatomy.tsv').set_index('sample')
nci=read('phase102_NCI_author_CAF_20261020/tables/author_CAF_cell_registry.tsv.gz')
collection=jread('phase93_expanded_admission_20261016/objects/COUNT_COLLECTION.json')
objects=[dict(path=collection['baseline'],source='baseline',ids=baseline.cell_id.tolist(),meta=baseline)]
for ent in collection['added_objects']:
 q=query[query['sample'].eq(ent['sample'])].copy();q['donor']=q['sample'].map(cross.author_corrected_patient_id)
 objects.append(dict(path=ent['path'],source='NCC',ids=q.cell_id.tolist(),meta=q))
for sm,q in nci.groupby('sample'):
 objects.append(dict(path=str(R/'phase95_NCI_expanded_reanalysis_20261016/inputs'/sm/'available_counts.h5ad'),source='NCI',ids=q.cell_id.tolist(),meta=q.copy()))
prespath=R/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz';tracked.append(prespath)
pres=sparse.load_npz(prespath).tocsr();si=read('phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv').set_index('sample').presence_row
shared=None;feature_maps=[]
for ent in objects:
 p=Path(ent['path']);tracked.append(p);a=ad.read_h5ad(p,backed='r')
 sy=a.var.symbol.astype(str);valid=~sy.duplicated(keep=False)&sy.ne('')&~sy.str.startswith('ENSG')
 if 'symbol_ambiguous' in a.var:valid &= ~a.var.symbol_ambiguous.astype(bool)
 if ent['source']=='baseline':
  rows=si.loc[baseline['sample'].unique()].to_numpy();valid &= np.asarray(pres[rows].sum(0)).ravel()==len(rows)
 mapping={s:i for i,s in enumerate(sy) if valid.iloc[i]}
 ent['symbol_map']=mapping;shared=set(mapping) if shared is None else shared&set(mapping)
 feature_maps.append(dict(path=str(p),source=ent['source'],unique_measured_symbols=len(mapping),features=a.n_vars,candidate_cells=len(ent['ids'])))
 a.file.close()
genes=sorted(shared);gindex={g:i for i,g in enumerate(genes)}
assert len(genes)>10000
dump('resources/PROGRAMS.json',programs);dump('resources/INPUT_OBJECTS.json',feature_maps)
coverage=[]
for k,gs in programs.items():coverage.append(dict(program=k,defined=len(gs),shared_measured=sum(g in shared for g in gs),coverage=sum(g in shared for g in gs)/len(gs),missing=';'.join(sorted(set(gs)-shared))))
pd.DataFrame(coverage).to_csv(O/'tables/program_coverage.tsv',sep='\t',index=False)
protocol=dict(created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Incremental functional evaluation of frozen F1-F5 and 33 states with NCC806 plus NCI1530; no figure production',baseline='8174 frozen cells; labels and raw counts unchanged',query='NCC author-corrected seven source-local patients; NCI17 accepted source-local patient groupings',no_cross_source_patient_pooling=True,no_Bone_Lung_inference=True,feature_alignment='Intersection of unique unambiguous gene symbols; baseline requires assayed presence in all baseline samples; no missing-feature zero filling',shared_genes=len(genes),expression='log1p(raw / full-library * 10000); full library before gene subset',mapping='Reference-only markers exclude union of all tested functional/lineage panel genes. Source-balanced donor-label templates and distinct-donor nearest-profile vote. Frozen labels remain unmodified.',validation='Deterministic source-stratified donor 60/20/20 train/calibration/test; thresholds from calibration only; heldout test uses identical cell-level acceptance. Final model uses training donors only so reported test and query rules match exactly. Additional full-reference nearest-template mapping is descriptive sensitivity only.',mapping_thresholds='Similarity and margin lower 5th percentile among calibration cells with correct dual-method label; per-label availability recorded; unknown retained. Thresholds not tuned to query counts.',statistics='Within-source patient/site label versus all other candidate cells; equal-weight site effects then donors; two-sided exact sign-flip for N<=16, MC65536 plus observed otherwise; mean log1pCP10k difference. N<2 or all-zero effects untested. Raw P<0.05 nominal; BH within source/level/arm and separate functional-program versus targeted-gene families.',arms='All candidates; singlet/non-model-positive; minimum5 target and rest cells per patient/site within singlet. Original-label baseline evaluated separately from exploratory transferred labels.',functional_scope='23 predefined programs including original11, collagen structure3, customer8, F2C3 support5gene. Tested score genes withheld from mapping; biological independence not guaranteed by feature exclusion.',interpretation='Transferred labels and resulting tests are exploratory conditional annotations, not validated subtype identities or independent mechanism validation. Unknown-label enrichment indicates reference mismatch, not proven new subtype.',no_new_cell_admission_threshold=True,original_results_modified=False)
dump('resources/PROTOCOL.json',protocol)
parts=[];metas=[]
truth=lambda x:x.astype(str).str.lower().isin(['true','1','1.0'])
for ent in objects:
 a=ad.read_h5ad(ent['path'],backed='r');ids=a.obs_names.get_indexer(ent['ids']);assert (ids>=0).all()
 # Read only the requested rows, then retain shared genes. Backed CSR slicing is row-wise.
 xx=a.X[ids,:].tocsr();den=np.asarray(xx.sum(1)).ravel();assert (den>0).all()
 m=ent['meta'].copy().reset_index(drop=True)
 if ent['source']=='baseline':
  assert np.allclose(den,m.full_library_counts)
  out=pd.DataFrame(dict(cell_id=m.cell_id,source=m.dataset,donor=m.donor_record,sample=m['sample'],site=m.site.fillna('Unknown'),cohort='baseline',class_label=m.expanded_class_id,state_label=m.subgroup_id,doublet=truth(m.doublet_positive_any),lineage=m.review_lineage_code,conservative_fibro=truth(m.conservative_fibro)))
 elif ent['source']=='NCC':
  assert np.allclose(den,m.n_counts)
  out=pd.DataFrame(dict(cell_id=m.cell_id,source='NCC',donor=m.donor,sample=m['sample'],site='Unknown',cohort='added',class_label='',state_label='',doublet=m.full_scDblFinder_class.eq('doublet'),lineage=m.lineage_phase3.fillna('Unknown'),conservative_fibro=truth(m.ambient_supported_reference_validated_candidate)))
 else:
  assert np.allclose(den,a.obs.iloc[ids].n_counts)
  out=pd.DataFrame(dict(cell_id=m.cell_id,source='NCI',donor=m.patient_id,sample=m['sample'],site='Lung',cohort='added',class_label='',state_label='',doublet=m.full_scDblFinder_class.eq('doublet'),lineage='Author_CAF',conservative_fibro=False))
 out['library_counts']=den;out['donor_key']=out.source.astype(str)+'::'+out.donor.astype(str)
 ix=[ent['symbol_map'][g] for g in genes];raw=xx[:,ix].astype(np.float64);norm=raw.multiply((10000/den)[:,None]).tocsr();norm.data=np.log1p(norm.data)
 parts.append(norm);metas.append(out);a.file.close();print('EXTRACTED',ent['source'],len(out),flush=True)
matrix=sparse.vstack(parts,format='csr');meta=pd.concat(metas,ignore_index=True);assert len(meta)==10510 and meta.cell_id.is_unique
assert meta[meta.source.eq('NCC')].donor.nunique()==7 and meta[meta.source.eq('NCI')].donor.nunique()==17
sparse.save_npz(O/'objects/log1p_cp10k_shared.npz',matrix);meta.to_csv(O/'tables/cell_registry.tsv.gz',sep='\t',index=False);pd.DataFrame({'symbol':genes}).to_csv(O/'tables/shared_genes.tsv',sep='\t',index=False)
scores=meta[['cell_id','source','donor','site']].copy()
for k,gs in programs.items():
 ix=[gindex[g] for g in gs if g in gindex];scores[k]=np.asarray(matrix[:,ix].mean(1)).ravel() if ix else np.nan
for k,gs in lineages.items():
 ix=[gindex[g] for g in gs if g in gindex];scores['lineage__'+k]=np.asarray(matrix[:,ix].mean(1)).ravel() if ix else np.nan
target=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','PCOLCE','LOX','LOXL2','POSTN','CD248','ACTA2','TAGLN','MYH11','FN1','COL4A1','COL4A2']
for g in target:scores['gene__'+g]=matrix[:,gindex[g]].toarray().ravel() if g in gindex else np.nan
scores.to_csv(O/'tables/cell_function_scores.tsv.gz',sep='\t',index=False)
excluded=set().union(*map(set,programs.values()),*map(set,lineages.values()),set(target))
dump('resources/MAPPING_EXCLUDED_GENES.json',sorted(excluded))
hashes=[dict(path=str(p),sha256=sha(p)) for p in dict.fromkeys(tracked)];dump('resources/INPUT_HASHES.json',hashes)
dump('EXTRACTION_COMPLETE.json',dict(status='PASS',cells=len(meta),genes=len(genes),source_counts=meta.groupby('source').size().to_dict(),source_donors=meta.groupby('source').donor.nunique().to_dict(),programs=len(programs),all_ids_unique=True,full_library_denominators_checked=True,patient_correction_applied=True))
print('EXTRACTION_COMPLETE',flush=True)
