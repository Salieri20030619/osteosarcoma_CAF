
from pathlib import Path
import pandas as pd,numpy as np,json,gzip,re
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1';idx=pd.read_csv(O/'tables/ORIGINAL_SOURCE_BARCODE_INDEX.tsv',sep='\t');clinical=pd.read_csv(P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv',sep='\t');cache={}
for r in idx.itertuples():
 with gzip.open(r.barcode_list,'rt') as f:cache[r.sample]=set(f.read().splitlines())
parts=[]
for kind in ['patient_prim','patient_mets']:
 m=pd.read_csv(O/'tables/OsteoCAR'/f'{kind}_author_metadata.tsv.gz',sep='\t',low_memory=False);m['object']=kind;parts.append(m)
m=pd.concat(parts,ignore_index=True)
for col in ['object','data_source','sample_name','model','location']:m[col]=m[col].fillna('not_reported').astype(str)
m['data_source_original']=m.data_source
m['data_source']=m.data_source.replace({'NCI_POB':'NCI-POB'})
m['author_CAF']=m.Ann_Level2.eq('CAFs')&m.Ann_Level0.eq('Host');m['strict_author_site']=np.where(m.location.str.lower().eq('lung'),'Lung',np.where(m.location.str.lower().isin(['tibia','femur','humerus','bone','pelvis','sacrum','vertebra']),'Bone','Other_or_unclear'))
roster=[];overlaps=[]
for key,g in m.groupby(['object','data_source','sample_name','model','location','strict_author_site']):
 kind,ds,sm,model,location,site=key;bb=set(g.cell_id.str.extract(r'([ACGT]{16})')[0].dropna());comparison=[]
 for r in idx.itertuples():
  n=len(bb&cache[r.sample]);comparison.append(dict(object=kind,author_source=ds,author_sample=sm,author_model=model,author_site=site,query_cells=len(g),query_unique_barcodes=len(bb),project_sample=r.sample,project_dataset=r.dataset,project_donor=r.donor_record,project_site=r.site,shared_barcodes=n,query_overlap_fraction=n/len(bb) if bb else np.nan))
 comparison.sort(key=lambda r:r['shared_barcodes'],reverse=True);overlaps.extend(comparison);best=comparison[0];second=comparison[1]
 patient_key=None;identity_status='unresolved_author_model_is_not_verified_patient_id';source_roles='candidate_external_source'
 if ds=='SJ' and re.match(r'SJOS\d+',model):
  patient_key=re.match(r'SJOS\d+',model).group(0);identity_status='verified_existing_patient_ID' if clinical.patient_id.astype(str).eq(patient_key).any() else 'author_patient_ID_not_in_existing_registry';source_roles='SCPCP000023_related'
 elif ds=='GEO':
  patient_key=best['project_donor'];identity_status='known_GEO_label_and_barcode_correspondence_to_existing_source';source_roles='existing_GEO_reanalysis'
 elif ds=='NCH':
  patient_key=model;identity_status='known_candidate_overlap_with_GSE270231;barcode_evidence_retained';source_roles='GSE270231_related'
 hypothesis=re.search(r'_(\d+)',model).group(1) if ds=='NCI-POB' and re.search(r'_(\d+)',model) else None
 roster.append(dict(object=kind,author_source=ds,author_sample=sm,author_model=model,author_location=location,site=site,cells=len(g),author_CAF_cells=int(g.author_CAF.sum()),author_fibrogenic_tumor_cells=int((g.Ann_Level2.eq('Fibrogenic')&g.Ann_Level0.eq('Tumor')).sum()),patient_key=patient_key,identity_status=identity_status,source_role=source_roles,unverified_NCI_name_based_patient_key=hypothesis,best_existing_sample=best['project_sample'],best_existing_dataset=best['project_dataset'],best_existing_donor=best['project_donor'],best_existing_site=best['project_site'],shared_barcodes=best['shared_barcodes'],query_overlap_fraction=best['query_overlap_fraction'],second_best_shared_barcodes=second['shared_barcodes']))
r=pd.DataFrame(roster);assert int(r.cells.sum())==len(m);assert int(r.author_CAF_cells.sum())==int(m.author_CAF.sum());r.to_csv(O/'tables/OsteoCAR/ALL_HUMAN_SAMPLE_CAF_ROSTER.tsv',sep='\t',index=False);pd.DataFrame(overlaps).to_csv(O/'tables/OsteoCAR/ALL_HUMAN_ORIGINAL_BARCODE_OVERLAPS.tsv.gz',sep='\t',index=False)
source=m.groupby(['object','data_source','strict_author_site']).agg(cells=('cell_id','size'),sample_records=('sample_name','nunique'),author_CAF_cells=('author_CAF','sum')).reset_index();source.to_csv(O/'tables/OsteoCAR/OBJECT_SOURCE_CAF_CAPACITY.tsv',sep='\t',index=False)
labels=m.groupby(['object','data_source','Ann_Level0','Ann_Level1','Ann_Level2','Ann_Level3']).size().reset_index(name='cells');labels.to_csv(O/'tables/OsteoCAR/AUTHOR_ANNOTATION_HIERARCHY.tsv',sep='\t',index=False)
caf=m[m.author_CAF].copy();caf.to_csv(O/'tables/OsteoCAR/ALL_AUTHOR_CAF_METADATA.tsv.gz',sep='\t',index=False)
sj=r[(r.author_source=='SJ')&r.site.isin(['Bone','Lung'])].copy();capacity=[]
for mincells in [1,5]:
 b=set(sj[(sj.site=='Bone')&(sj.author_CAF_cells>=mincells)].patient_key);l=set(sj[(sj.site=='Lung')&(sj.author_CAF_cells>=mincells)].patient_key);pair=sorted(b&l);capacity.append(dict(minimum_author_CAF_cells_per_sample=mincells,Bone_patients=len(b),Lung_patients=len(l),paired_patients=len(pair),independent_Bone=len(b-set(pair)),independent_Lung=len(l-set(pair)),paired_IDs=';'.join(pair),note='descriptive capacity under author annotation; original analysis not replaced'))
pd.DataFrame(capacity).to_csv(O/'tables/OsteoCAR/SJ_AUTHOR_CAF_PATIENT_CAPACITY.tsv',sep='\t',index=False)

conflicts=json.loads((O/'resources/CONFIRMED_SAMPLE_SITE_CONFLICTS.json').read_text());blocked_samples={c['author_sample'] for c in conflicts};safe=sj[~sj.author_sample.isin(blocked_samples)];safe_capacity=[]
for mincells in [1,5]:
 b=set(safe[(safe.site=='Bone')&(safe.author_CAF_cells>=mincells)].patient_key);l=set(safe[(safe.site=='Lung')&(safe.author_CAF_cells>=mincells)].patient_key);pair=sorted(b&l);safe_capacity.append(dict(minimum_author_CAF_cells_per_sample=mincells,Bone_patients=len(b),Lung_patients=len(l),paired_patients=len(pair),independent_Bone=len(b-set(pair)),independent_Lung=len(l-set(pair)),paired_IDs=';'.join(pair),excluded_author_samples=';'.join(sorted(blocked_samples)),note='Capacity after excluding unresolved provider-versus-author anatomical-site conflicts; not new biological inference'))
pd.DataFrame(safe_capacity).to_csv(O/'tables/OsteoCAR/SJ_SITE_CONFLICT_EXCLUDED_CAPACITY.tsv',sep='\t',index=False)

summary={'status':'complete','total_cells':len(m),'author_CAF_cells':int(m.author_CAF.sum()),'observed_sample_records':m.sample_name.nunique(),'NCI_lung_records':int(((r.author_source=='NCI-POB')&(r.site=='Lung')).sum()),'NCI_lung_records_with_author_CAF':int(((r.author_source=='NCI-POB')&(r.site=='Lung')&r.author_CAF_cells.gt(0)).sum()),'NCI_independent_new_patients':None,'SJ_capacity':capacity,'SJ_site_conflict_excluded_capacity':safe_capacity,'no_pooled_cross_source_inference':True,'baseline_unchanged':True}
(O/'resources/OSTEOCAR_PATIENT_INVESTIGATION_COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print('SOURCE_CAPACITY');print(source.to_string(index=False));print('NCI_RECORDS');print(r[r.author_source=='NCI-POB'].to_string(index=False));print('SJ_CAPACITY');print(pd.DataFrame(capacity).to_string(index=False))
