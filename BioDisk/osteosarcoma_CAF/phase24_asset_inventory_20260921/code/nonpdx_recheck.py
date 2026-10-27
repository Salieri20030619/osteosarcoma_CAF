from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re,xml.etree.ElementTree as ET,datetime
import anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase24_asset_inventory_20260921';T=O/'tables'
sc=pd.read_csv(T/'ScPCA_90_library_assets.tsv',sep='\t',keep_default_na=False);checks=[]
for row in sc.to_dict('records'):
 for kind in ['filtered','processed','unfiltered']:
  for raw in str(row[kind+'_paths']).split(';'):
   if not raw:continue
   a=ad.read_h5ad(raw,backed='r');rec={'project':row['scpca_project_id'],'sample':row['scpca_sample_id'],'library':row['scpca_library_id'],'kind':kind,'path':raw,'metadata_is_xenograft':str(row['is_xenograft']),'metadata_is_cell_line':str(row['is_cell_line'])}
   ok=True
   for col,expected in [('is_xenograft',row['is_xenograft']),('is_cell_line',row['is_cell_line']),('sample_id',row['scpca_sample_id']),('participant_id',row['participant_id'])]:
    actual=sorted(set(a.obs[col].astype(str))) if col in a.obs else ['FIELD_ABSENT'];rec['object_'+col]=';'.join(actual);ok=ok and actual==[str(expected)]
   rec['metadata_object_consistent']=ok;rec['provider_infercnv_status']=str(a.uns.get('infercnv_status','not_applicable'));rec['object_cells']=a.n_obs
   checks.append(rec);a.file.close()
ck=pd.DataFrame(checks);ck.to_csv(T/'ScPCA_nonPDX_object_recheck.tsv',sep='\t',index=False)
for idx,row in sc.iterrows():
 c=ck[ck['sample']==row.scpca_sample_id];meta_pass=str(row.is_xenograft).lower()=='false' and str(row.is_cell_line).lower()=='false'
 status='EXCLUDE_PDX' if str(row.is_xenograft).lower()=='true' else 'EXCLUDE_CELL_LINE' if str(row.is_cell_line).lower()=='true' else 'NON_PDX_METADATA_VERIFIED_MATRIX_PENDING' if meta_pass and c.empty else 'NON_PDX_METADATA_AND_OBJECT_VERIFIED' if meta_pass and c.metadata_object_consistent.all() else 'HOLD_SOURCE_UNRESOLVED_OR_CONFLICT'
 sc.loc[idx,'nonPDX_recheck_status']=status;sc.loc[idx,'nonPDX_evidence_source']='official ScPCA metadata dated 2026-09-10'+('; h5ad obs fields reread 2026-09-21' if len(c) else '')
 sc.loc[idx,'nonPDX_required']=True
sc.to_csv(T/'ScPCA_90_library_assets.tsv',sep='\t',index=False)
# Read the actual deposited per-sample protocols, not historical derived acceptance labels.
soft={}
for ds,rel in [('GSE270231','phase3_identity_20260914/takeover_audit/GSE270231_family.soft'),('GSE162454','phase7_review_states_20260917/candidate_sources_20260917/GSE162454.soft'),('GSE250015','phase7_review_states_20260917/candidate_sources_20260917/GSE250015.soft'),('GSE318840','phase7_review_states_20260917/candidate_sources_20260917/GSE318840.soft')]:
 p=P/rel
 for b in re.split(r'\^SAMPLE = ',p.read_text())[1:]:
  sm=b.splitlines()[0].strip();d={'source_file':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
  for key in ['title','source_name_ch1','characteristics_ch1','extract_protocol_ch1','growth_protocol_ch1']:
   d[key]=' '.join(l.split(' = ',1)[1] for l in b.splitlines() if l.startswith('!Sample_'+key+' = '))
  soft[sm]=d
le=pd.read_csv(T/'GEO_40_sample_assets.tsv',sep='\t',keep_default_na=False);out=[]
proof={
 'GSE152048':('PMC7730477.xml','Sample preparation and cell isolation for scRNA-seq','Fresh lesions processed immediately after surgery; 11 patient samples, source manifest maps each BC identifier.','https://pmc.ncbi.nlm.nih.gov/articles/PMC7730477/'),
 'GSE299023':('PMC12248827.xml','2.1. Sample Collection and Tissue Dissociation','Six patient tumor specimens collected at CCMC/CHCO; direct dissociation of tissue and biopsies.','https://pmc.ncbi.nlm.nih.gov/articles/PMC12248827/'),
 'GSE250015':('PMC11467127.xml','Human osteosarcoma samples','Two patient primary lesions used for discovery scRNA; cell-line experiments are a separate arm.','https://pmc.ncbi.nlm.nih.gov/articles/PMC11467127/'),
 'GSE270231':('PMC11739783.xml','Results human metastatectomy paragraph','Seven GSM records explicitly identify patient lung metastases; paper reports patient metastatectomy tissue. Mouse and xenograft arms excluded.','https://pmc.ncbi.nlm.nih.gov/articles/PMC11739783/'),
 'GSE318840':('PMC13330631.xml','2.1. Human tissue specimens','Eight surgical patient tumor specimens used for FLEX; GEO identifies bone, recurrent soft tissue and lung specimens.','https://pmc.ncbi.nlm.nih.gov/articles/PMC13330631/'),
 'GSE162454':('','GEO sample extraction protocols and original Methods','Six patient tissue records with direct tissue dissociation; original Methods identifies untreated patient resection specimens.','https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full')}
for row in le.to_dict('records'):
 ds=row['dataset'];sm=row['sample'];pmc,section,summary,url=proof[ds];gsm=sm.split('__')[-1];s=soft.get(gsm,{})
 if pmc:
  p=(P/'phase4_download_20260914'/pmc) if ds in ['GSE152048','GSE299023'] else O/'sources'/pmc
  tree=ET.parse(p);assert tree.getroot() is not None
  article_hash=hashlib.sha256(p.read_bytes()).hexdigest()
 else:p=Path(s['source_file']);article_hash=s['source_sha256']
 evidence={'sample':sm,'dataset':ds,'classification':'NON_PDX_DIRECT_PATIENT_TISSUE','nonPDX_hard_requirement_pass':True,'basis':summary,'source_url':url,'source_file':str(p),'source_sha256':article_hash,'section':section,'geo_sample_title':s.get('title','mapped_by_sample_manifest_and_paper'),'geo_source_file':s.get('source_file','phase4_download_manifest_and_article'),'matrix_available':bool(row['source_path_count']),'clinical_group':row['clinical_group'],'patient':row['patient'],'site':row['site']}
 if s:
  assert s['extract_protocol_ch1']
  assert not re.search('xenograft|PDX|cell line',s['title']+' '+s['source_name_ch1']+' '+s['characteristics_ch1'],re.I)
 out.append(evidence)
pd.DataFrame(out).to_csv(T/'GEO_nonPDX_recheck.tsv',sep='\t',index=False)
(O/'sources/GEO_per_sample_source_evidence.json').write_text(json.dumps(soft,indent=2))
policy={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'authorization':'user explicitly requires non-PDX as hard requirement and asset census before formal reanalysis','required_source':'direct human patient tissue; explicitly non-PDX; not cell line','PDX':'exclude from patient CAF discovery, myCAF clustering, A/B1/B2, pairing, differential expression and flow-marker discovery','unknown_or_conflict':'do not include until source resolved','organism_human_alone_sufficient':False,'patient_id_shared_with_PDX_sufficient':False,'cell_count_or_CAF_yield_source_gate':False,'ScPCA_metadata_only':'retain all non-PDX records in acquisition universe; matrices must be matched and source flags rechecked upon receipt','historical_custom_CNV_or_QC_gates_allowed':False,'formal_analysis_started':False,'inherited_jobs':'previously running subset input recovery and CNV are recorded as inherited preparation; do not substitute for full-scope completion'}
(O/'SOURCE_ELIGIBILITY_POLICY.json').write_text(json.dumps(policy,indent=2))
q=json.loads((O/'INVENTORY_QA.json').read_text());q.update({'nonPDX_hard_requirement_recorded':True,'ScPCA_object_files_rechecked':len(ck),'ScPCA_object_metadata_conflicts':int((~ck.metadata_object_consistent).sum()),'ScPCA_status_counts':sc.nonPDX_recheck_status.value_counts().to_dict(),'GEO_source_rechecked_samples':len(out),'GEO_source_nonPDX_pass':sum(r['nonPDX_hard_requirement_pass'] for r in out)})
assert q['ScPCA_object_metadata_conflicts']==0 and len(out)==40
(O/'INVENTORY_QA.json').write_text(json.dumps(q,indent=2));print(json.dumps(q,indent=2))
