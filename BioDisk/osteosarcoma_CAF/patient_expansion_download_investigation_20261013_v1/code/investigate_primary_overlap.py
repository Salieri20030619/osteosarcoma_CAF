
from pathlib import Path
import json,re,h5py,numpy as np,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1'
m=pd.read_csv(O/'tables/OsteoCAR/patient_prim_author_metadata.tsv.gz',sep='\t',low_memory=False);clinical=pd.read_csv(P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv',sep='\t');source=pd.read_csv(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',sep='\t');cached={};bc=re.compile(r'([ACGT]{16})')
def barcodes(s):return [bc.search(v).group(1) if bc.search(v) else None for v in map(str,s)]
for rr in source[source.dataset.isin(['GSE152048','GSE162454','SCPCP000023'])].itertuples():
 f=P/'phase26_client_reanalysis_20260922/inputs'/rr.sample/'full_called_counts.h5ad'
 with h5py.File(f,'r') as h:
  obs=h['obs'];key=obs.attrs['_index'];v=obs[key].asstr()[:];ids=barcodes(v)
 cached[rr.sample]=set(b for b in ids if b is not None)
rows=[];mapping=[]
for (src,sm,model,site),g in m.groupby(['data_source','sample_name','model','location']):
 ids=barcodes(g.cell_id);assert all(v is not None for v in ids);g=g.copy();g['barcode16']=ids;bb=set(ids)
 if src=='SJ':
  donor=re.match(r'(SJOS\d+)',model).group(1);c=clinical[(clinical.dataset=='SCPCP000023')&(clinical.patient_id.astype(str)==donor)&(clinical.strict_site_group=='Bone')]
 else:
  donor=None;c=clinical[clinical.dataset.isin(['GSE152048','GSE162454'])&(clinical.strict_site_group=='Bone')]
 candidates=[]
 for r in c.itertuples():
  overlap=len(bb&cached[r.sample]);rec=dict(author_sample=sm,author_model=model,author_source=src,author_cells=len(g),project_sample=r.sample,project_dataset=r.dataset,project_donor=r.donor_record,shared_barcodes=overlap,query_barcode_fraction=overlap/len(bb),project_barcodes=len(cached[r.sample]));candidates.append(rec)
 candidates=sorted(candidates,key=lambda z:z['shared_barcodes'],reverse=True);rows.extend(candidates)
 if candidates:
  best=candidates[0];best=best.copy();best['second_best_overlap']=candidates[1]['shared_barcodes'] if len(candidates)>1 else 0;best['identity_interpretation']='known SJ patient ID with same-site barcode comparison' if src=='SJ' else 'barcode overlap evidence; exact count fingerprint checked separately';mapping.append(best)
pd.DataFrame(rows).to_csv(O/'tables/OsteoCAR/patient_prim_all_overlap_comparisons.tsv',sep='\t',index=False);out=pd.DataFrame(mapping);out.to_csv(O/'tables/OsteoCAR/patient_prim_best_overlap_candidates.tsv',sep='\t',index=False);print(out[['author_sample','author_model','project_sample','shared_barcodes','author_cells','query_barcode_fraction','second_best_overlap']].to_string(index=False))
(O/'resources/OSTEOCAR_PRIM_OVERLAP_COMPLETE.json').write_text(json.dumps({'status':'completed','author_samples':m.sample_name.nunique(),'comparison_rows':len(rows),'note':'No new biological patients claimed from barcode matching alone'},indent=2))
