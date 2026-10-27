from pathlib import Path
import json,hashlib,datetime,requests,gzip,shutil
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase8_original_plan_20260919';D=P/'phase7_review_states_20260917'
for s in ['tables','figures','bulk','validation','resources','objects','delivery']:(O/s).mkdir(exist_ok=True)
panels=json.loads((D/'analysis/gene_panels.json').read_text())
panels['Osteogenic_full_customer']=['RUNX2','SP7','DLX5','MSX2','ALPL','BGLAP','IBSP','SPP1','SPARC','COL1A1','SOX9','BMP2','BMP4','PHEX','DMP1','SOST']
panels['Calcium_regulation_full_customer']=['MGP','ENPP1','STC1','ANXA2','ANXA5','CALM1','CALM2','CACNA1C','SLC8A1','PTH1R','VDR','CASR','TRPM7','TNFRSF11B','TNFSF11']
panels['Myeloid_calcium_caution']=['S100A8','S100A9']
top=pd.read_csv(D/'integration_review/validation/cross_patient_top30_genes.tsv',sep='\t',dtype={'state':str})
for st in ['0','2']:
 # These were already ranked using discovery donors; no external outcome is inspected here.
 panels['Frozen_H'+st+'_top30']=top[top.state==st].symbol.tolist()
(O/'resources/programs_frozen.json').write_text(json.dumps(panels,indent=2))
surface=['PDGFRA','FAP','PDPN','THY1','DPP4','CD34','ITGA11','ALPL','GPC3','CD302','SCARA5','JAG1','GPRC5A','MCAM','PDGFRB','CD36']
checks=['PTPRC','PECAM1','EPCAM','S100A4','ACTA2','TAGLN','MYH11','RGS5','RUNX2','SP7','IBSP','SPP1']
(O/'resources/surface_candidates_frozen.json').write_text(json.dumps({'surface_candidates':surface,'controls':checks,'candidates_source':'discovery H0/H2 marker table plus predefined CAF surface panels; molecular localization and antibody evidence require separate review'},indent=2))
protocol={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'objective':'Finish original 10-day plan before evaluating CD36/MYH11/C7/FNDC1 reference-paper grouping',
 'discovery':'GSE152048 strict CAF 2047 cells from 6 donors; reviewed Harmony H0/H2 and existing myCAF subdivision fallback',
 'validation':'GSE318840 and GSE162454 analyzed separately; GSE250015 ancillary with overlap-unresolved flag',
 'patient_policy':'No cross-cohort pooling to inflate independent N; no unverified A/B1 or paired labeling',
 'bulk':'GSE21257 original metastasis-within-5-years endpoint; tissue expression, not CAF proportion',
 'bulk_primary_family':['Frozen_H0_top30','Frozen_H2_top30','Contractile','Inflammatory','Antigen_presentation','Matrix_remodeling','Osteogenic_full_customer','Calcium_regulation_full_customer'],
 'bulk_test':'patient-level standardized gene-set score; logistic OR per SD with Wald CI and P; BH over the 8 predefined program endpoints; no outcome-tuned signatures or cutpoints',
 'validation_expression':'frozen signatures/programs and marker expression; no new cluster discovery driven by validation outcomes; report donor-specific values and technique stratification',
 'flow_scope':'RNA-based candidate coverage and off-target risks plus protein localization and antibody evidence; channel-agnostic design without inventing instrument details',
 'data_dependent_enhancements':['A/B1 history and simultaneous paired samples','late ScPCA expression expansion','instrument-specific fluorochrome configuration'],
 'original_labels_changed':False}
(O/'PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2))
hashes=[]
for f in list((O/'sources').glob('*.docx'))+[D/'analysis/DISCOVERY_FROZEN.json',D/'integration_review/validation/cross_patient_top30_genes.tsv',D/'analysis/gene_panels.json',P/'geo_identity_20260919/STATUS.json']:
 hashes.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
pd.DataFrame(hashes).to_csv(O/'resources/input_state_hashes.tsv',sep='\t',index=False)
(O/'STATUS.json').write_text(json.dumps({'status':'ORIGINAL_PLAN_IN_PROGRESS','protocol':'PROTOCOL_FROZEN.json','next':'bulk data acceptance; external validation; marker and flow design; full report'},indent=2))
urls={
 'GSE21257_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21257/matrix/GSE21257_series_matrix.txt.gz',
 'GPL10295.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10295/annot/GPL10295.annot.gz',
 'GPL10295_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10295/soft/GPL10295_family.soft.gz'}
rows=[]
for name,url in urls.items():
 dest=O/'bulk'/name;error=''
 for attempt in range(2):
  try:
   with requests.get(url,stream=True,timeout=(20,120)) as r:
    r.raise_for_status()
    with open(dest.with_suffix(dest.suffix+'.part'),'wb') as f:
     for chunk in r.iter_content(1024*1024):f.write(chunk)
   temp=dest.with_suffix(dest.suffix+'.part')
   with gzip.open(temp,'rb') as f:
    while f.read(1024*1024):pass
   temp.replace(dest);error='';break
  except Exception as e:error=str(e);print('DOWNLOAD_RETRY',name,error,flush=True)
 rows.append({'file':name,'url':url,'complete':not bool(error),'error':error,'bytes':dest.stat().st_size if dest.exists() else 0})
 print('DOWNLOAD',rows[-1],flush=True)
pd.DataFrame(rows).to_csv(O/'bulk/download_manifest.tsv',sep='\t',index=False)
assert (O/'bulk/GSE21257_series_matrix.txt.gz').exists()
assert (O/'bulk/GPL10295.annot.gz').exists() or (O/'bulk/GPL10295_family.soft.gz').exists()
(O/'bulk/DOWNLOAD_COMPLETE').write_text('GEO_BULK_FILES_READABLE\n')
print('BULK_DOWNLOAD_COMPLETE',flush=True)
