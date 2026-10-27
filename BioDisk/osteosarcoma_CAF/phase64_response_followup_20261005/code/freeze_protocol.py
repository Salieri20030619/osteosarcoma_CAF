from pathlib import Path
import csv,gzip,json,datetime,re
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase64_response_followup_20261005')
protocol={'frozen_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary':'fibrillar_I_III_V, poor minus good pathological response; each cohort separately; no FDR required for one prespecified primary endpoint per cohort','secondary':'other 13 fixed phase63 programs; BH within cohort/model; all 14 program BH also provided','normalization':'Use submitted lumi VST normalized expression unchanged. Multiple uniquely mapped gene probes median. Per-gene z scores across all submitted samples, equally averaged; coverage>=80% and>=3 genes. No response-directed gene selection. Gene tests use all measured fixed-program genes; BH separate family.','test':'Unadjusted two-sided Welch test (poor minus good); Hedges g and raw effect; primary score also 99999 label permutations. Adjusted OLS score~poor+age_years+male (+metastases_at_diagnosis only GSE21257), HC3 t inference. Age/sex-complete cases; no imputation.','response':'GSE39055>=90% necrosis good; <90% poor; >80% ambiguous excluded primary and both extreme assignments sensitivity. GSE21257 Huvos1/2 poor,3/4 good,Unknown excluded.','limits':'Independent biopsies as units. Tissue bulk associations do not identify CAF-specific source, drug delivery or causality. GPL assay coverage documented. GSE21257 selected by future metastasis; no future-metastasis covariate.','frozen_program_source':'phase63 PROGRAMS.json plus POSTHOC_COLLAGEN_STRUCTURE_PROTOCOL.json'}
p=B/'resources'/'CLINICAL_PROTOCOL_FROZEN.json'
if not p.exists():p.write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
for g in ['GSE39055','GSE21257']:
 rows=[]
 with gzip.open(B/'sources'/f'{g}_series_matrix.txt.gz','rt') as f:
  for l in f:
   if l.startswith('!series_matrix_table_begin'):break
   if l.startswith('!Sample'):rows.append(next(csv.reader([l],delimiter='\t')))
 ids=next(r[1:] for r in rows if r[0]=='!Sample_geo_accession')
 data=[{'gsm':x} for x in ids]
 for r in rows:
  for d,v in zip(data,r[1:]):
   if r[0]=='!Sample_characteristics_ch1' and ': ' in v:
    k,val=v.split(': ',1); d[k]=val
   elif r[0] in ['!Sample_title','!Sample_source_name_ch1','!Sample_data_processing']:d[r[0][8:]]=v
 (B/'resources'/f'{g}_clinical_raw.json').write_text(json.dumps(data,indent=2))
 print(g,len(data),list(data[0]),flush=True)
 for k in data[0]:
  if any(w in k.lower() for w in ['necrosis','huvos','gender','age','group']):print(k,dict(__import__('collections').Counter(d.get(k) for d in data)),flush=True)
for f in ['GSE39055_paper_plain.txt','Zheng2024_fulltext_plain.txt']:
 s=(B/'sources'/f).read_text()
 for pat in ['90%','Data availability','Availability of data','github','BS3']:
  matches=list(re.finditer(re.escape(pat),s,re.I))
  for m in matches[:3]:print(f,pat,s[max(0,m.start()-200):m.end()+500].replace('\n',' '))
