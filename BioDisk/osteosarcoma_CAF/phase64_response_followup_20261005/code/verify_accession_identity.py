from pathlib import Path
import json,datetime,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
h=pd.read_csv(B/'phase41_full_transcriptome_setup_20261001/resources/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
print(h[h.symbol.isin(['PLPPR3','PXDN','ACAN'])][['symbol','refseq_accession']].to_string(index=False))
p=O/'code/verify_response.py';s=p.read_text().replace("mapping.loc['9ilipx9ZCOSmZZ9Vhw','canonical_symbol']=='PXDN'","mapping.loc['9ilipx9ZCOSmZZ9Vhw','canonical_symbol']=='PLPPR3'")
s=s.replace("assert mapping.loc['iVWUrXuSKJXuieoXkI','canonical_symbol']=='SLC25A12'","assert pd.isna(mapping.loc['iVWUrXuSKJXuieoXkI','canonical_symbol'])")
p.write_text(s)
p=O/'resources/ANNOTATION_QA_AMENDMENT.json';old=json.loads(p.read_text())
(O/'resources/ANNOTATION_QA_INITIAL_UNVERIFIED_DRAFT.json').write_text(json.dumps(old,indent=2))
old.update({'verified_correction_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'GPL10295 historical PRG2 symbols include both PLPPR3 (NM_024888) and PRG2 (NM_002728); HGNC RefSeq resolves these. AGC1 is an ambiguous historical symbol; NM_001135 resolves ACAN. PXDN not present among these probes.','status':'Initial tentative attribution of NM_024888 to PXDN and NM_013227 to SLC25A12 was incorrect, corrected by accession verification before final analysis. Final mapping uses verified HGNC RefSeq/Entrez, not these tentative guesses.','unresolved':'ACAN alternate isoform NM_013227 is not in HGNC selected RefSeq field and remains unmapped by current conservative rule; ACAN NM_001135 is mapped.'})
p.write_text(json.dumps(old,indent=2));print('ANNOTATION_IDENTITY_CORRECTED')
