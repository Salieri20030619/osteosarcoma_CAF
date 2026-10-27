from pathlib import Path
import pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
h=pd.read_csv(B/'phase41_full_transcriptome_setup_20261001/resources/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
print(h[h.symbol.isin(['PXDN','ACAN','PRG2','SLC25A12'])][['symbol','entrez_id','refseq_accession']].to_string(index=False))
p=O/'code/verify_response.py';s=p.read_text()
needle="assert cov[cov.program=='fibrillar_I_III_V'].evaluable.all()"
replacement="""mapping=pd.read_csv(O/'tables/GSE21257_probe_mapping.tsv',sep='\\t').set_index('ID')
assert mapping.loc['9ilipx9ZCOSmZZ9Vhw','canonical_symbol']=='PXDN'
assert mapping.loc['NSjRKdq2eSGf0ur4aQ','canonical_symbol']=='PRG2'
assert mapping.loc['iVWUrXuSKJXuieoXkI','canonical_symbol']=='SLC25A12'
assert mapping.loc['Bbg0g0HnkhSis7iip4','canonical_symbol']=='ACAN'
"""+needle
assert needle in s
p.write_text(s.replace(needle,replacement))
p=O/'code/fix_annotation_ids.py';s=p.read_text();s=s[:s.index("print(h[h.symbol")];p.write_text(s)
print('QA_ID_ASSERTIONS_ADDED')
