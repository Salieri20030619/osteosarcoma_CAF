from pathlib import Path
import pandas as pd,io
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
h=pd.read_csv(B/'phase41_full_transcriptome_setup_20261001/resources/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
print(h[h.symbol.isin(['ACAN','PXDN'])][['symbol','prev_symbol','alias_symbol']].to_string(index=False))
t=(O/'sources/GPL10295_platform_only.soft').read_text().split('!platform_table_begin')[1].split('!platform_table_end')[0].strip()
a=pd.read_csv(io.StringIO(t),sep='\t',dtype=str,keep_default_na=False)
for g in ['ACAN','PXDN']:
 r=h[h.symbol==g].iloc[0];ss={g}|set(r.prev_symbol.split('|'))|set(r.alias_symbol.split('|'));ss.discard('')
 print('PLATFORM',g,a[a.Symbol.isin(ss)][['ID','Symbol','GB_ACC']].to_string(index=False))
 for sym in ss:
  print('SYMBOL_CONFLICT',sym,h[h.prev_symbol.str.split('|').map(lambda z:sym in z)|h.alias_symbol.str.split('|').map(lambda z:sym in z)][['symbol','prev_symbol','alias_symbol']].to_string(index=False))
