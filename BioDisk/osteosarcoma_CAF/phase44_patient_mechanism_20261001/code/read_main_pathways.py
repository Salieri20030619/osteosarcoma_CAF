from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase43_full_transcriptome_discovery_20261001/delivery')
for ctx in ['supported__SCPCP000017__F3_vs_rest','expanded__SCPCP000017__F4_vs_rest']:
 for lib in ['hallmark','matrisome']:
  x=pd.read_csv(O/f'{ctx}__{lib}__nominal_GSEA.tsv',sep='\t');print('RESULT',ctx,lib);print(x[['Term','NES','NOM p-val','FDR q-val','GSEA_model_P','P_type']].to_string(index=False))