from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase62_F1F5_internal_20261005');jobs=json.loads((O/'resources/GSEA_JOBS.json').read_text());files=[j for j in jobs if Path(j['output']).exists()];print('COMPLETED_FILES',len(files),'TOTAL',len(jobs))
for j in files:
 if j['scope']!='all' or j['library']=='reactome':continue
 g=pd.read_csv(j['output'],sep='\t');g=g[(g.NES>0)&(g.GSEA_model_P<.05)].sort_values(['model_BH_FDR','NES'],ascending=[True,False]).head(4)
 print(j['contrast'],j['library'],'; '.join(f'{r.Term}: {r.NES:.2f},P={r.GSEA_model_P:.3g},q={r.model_BH_FDR:.3g}' for r in g.itertuples()))
