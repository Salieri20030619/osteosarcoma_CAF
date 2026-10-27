from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_identity_20260919')
s='GSE318840__GSM9504304'
a=pd.read_csv(O/'tables'/(s+'_identity.tsv.gz'),sep='\t')
g=a[a.raw_strict_candidate]
print(g.groupby(['caf_supported_candidate','caf_supported_positive_sensitivity']).agg(n=('cell_id','size'),raw_genes=('n_genes','median'),contamination=('decontx_contamination','median'),raw_Fib_genes=('raw_det_Fib','median'),corrected_Fib_ge1=('corrected_det_Fib','median'),corrected_Fib_gt0=('corrected_pos_Fib','median'),raw_Fib_pct=('raw_pct_Fib','median'),corrected_Fib_pct=('corrected_pct_Fib','median')).to_string())
print(a[a.lineage_phase3=='FibLike'].exclusion_reasons.fillna('').str.split(';').explode().value_counts().to_string())
