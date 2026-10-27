from pathlib import Path
import anndata as ad,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for p in [P/'phase7_review_states_20260917/integration_review/strict_CAF_harmony_review.h5ad',P/'phase7_review_states_20260917/analysis/mycaf_subclusters.h5ad']:
 a=ad.read_h5ad(p,backed='r')
 print('FILE',p,'SHAPE',a.shape,'OBS',a.obs.columns.tolist(),'VAR',a.var.head().to_dict(),'LAYERS',list(a.layers))
 print(a.obs.head(2).to_string());a.file.close()
print(pd.read_csv(P/'phase7_review_states_20260917/analysis/mycaf_subcluster_patient_counts.tsv',sep='\t').head(12).to_string())
