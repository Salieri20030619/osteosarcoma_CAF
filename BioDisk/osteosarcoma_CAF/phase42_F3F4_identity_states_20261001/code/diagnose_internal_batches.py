from pathlib import Path
import anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
a=ad.read_h5ad(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',backed='r')
for layer,label in [('supported','supported_class_id'),('expanded','expanded_class_id')]:
 for cls in ['F3','F4']:
  g=a.obs[a.obs[label].astype(str).eq(cls)]
  print('SUBSET',layer,cls,len(g),g.groupby('study_technology',observed=True).size().to_dict())
a.file.close()