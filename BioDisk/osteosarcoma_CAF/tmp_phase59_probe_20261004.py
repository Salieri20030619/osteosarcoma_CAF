from pathlib import Path
import json, pandas as pd, anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
S=P/'phase58_F1F5_round1_20261004'
a=ad.read_h5ad(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',backed='r')
print('OBJECT',a.shape,a.obs.columns.tolist(),a.var.columns.tolist());a.file.close()
r=pd.read_csv(P/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz',sep='\t',low_memory=False)
for col in ['review_lineage_code','CNV_evidence_current','candidate_tier','phase54_review_decision','doublet_assessed_any','doublet_positive_any','lineage_confidence']:
 print(col,r[col].fillna('NA').value_counts().to_dict())
ob=pd.read_csv(S/'tables/frozen_expanded_cells.tsv.gz',sep='\t');q=r[r.original_cell_id.isin(ob.cell_id)].merge(ob[['cell_id','class_id']],left_on='original_cell_id',right_on='cell_id')
print('CLASS_REVIEW',q.groupby(['class_id','review_lineage_code','CNV_evidence_current']).size().to_string())
libs=json.loads((S/'resources/libraries_canonical.json').read_text())
for lib,sets in libs.items():
 for k,v in sets.items():
  if any(t in k.lower() for t in ['osteoblast','ossification','mineralization','smooth_muscle_contraction','inflammatory_response','il6_jak','core_matrisome','ecm_regulators','collagen_formation','degradation_of_the_extracellular','tgf_beta','myogenesis']):print('SET',lib,k,len(v))
print('SAMPLE',r[['sample','dataset','original_donor_record','site','sample_id']].drop_duplicates().head(12).to_string(index=False))
print('RESOURCES', [str(p) for p in (P/'phase53_de_novo_reannotation_20261002/blind_review/raw_counts').glob('*.h5ad')][:3])
