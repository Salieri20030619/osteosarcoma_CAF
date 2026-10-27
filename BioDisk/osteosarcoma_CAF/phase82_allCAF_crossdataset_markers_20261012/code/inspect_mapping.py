from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase82_allCAF_crossdataset_markers_20261012';B=P/'phase81_allCAF_bone_lung_20261012'
def r(p):return pd.read_csv(p,sep='\t',low_memory=False)
primary=r(O/'tables/primary_same_layer_crossdataset_markers.tsv');print('PRIMARY\n'+primary.to_string(index=False))
allstats=r(B/'tables/all_genes_statistics.tsv.gz');d=allstats[(allstats.sensitivity=='all_cells')&(allstats.layer=='expanded')&(allstats.design=='independent')&(allstats.P<.05)]
z=d[d.dataset.eq('SCPCP000017')].merge(d[d.dataset.eq('SCPCP000023')],on='feature',suffixes=('_17','_23'));z=z[z.mean_lung_minus_bone_17*z.mean_lung_minus_bone_23>0]
can=r(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv')
print('MAPPING\n'+can[(can.symbol.isin(z.symbol_17)|can.feature.isin(z.feature))&can.dataset.isin(['SCPCP000017','SCPCP000023'])].to_string(index=False))
print('RAW\n'+z[['feature','symbol_17','symbol_23']].to_string(index=False))
reg=r(O/'tables/MARKER_EVIDENCE_REGISTRY.tsv');print('EVIDENCE\n'+reg.to_string(index=False,max_colwidth=20))
surf=r(O/'tables/candidate_existing_membrane_annotation.tsv');print('SURFACE_COLS',surf.columns.tolist());print(surf.head().to_string(index=False,max_colwidth=30))
