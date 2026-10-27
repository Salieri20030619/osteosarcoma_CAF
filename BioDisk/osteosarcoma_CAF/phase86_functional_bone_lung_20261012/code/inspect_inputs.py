from pathlib import Path
import json,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
Q=P/'phase58_F1F5_round1_20261004/resources'
sets=json.loads((Q/'programs_frozen.json').read_text())['sets']
keys=[k for k in sets if any(s in k.upper() for s in ['CUSTOMER','COLLAGEN','MATRISOME','TGF_BETA','TGFB','IL6_JAK','INFLAMMATORY_RESPONSE','SMOOTH_MUSCLE','CROSSLINK'])]
print('PROGRAMS',json.dumps({k:sets[k] for k in keys},ensure_ascii=False))
libs=json.loads((Q/'libraries_canonical.json').read_text())
print('LIB_TERMS',json.dumps({lib:[k for k in terms if any(s in k.upper() for s in ['COLLAGEN','CROSSLINK','MATRISOME','SMOOTH_MUSCLE'])] for lib,terms in libs.items()}))
ob=pd.read_csv(P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz',sep='\t',keep_default_na=False)
print('CELL_COLUMNS',ob.columns.tolist())
print('IDENTITY',ob.review_lineage_code.value_counts().to_dict())
print('BONE_LUNG_COUNTS',ob[ob.dataset.isin(['SCPCP000017','SCPCP000023'])].groupby(['dataset','strict_site','review_lineage_code']).size().to_string())
cm=pd.read_csv(P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv',sep='\t',keep_default_na=False)
print('CLINICAL',cm[cm.dataset.isin(['SCPCP000017','SCPCP000023'])].groupby(['dataset','provider_disease_timing','treatment']).size().to_string())
print('PHASES',[f.name for f in P.iterdir() if f.is_dir() and f.name.startswith(('phase67','phase68','phase69','phase70'))])
print('STATE',json.loads((P/'CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json').read_text()))
