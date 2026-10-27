from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026')
x=pd.read_csv(O/'tables/historical_mapping_overlap.tsv.gz',sep='\t')
def broad(v):
 if 'fibroblast' in v.lower():return 'Fibroblast'
 if v=='Pericytes':return 'Pericyte'
 if v.startswith(('Smooth muscle','SM activated')):return 'Smooth_muscle'
 return 'Other_lineage'
x['old_broad']=x.ann_finest_level_pred.map(broad)
rows=[]
for n in ['ALL','Lung_04','Bone_07','Bone_07_s1','Lung_04_s2']:
 g=x if n=='ALL' else x[x['leaf' if '_s' in n else 'parent'].eq(n)]
 rows.append(dict(node=n,cells=len(g),fine_agreement=float(g.finest_pred.eq(g.ann_finest_level_pred).mean()),broad_agreement=float(g.broad_pred.eq(g.old_broad).mean()),old_types=g.old_broad.value_counts().to_dict(),current_types=g.broad_pred.value_counts().to_dict()))
pd.crosstab(x.ann_finest_level_pred,x.finest_pred).to_csv(O/'tables/historical_mapping_confusion.tsv',sep='\t')
(O/'resources/HISTORICAL_MAPPING_REVIEW.json').write_text(json.dumps(rows,indent=2),encoding='utf-8');print(json.dumps(rows))
for p in (O/'tables').glob('mapping_history*elbo_train*'):
 t=pd.read_csv(p,sep='\t');print(p.name,t.head(2).to_dict('records'),t.tail(10).to_dict('records'))
print('REFERENCE_LOG', (O/'logs/reference_extract_attempt1.log').read_text()[-1800:])
