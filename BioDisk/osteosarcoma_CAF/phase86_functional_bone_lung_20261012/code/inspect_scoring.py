from pathlib import Path
import pandas as pd,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
for prefix in ['phase63','phase64','phase65','phase67','phase69','phase70']:
 for d in P.glob(prefix+'*'):
  print('FILES',d.name,[str(f.relative_to(d)) for f in d.rglob('*') if f.is_file() and ('resources/' in str(f) or 'code/' in str(f))][:35])
ob=pd.read_csv(P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz',sep='\t',keep_default_na=False)
q=ob[ob.dataset.isin(['SCPCP000017','SCPCP000023'])&ob.strict_site.isin(['Bone','Lung'])]
print('CONSERVATIVE',q.groupby(['dataset','strict_site','conservative_fibro']).size().to_string())
print('QUALITY',q[['raw__collagen_biosynthesis','corrected__collagen_biosynthesis','mitochondrial_fraction','detected_features']].describe().to_string())
for d in P.glob('phase6*'):
 for f in (d/'code').glob('*.py'):
  lines=f.read_text(errors='replace').splitlines()
  hits=[i for i,s in enumerate(lines) if 'collagen_biosynthesis' in s and ('=' in s or ':' in s)]
  if hits:
   print('SCORING_CODE',str(f))
   for i in hits[:2]:print('\n'.join(lines[max(0,i-4):i+15]))
