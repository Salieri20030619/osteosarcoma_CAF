from pathlib import Path
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
p=B/'review_intake_20261004/freeze_annotation_phase55.py'
ss=p.read_text().splitlines()
inds=set()
for i,l in enumerate(ss):
 if any(k in l for k in ['doublet','candidate_tier','lineage_confidence']):inds.update(range(max(0,i-3),min(len(ss),i+5)))
print('\n'.join(f'{i+1}: {ss[i]}' for i in sorted(inds)))
