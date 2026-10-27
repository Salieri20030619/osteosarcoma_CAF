from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
q=P/'geo_preprocess_20260918'
pack=pd.read_csv(q/'review/cluster_review_packets.tsv',sep='\t',dtype={'cluster':str})
rows=[]
for p in sorted((q/'results').glob('*/cluster_markers.tsv')):
 t=pd.read_csv(p,sep='\t',dtype={'group':str})
 t['rank_within_cluster']=t.groupby('group').cumcount()+1
 hit=t[t.symbol.str.upper().isin(['SR1004','S100A4'])].copy()
 for r in hit.to_dict('records'):
  lab=pack[(pack['sample']==p.parent.name)&(pack.cluster==r['group'])]
  r['sample']=p.parent.name;r['lineage']=lab.proposed_lineage.iloc[0] if len(lab) else 'unknown';rows.append(r)
print(pd.DataFrame(rows).to_string(index=False) if rows else 'NO_EXACT_MARKER_HITS')
