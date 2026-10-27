import pathlib,json,hashlib
import pandas as pd,numpy as np
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');q=root/'phase95_NCI_expanded_reanalysis_20261016';p=root/'phase100_NCI_completion_20261020'
for d in ['code','resources','tables','figures','logs','objects','delivery']:(p/d).mkdir(parents=True,exist_ok=True)
run=json.loads((q/'COMPUTATION_FINISHED.json').read_text(encoding='utf-8'))
rows=[]
for rr in sorted(run['samples'],key=lambda z:z['sample']):
 sm=rr['sample'];c=pd.read_csv(q/'tables'/(sm+'_all_cell_evidence.tsv.gz'),sep='\t',encoding='utf-8');mk=pd.read_csv(q/'tables'/(sm+'_cluster_markers.tsv'),sep='\t',encoding='utf-8')
 for cl,z in c.groupby('cluster'):
  m=mk[mk.cluster.astype(str).eq(str(cl))].sort_values('rank');top=';'.join(m.symbol.astype(str).head(16))
  k=z[z.mechanical_expanded_candidate]
  scores={a:float(z[a].median()) for a in z if a.startswith('phase27_')}
  top_panels=';'.join(a.replace('phase27_','')+':'+str(round(v,2)) for a,v in sorted(scores.items(),key=lambda a:-a[1])[:3])
  rows.append(dict(sample=sm,cluster=str(cl),cells=len(z),mechanical_candidates=len(k),reference_pass=rr['reference_pass'],author_CAF=int(z.author_CAF.sum()),author_tumor=int(z.Ann_Level0.eq('Tumor').sum()),fibro_winner=int(z.provisional_lineage.eq('Fibroblast').sum()),native_diploid=int(z.native_copykat_call.eq('diploid').sum()),native_aneuploid=int(z.native_copykat_call.eq('aneuploid').sum()),candidate_author_CAF=int(k.author_CAF.sum()),candidate_author_tumor=int(k.Ann_Level0.eq('Tumor').sum()),candidate_ambient_fibro=int(k.decontX_lineage_winner.eq('Fibroblast').sum()) if 'decontX_lineage_winner' in k else 0,top_panels=top_panels,top_markers=top))
df=pd.DataFrame(rows);df.to_csv(p/'tables/cluster_review_evidence.tsv',sep='\t',encoding='utf-8',index=False)
for sm,z in df.groupby('sample'):
 print(sm,'clusters',len(z),'mechanical',z.mechanical_candidates.sum(),'authorCAF',z.author_CAF.sum(),'candidate_clusters',(z.mechanical_candidates>0).sum())
print('CELL_COLUMNS',c.columns.tolist())
print('manual first', (q/'tables/NCI__CCDI0048_5864_manual_cluster_review.tsv').read_text(encoding='utf-8'))
tree=json.loads((q/'resources/author_repo_tree.json').read_text())
for x in tree['tree']:
 if any(w in x['path'].lower() for w in ['sample','patient','clinical','nci','metadata','supplement']):print('REPOFILE',x['path'])
