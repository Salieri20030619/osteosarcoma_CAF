from pathlib import Path
import json,pandas as pd,numpy as np,anndata as ad
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase122_HLCA_supplement_20261026');q=ad.read_h5ad(O/'objects/query_program_counts.h5ad');r=ad.read_h5ad(O/'objects/reference_http_program_counts.h5ad');p=json.loads((O/'resources/SUPPLEMENT_PROGRAMS.json').read_text());genes=list(dict.fromkeys(p['HLCA_IPF_reference']+p['mCAF_matrix']+p['myCAF_customer']));assert all(g in q.var_names and g in r.var_names for g in genes);rows=[];noncol=[]
for kind,a,den in [('query',q,q.obs.library_counts.to_numpy()),('reference',r,r.obs.raw_sum.to_numpy())]:
 raw=a.X[:,a.var_names.get_indexer(genes)].toarray();norm=np.log1p(raw*10000/den[:,None]);meta=a.obs.reset_index(drop=True).copy();groupings=[]
 if kind=='query':
  for node in ['Lung_04','Bone_07','Bone_07_s1','Lung_04_s2']:
   ix=meta.index[meta['leaf' if '_s' in node else 'parent'].eq(node)];groupings.append((node,meta.loc[ix].groupby(['source','donor'],observed=True)))
 else:
  for (condition,typ),g in meta.groupby(['lung_condition','ann_finest_level'],observed=True):groupings.append((condition+' | '+typ,g.groupby(['study','donor_id'],observed=True)))
 for group,gb in groupings:
  for (study,donor),g in gb:
   ix=g.index.to_numpy();means=norm[ix].mean(0);fracs=(raw[ix]>0).mean(0)
   for j,gene in enumerate(genes):rows.append(dict(kind=kind,group=group,study=study,donor=donor,cells=len(ix),gene=gene,mean_logCP10k=float(means[j]),fraction_detected=float(fracs[j])))
   five=[genes.index(v) for v in p['HLCA_IPF_reference'] if v!='COL1A1'];noncol.append(dict(kind=kind,group=group,study=study,donor=donor,cells=len(ix),IPF_without_COL1A1=float(means[five].mean())))
df=pd.DataFrame(rows);df.to_csv(O/'tables/marker_patient_evidence.tsv',sep='\t',index=False);nc=pd.DataFrame(noncol);nc.to_csv(O/'tables/IPF_without_COL1A1_patient_scores.tsv',sep='\t',index=False);summ=[]
for (kind,group,gene),g in df.groupby(['kind','group','gene']):
 g=g[g.cells.ge(5 if kind=='query' else 20)]
 if not len(g):continue
 v=g[['mean_logCP10k','fraction_detected']].mean() if kind=='query' else g.groupby('study')[['mean_logCP10k','fraction_detected']].mean().mean()
 summ.append(dict(kind=kind,group=group,gene=gene,patients=len(g),studies=g.study.nunique(),mean_logCP10k=float(v.mean_logCP10k),fraction_detected=float(v.fraction_detected)))
pd.DataFrame(summ).to_csv(O/'tables/marker_plot_summary.tsv',sep='\t',index=False);comparisons=[]
for node,a in nc[(nc.kind=='query')&nc.cells.ge(5)].groupby('group'):
 for group,b in nc[(nc.kind=='reference')&nc.cells.ge(20)].groupby('group'):
  pct=[np.mean([(g.IPF_without_COL1A1<=v).mean() for _,g in b.groupby('study')]) for v in a.IPF_without_COL1A1]
  comparisons.append(dict(node_id=node,reference_group=group,query_patients=len(a),reference_patients=len(b),difference=float(a.IPF_without_COL1A1.median()-b.groupby('study').IPF_without_COL1A1.mean().mean()),median_reference_percentile=float(np.median(pct))))
c=pd.DataFrame(comparisons);c.to_csv(O/'tables/IPF_without_COL1A1_background.tsv',sep='\t',index=False);assert c.median_reference_percentile.between(0,1).all();print(c[c.reference_group.isin(['Healthy | Alveolar fibroblasts','IPF | Alveolar fibroblasts'])].to_string(index=False));(O/'resources/GENE_EVIDENCE_COMPLETE.json').write_text(json.dumps(dict(status='PASS',genes=genes,patient_rows=len(df),purpose='Check individual IPF markers and exclude COL1A1-driven score interpretation')),encoding='utf-8')
