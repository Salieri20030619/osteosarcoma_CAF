from pathlib import Path
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918');O=P/'result_review_20260919';rows=[]
s=pd.read_csv(P/'review/input_summary.tsv',sep='\t')
for sm in s[s.cnv_eligible]['sample']:
 a=pd.read_csv(P/'inputs'/sm/'cells.tsv',sep='\t');n=a[a.cnv_role.str.startswith(('ref_','holdout_'))]
 ds=a.scDblFinder_score if a.scDblFinder_score.notna().all() else a.doublet_score
 assert (n.n_genes>=500).all() and (n.n_counts>=750).all() and (ds.loc[n.index]<=ds.quantile(.9)).all()
 assert (n.pctUMI_competing_mes<=.2).all() and ((n.pctUMI_Epi+n.pctUMI_Myogenic)<=.2).all() and (n.pctUMI_OC_specific<=.1).all() and (~n.CTA_context).all()
 assert n.lineage_phase3.isin(['T_NK','B_Plasma','Myeloid','DC','Endo']).all()
 for typ,g in n.groupby('lineage_phase3'):
  assert (g['det_'+typ+'_review']>=3).all() and (g['pctUMI_'+typ]>=.5).all()
  if typ!='Endo':assert (g.pctUMI_Endo<=.2).all()
  else:assert (g[['pctUMI_'+t for t in ['T_NK','B_Plasma','Myeloid','DC']]].sum(axis=1)<=.2).all()
 ann=pd.read_csv(P/'inputs'/sm/'annotations.tsv',sep='\t',header=None,names=['tool_id','annotation'])
 expected=set(a[a.cnv_role.str.startswith('ref_')].tool_id);observed=set(ann[ann.annotation.str.startswith('ref_')].tool_id)
 assert expected==observed and not set(a[a.cnv_role.str.startswith('holdout_')].tool_id)&observed
 rows.append(dict(sample=sm,screened_normals=len(n),reference_rules_pass=True,no_holdout_baseline_overlap=True))
pd.DataFrame(rows).to_csv(O/'tables/reference_gate_checks.tsv',sep='\t',index=False)
print('REFERENCE_GATE_PASS',len(rows))
