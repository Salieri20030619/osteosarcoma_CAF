from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase8_original_plan_20260919';D=P/'phase7_review_states_20260917'
panels={'Gong_human_PTGS2_coexpression_reference':['PTGS2','IL6','CCL2','CXCL1','PTX3','TNFAIP6'],'Alveolar_reference':['SPINT2','FGFR4','GPC3','NPNT','CES1'],'Adventitial_reference':['SFRP2','PI16','SERPINF1']}
(O/'resources/lung_reference_protocol.json').write_text(json.dumps({'panels':panels,'Gong_source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9830653/','Gong_scope':'six explicitly named human lung PTGS2 co-expression genes; not the complete murine signature','other_source':str(D/'resources/lung_reference_signatures.json'),'rule':'descriptive fraction detecting >=2 measured genes; overlap retained; no inferred cell origin or migration; no normal atlas classifier'},indent=2))
files=[D/'integration_review/strict_CAF_harmony_review.h5ad']+sorted((P/'geo_identity_20260919/objects').glob('*.h5ad'))
rows=[]
for p in files:
 a=ad.read_h5ad(p);symbols=a.var.symbol.astype(str).to_numpy();x=a.layers['counts'].tocsc();tot=np.asarray(x.sum(1)).ravel()
 for sample,ii in a.obs.groupby('sample',observed=True).indices.items():
  for name,genes in panels.items():
   detected=np.zeros(len(ii),int);expr=np.zeros(len(ii));covered=0
   for g in genes:
    ix=np.flatnonzero(symbols==g)
    if not len(ix):continue
    covered+=1;v=np.asarray(x[ii][:,ix].sum(1)).ravel();detected+=v>0;expr+=np.log1p(v/np.maximum(tot[ii],1)*1e4)
   rows.append({'sample':sample,'dataset':str(a.obs.iloc[ii[0]].dataset),'site':str(a.obs.iloc[ii[0]].site),'cells':len(ii),'program':name,'covered':covered,'requested':len(genes),'two_gene_positive_cells':int((detected>=2).sum()),'two_gene_positive_pct':float((detected>=2).mean()*100),'mean_log1p_CP10K':float(expr.mean()/max(covered,1)),'interpretation':'reference-marker resemblance; not confirmed lung-resident origin'})
 del a,x
pd.DataFrame(rows).to_csv(O/'tables/lung_reference_phenotype.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
