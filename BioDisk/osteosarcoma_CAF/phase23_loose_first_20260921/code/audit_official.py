from pathlib import Path
import json
import pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase23_loose_first_20260921'
(O/'tables').mkdir(exist_ok=True);(O/'resources').mkdir(exist_ok=True)
rows=[];cells=[]
for r in pd.read_csv(P/'phase16_scpca_supplement_20260920/manifest.tsv',sep='\t').to_dict('records'):
 a=ad.read_h5ad(r['processed'],backed='r');sm=r['scpca_sample_id'];ob=a.obs.copy()
 fields=[c for c in ob if c.startswith(('consensus_','singler_','cellassign_','scimilarity_','infercnv_','is_infercnv'))]+['scDblFinder_class','scDblFinder_score','miQC_pass']
 ob=ob[fields];ob['sample']=sm;ob['cell_id']=[sm+'_'+str(z) for z in ob.index];cells.append(ob)
 row={'sample':sm,'processed_cells':a.n_obs,'reference_cells':str(a.uns.get('infercnv_num_reference_cells','missing')),'infercnv_status':str(a.uns.get('infercnv_status','missing')),'CNV_cells':int(pd.to_numeric(ob.get('infercnv_total_cnv',pd.Series(dtype=float)),errors='coerce').notna().sum())}
 rows.append(row)
 meta={k:str(v) for k,v in a.uns.items() if k.startswith(('infercnv_','min_gene_','prob_compromised_','celltype_','consensus_','scpca_filter_'))}
 (O/'resources'/f'{sm}_official_provenance.json').write_text(json.dumps(meta,indent=2));a.file.close()
pd.DataFrame(rows).to_csv(O/'tables/official_ScPCA_CNV_availability.tsv',sep='\t',index=False)
pd.concat(cells,ignore_index=True).to_csv(O/'tables/official_ScPCA_cell_evidence.tsv.gz',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
