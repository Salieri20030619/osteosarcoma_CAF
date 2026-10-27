from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase77_C3_external_function_20261010';cells=pd.read_csv(O/'tables/query_cell_scores.tsv.gz',sep='\t').set_index('cell_id');gs=['COL3A1','COL5A1','COL5A2','COL16A1','EMILIN1','ADAMTS2','PCOLCE','LOXL2'];parts=[]
(O/'resources/POSTHOC_GENE_DECOMPOSITION.json').write_text(json.dumps({'reason':'Diagnose opposite processing-module directions in the two external primary patients','created_after_main_results':True,'scope':'Display all eight prespecified module genes, per-patient expression and detection, no new gene-level P tests, no mapping or thresholds changed'},indent=2))
a=ad.read_h5ad(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad');mp=pd.read_csv(P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',sep='\t')
for ds,d in cells[cells.dataset!='GSE237070'].groupby('dataset'):
 ii=a.obs_names.get_indexer(d.index);assert (ii>=0).all();gm=mp[mp.dataset==ds].drop_duplicates('symbol').set_index('symbol').feature;js=a.var_names.get_indexer(gm.loc[gs]);assert (js>=0).all();counts=a.X[ii][:,js].toarray();out=d.copy()
 for j,g in enumerate(gs):out[g+'_count']=counts[:,j];out[g+'_score']=np.log1p(counts[:,j]/out.library.to_numpy()*1e4)
 parts.append(out)
for path in sorted((P/'phase37_requested_source_search_20260927/objects/GSE237070_reconstructed').glob('*/full_native_called_counts.h5ad')):
 d=cells[cells['sample']==path.parent.name];q=ad.read_h5ad(path);ii=q.obs_names.get_indexer(d.index);assert (ii>=0).all();gm=pd.Series(q.var_names,index=q.var.symbol);assert not gm.loc[gs].index.duplicated().any();js=q.var_names.get_indexer(gm.loc[gs]);counts=q.X[ii][:,js].toarray();out=d.copy()
 for j,g in enumerate(gs):out[g+'_count']=counts[:,j];out[g+'_score']=np.log1p(counts[:,j]/out.library.to_numpy()*1e4)
 parts.append(out)
d=pd.concat(parts);assert np.allclose(d.processing2,(d.PCOLCE_score+d.LOXL2_score)/2);assert np.allclose(d.collagen_core5,d[[g+'_score' for g in ['COL3A1','COL5A1','COL5A2','PCOLCE','LOXL2']]].mean(1));d.to_csv(O/'tables/posthoc_cell_gene_scores.tsv.gz',sep='\t');rows=[]
for (ds,tissue,donor,site),z in d.groupby(['dataset','tissue','donor_record','site']):
 l=z[z.C3like];r=z[~z.C3like]
 if not len(l) or not len(r):continue
 for g in gs:
  rows.append({'dataset':ds,'tissue':tissue,'donor':donor,'site':site,'gene':g,'target_cells':len(l),'reference_cells':len(r),'target_mean_score':l[g+'_score'].mean(),'reference_mean_score':r[g+'_score'].mean(),'score_effect':l[g+'_score'].mean()-r[g+'_score'].mean(),'target_detected':int((l[g+'_count']>0).sum()),'reference_detected':int((r[g+'_count']>0).sum())})
r=pd.DataFrame(rows);r.to_csv(O/'tables/posthoc_patient_site_gene_decomposition.tsv',sep='\t',index=False);print(r[(r.dataset=='GSE237070')&(r.tissue=='primary')][['donor','gene','target_cells','reference_cells','score_effect','target_detected','reference_detected']].to_csv(index=False));(O/'POSTHOC_COMPLETE.json').write_text(json.dumps({'status':'completed','records':len(r),'module_scores_reproduced':True},indent=2))
