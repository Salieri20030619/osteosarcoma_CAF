
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'patient_expansion_download_investigation_20261013_v1';d=pd.read_csv(O/'tables/S23_additional_Fibro_expression_cells.tsv.gz',sep='\t');lineage=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text());all14=json.loads((P/'CAF_ALL14_RESEARCH_CURRENT.json').read_text())['retained_genes'];genes=sorted(set(sum(lineage.values(),[])+all14+['ACTA2','TAGLN','FAP','PDPN','CD248','GPC1','THY1','POSTN','COL11A1']));rows=[];cells=[];refs=[]
for sm,g in d.groupby('sample'):
 ev=pd.read_csv(P/'phase27_identity_review_20260926/tables'/f'{sm}_review.tsv.gz',sep='\t',low_memory=False).set_index('cell_id');ids=g.original_cell_id.astype(str).tolist();assert set(ids)<=set(ev.index)
 f=P/'phase26_client_reanalysis_20260922/inputs'/sm/'full_called_counts.h5ad';a=ad.read_h5ad(f,backed='r');assert set(ids)<=set(a.obs_names);symbols=a.var.symbol.astype(str);mapping={}
 for j,sy in enumerate(symbols):
  for s in sy.split(';'):
   if s in genes:mapping.setdefault(s,[]).append(j)
 js=sorted(set(j for v in mapping.values() for j in v));ix=a.obs_names.get_indexer(ids);order=np.argsort(ix);rev=np.argsort(order);xx=a[ix[order],js].X.tocsr()[rev];jmap={v:i for i,v in enumerate(js)};lib=ev.loc[ids,'n_counts'].to_numpy(float);a.file.close()
 gg=g.copy();gg['known_doublet']=ev.loc[ids,'known_doublet'].astype(str).str.lower().isin(['true','1','1.0']).to_numpy();gg['provisional_lineage_stage27']=ev.loc[ids,'provisional_lineage'].to_numpy();gg['copykat_call']=ev.loc[ids,'copykat_call'].to_numpy();gg['infercnv_call']=ev.loc[ids,'infercnv_call'].to_numpy();gg['full_library']=lib
 for gene in genes:
  indices=[jmap[j] for j in mapping.get(gene,[])];raw=np.asarray(xx[:,indices].sum(1)).ravel() if indices else np.full(len(g),np.nan);val=np.log1p(raw/np.maximum(lib,1)*1e4)
  gg['RNA_log1pCP10k__'+gene]=val
  rows.append(dict(sample=sm,donor_record=g.donor_record.iloc[0],site=g.strict_site_group.iloc[0],gene=gene,cells=len(g),assayed=bool(indices),detected_fraction=float(np.nanmean(raw>0)) if indices else None,mean_log1pCP10k=float(np.nanmean(val)) if indices else None,known_doublet_cells=int(gg.known_doublet.sum())))
 cells.append(gg)
 print('S23_MARKER_REVIEW',sm,len(g),flush=True)
pd.DataFrame(rows).to_csv(O/'tables/S23_expression_Fibro_marker_evidence.tsv',sep='\t',index=False);cc=pd.concat(cells);cc.to_csv(O/'tables/S23_expression_Fibro_cells_with_markers.tsv.gz',sep='\t',index=False)
su=cc.groupby(['strict_site_group','donor_record']).agg(expression_Fibro_cells=('cell_id','size'),known_doublet_cells=('known_doublet','sum')).reset_index();su['not_known_doublet']=su.expression_Fibro_cells-su.known_doublet_cells;su.to_csv(O/'tables/S23_expression_Fibro_patient_capacity.tsv',sep='\t',index=False)
r={'status':'completed','cells_reviewed':len(cc),'bone_patients_with_expression_Fibro':int(su[su.strict_site_group=='Bone'].donor_record.nunique()),'lung_patients_with_expression_Fibro':int(su[su.strict_site_group=='Lung'].donor_record.nunique()),'lung_patients_with_nonflagged_expression_Fibro':int(su[(su.strict_site_group=='Lung')&su.not_known_doublet.gt(0)].donor_record.nunique()),'no_new_identity_gates_or_P_values':True,'baseline_changed':False}
(O/'resources/S23_MARKER_INVESTIGATION_COMPLETE.json').write_text(json.dumps(r,indent=2));print(json.dumps(r));print(su.to_string(index=False))
