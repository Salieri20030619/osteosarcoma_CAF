from pathlib import Path
import json,pandas as pd,anndata as ad,numpy as np,hashlib
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'review_intake_20261005_checklist'
a=pd.read_csv(B/'phase63_matrix_resistance_20261005/tables/cell_mechanism_scores_with_confidence.tsv.gz',sep='\t',low_memory=False)
out=[];sources=[]
for sample,t in a.groupby('sample',sort=True):
 p=B/'phase26_client_reanalysis_20260922/inputs'/sample/'full_called_counts.h5ad'
 if not p.exists():sources.append({'sample':sample,'status':'path_missing','path':str(p)});continue
 obj=ad.read_h5ad(p,backed='r');obs=obj.obs.copy();obj.file.close()
 assert set(t.cell_id)<=set(obs.index),(sample,list(set(t.cell_id)-set(obs.index))[:3])
 obs=obs.loc[t.cell_id];methods=[]
 if 'scDblFinder_class' in obs:
  val=obs.scDblFinder_class.astype(str).str.lower();methods.append(('provider_scDblFinder',val.isin(['singlet','doublet']).to_numpy(),val.eq('doublet').to_numpy()))
 if 'predicted_doublet' in obs:
  val=obs.predicted_doublet.astype(str).str.lower();valid=val.isin(['true','false','1','0'])
  if 'doublet_assessment_status' in obs:valid &= ~obs.doublet_assessment_status.astype(str).str.lower().str.contains('unassess|not_|not |missing|fail|unavail|skip',regex=True)
  methods.append(('upstream_predicted_doublet',valid.to_numpy(),val.isin(['true','1']).to_numpy()))
 native=B/'phase31_customer_completion_20260926/doublet'/sample
 if (native/'COMPLETE').is_file() and (native/'exit_code').read_text().strip()=='0':
  d=pd.read_csv(native/'doublet_results.tsv',sep='\t').set_index('cell_id').loc[t.cell_id]
  methods.append(('phase31_scDblFinder',np.ones(len(t),bool),d.scDblFinder_class.eq('doublet').to_numpy()))
 assessed=np.zeros(len(t),bool);positive=np.zeros(len(t),bool)
 for name,valid,pos in methods:assessed|=valid;positive|=valid&pos
 for j,row in enumerate(t.itertuples()):out.append({'cell_id':row.cell_id,'sample':sample,'dataset':row.dataset,'class_id':row.class_id,'subgroup_id':row.subgroup_id,'exported_assessed':bool(pd.notna(row.doublet_assessed_any) and row.doublet_assessed_any==True),'source_assessed':bool(assessed[j]),'source_positive':bool(positive[j]),'source_methods':';'.join(name for name,v,p in methods if v[j])})
 sources.append({'sample':sample,'status':'checked','path':str(p),'target_cells':len(t),'assessed':int(assessed.sum()),'positive':int(positive.sum()),'methods':';'.join(x[0] for x in methods)})
d=pd.DataFrame(out);assert len(d)==8174 and d.cell_id.is_unique
d.to_csv(O/'tables/doublet_source_reconciliation.tsv.gz',sep='\t',index=False);pd.DataFrame(sources).to_csv(O/'resources/doublet_source_inputs.tsv',sep='\t',index=False)
rows=[]
for level in ['class_id','subgroup_id']:
 for group,t in d.groupby(level):
  rows.append({'level':level,'group':group,'cells':len(t),'exported_assessed':int(t.exported_assessed.sum()),'source_assessed':int(t.source_assessed.sum()),'source_unknown':int((~t.source_assessed).sum()),'source_positive':int(t.source_positive.sum()),'recovered_from_existing_source':int((~t.exported_assessed&t.source_assessed).sum())})
r=pd.DataFrame(rows);r.to_csv(O/'tables/doublet_coverage_reconciled.tsv',sep='\t',index=False)
print('TOTAL',d[['exported_assessed','source_assessed','source_positive']].sum().to_dict(),flush=True)
print(r[(r.level=='class_id')|r['group'].isin(['F2_C1','F2_C3','F4_C8','F4_C10'])].to_string(index=False),flush=True)
(O/'resources/DOUBLET_RECONCILIATION_STATUS.json').write_text(json.dumps({'cells':len(d),'exported_assessed':int(d.exported_assessed.sum()),'source_assessed':int(d.source_assessed.sum()),'source_positive':int(d.source_positive.sum()),'source_unknown':int((~d.source_assessed).sum()),'analysis_labels_changed':False,'new_doublet_detection_run':False},indent=2))
