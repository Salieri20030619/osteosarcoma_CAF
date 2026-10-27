from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,anndata as ad,json,gzip
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase53_de_novo_reannotation_20261002';O=P/'phase54_raw_review_20261004';V=P/'review_intake_20261004/received_review'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
key=read(S/'unblinding/cell_identifier_key.tsv.gz').set_index('original_cell_id');roster=read(S/'blind_review/anonymous_cell_cluster_roster.tsv.gz').set_index('cell_id');target=set(read(O/'tables/review_target_registry.tsv').cluster_id);out=[];summ=[];status=[]
for sample,k in key.groupby('original_sample',sort=True):
 a=ad.read_h5ad(P/'phase26_client_reanalysis_20260922/inputs'/sample/'full_called_counts.h5ad',backed='r');obs=a.obs.copy();a.file.close();assert obs.index.tolist()==k.index.tolist()
 cols=[c for c in obs if any(v in c.lower() for v in ['doublet','scrublet','scdblfinder'])]
 r=roster.loc[k.cell_id].copy();r['original_cell_id']=k.index;r['original_sample']=sample
 for c in cols:r['input_'+c]=obs[c].to_numpy()
 # Keep model results separate. Do not interpret a populated False flag without a valid assessment status as a negative result.
 methods=[]
 if 'scDblFinder_class' in obs:
  val=obs.scDblFinder_class.astype(str).str.lower();valid=val.isin(['singlet','doublet']);methods.append(('provider_scDblFinder',valid.to_numpy(),val.eq('doublet').to_numpy()))
 if 'predicted_doublet' in obs:
  val=obs.predicted_doublet.astype(str).str.lower();valid=val.isin(['true','false','1','0'])
  if 'doublet_assessment_status' in obs:valid &= ~obs.doublet_assessment_status.astype(str).str.lower().str.contains('unassess|not_|not |missing|fail|unavail|skip',regex=True)
  methods.append(('upstream_predicted_doublet',valid.to_numpy(),val.isin(['true','1']).to_numpy()))
 native=P/'phase31_customer_completion_20260926/doublet'/sample
 if (native/'COMPLETE').is_file() and (native/'exit_code').read_text().strip()=='0':
  d=read(native/'doublet_results.tsv').set_index('cell_id').loc[obs.index];r['phase31_scDblFinder_class']=d.scDblFinder_class.to_numpy();r['phase31_scDblFinder_score']=d.scDblFinder_score.to_numpy();methods.append(('phase31_scDblFinder',np.ones(len(r),dtype=bool),d.scDblFinder_class.eq('doublet').to_numpy()))
 for method,valid,dbl in methods:
  for cl in r.r10.unique():
   m=r.r10.eq(cl).to_numpy();summ.append({'source_id':r.source_id.iloc[0],'sample_id':r.sample_id.iloc[0],'patient_id':r.patient_id.iloc[0],'cluster_id':cl,'method':method,'cells':int(m.sum()),'assessed':int((valid&m).sum()),'doublet_calls':int((valid&dbl&m).sum())})
 status.append({'sample':sample,'input_fields':';'.join(cols),'upstream_status_values':';'.join(sorted(obs.doublet_assessment_status.astype(str).unique())) if 'doublet_assessment_status' in obs else 'provider_native_fields_or_unavailable','methods':';'.join(x[0] for x in methods)})
 out.append(r[r.r10.isin(target)].reset_index());print(sample,flush=True)
write(pd.concat(out,ignore_index=True),'target_cell_doublet_evidence.tsv.gz');write(pd.DataFrame(summ),'doublet_sample_cluster_evidence.tsv');write(pd.DataFrame(status),'doublet_provenance.tsv')
(O/'DOUBLET_METADATA_COMPLETE').write_text('SUCCESS\n')