from pathlib import Path
import json,sys,datetime
import pandas as pd,numpy as np,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase32_lung_sensitivity_20260926';M=P/'phase30_two_layer_classification_20260926';kind=sys.argv[1]
def writej(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str))
def now():return datetime.datetime.now().astimezone().isoformat()
def read(f):return pd.read_csv(f,sep='\t',low_memory=False)
if kind=='hlca':
 qa=json.loads((O/'delivery/HLCA_QA.json').read_text());assert (O/'HLCA_COMPLETE').exists()
 calls=read(O/'tables/HLCA_cell_predictions.tsv.gz').set_index('cell_id');ev=read(M/'tables/expanded_classified_cells.tsv.gz').set_index('cell_id');su=read(M/'tables/supported_classified_cells.tsv.gz').set_index('cell_id')
 assert len(calls)==8174 and calls.index.equals(ev.index) and calls.index.is_unique
 probs=read(O/'tables/HLCA_finest_label_weights.tsv.gz').set_index('cell_id')
 assert probs.index.equals(calls.index) and np.isfinite(probs.to_numpy()).all() and np.allclose(probs.sum(1),1)
 assert probs.stack().between(0,1+1e-10).all()
 cov=read(O/'tables/HLCA_query_gene_coverage.tsv')
 assert cov['cells'].sum()==8174 and (cov.measured+cov.missing).eq(2000).all()
 z=np.load(O/'objects/HLCA_query_latent.npy');assert z.shape==(8174,30) and np.isfinite(z).all()
 rc=(O/'hlca.exit_code').read_text().strip();accepted_bookkeeping=False
 if rc!='0':
  log=(O/'logs/hlca_mapping.log').read_text()
  assert "dict() got multiple values for keyword argument 'at'" in log and (O/'objects/HLCA_query_model/model.pt').exists()
  accepted_bookkeeping=True
 summaries=[]
 for layer,labels in [('supported',su),('expanded',ev)]:
  tab=calls.loc[labels.index].copy();tab['class_id']=labels.class_id
  for (study,cls),g in tab.groupby(['dataset','class_id']):
   for ph in sorted(probs.columns):
    weights=probs.loc[g.index,ph];don=weights.groupby(g.donor_record).mean()
    summaries.append(dict(layer=layer,dataset=study,class_id=cls,phenotype=ph,cells=len(g),donor_records=len(don),predicted_cells=int(g.ann_finest_level_pred.eq(ph).sum()),unfiltered_label_fraction=float(g.ann_finest_level_pred.eq(ph).mean()),median_patient_weight=float(don.median()),mean_patient_weight=float(don.mean()),mean_prediction_uncertainty=float(g.ann_finest_level_uncertainty.mean())))
 pd.DataFrame(summaries).to_csv(O/'tables/HLCA_class_study_summary.tsv',sep='\t',index=False)
 acceptance={'at':now(),'status':'ACCEPTED_MAPPING_OUTPUTS_INTERPRETATION_PENDING','original_exit_code':rc,'original_status_write_error_only':accepted_bookkeeping,'original_exit_and_log_preserved':True,'cells':8174,'supported':1224,'reference_core_cells':qa['reference_core_cells'],'minimum_measured_model_genes':int(cov.measured.min()),'maximum_measured_model_genes':int(cov.measured.max()),'no_coverage_or_confidence_cell_filter':True,'cell_ids_and_probabilities_validated':True}
 writej(O/'delivery/HLCA_ACCEPTANCE_QA.json',acceptance);(O/'hlca.acceptance.exit_code').write_text('0')
elif kind=='sensitivity':
 S=O/'sensitivity_fixed_features';qa=json.loads((S/'delivery/SENSITIVITY_QA.json').read_text())
 assert (S/'sensitivity.exit_code').read_text().strip()=='0' and qa['complete']==20 and not qa['failures']
 summary=read(S/'tables/sensitivity_summary.tsv');assert len(summary)==20
 for layer,expected in [('supported',1224),('expanded',8174),('supported_myCAF',90),('expanded_myCAF',225)]:
  orig=ad.read_h5ad(M/f'objects/{layer}_integrated.h5ad',backed='r');ids=orig.obs_names;features=orig.var_names[orig.var.highly_variable.to_numpy(bool)]
  for r in summary.loc[summary.layer.eq(layer)].to_dict('records'):
   a=ad.read_h5ad(S/'objects'/f"{layer}__{r['scenario']}.h5ad",backed='r')
   assert a.n_obs==expected and a.obs_names.equals(ids) and a.var_names[a.var.highly_variable.to_numpy(bool)].equals(features)
   assert np.isfinite(a.obsm['X_umap']).all() and np.isfinite(a.obsm['X_harmony']).all();a.file.close()
  orig.file.close()
 acceptance={'at':now(),'status':'ACCEPTED_ALL20_MODELS','models':20,'exact_baseline_features_and_cells_retained':True,'finite_embeddings':True,'entry_threshold_changed':False}
 writej(S/'delivery/ACCEPTANCE_QA.json',acceptance);(S/'acceptance.exit_code').write_text('0')
else:raise ValueError(kind)

if kind=='hlca':
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 gong=read(O/'tables/Gong_human6_cell_scores.tsv.gz').set_index('cell_id')
 assert gong.index.equals(calls.index)
 joined=calls.join(gong[['raw_Gong_human6_mean','corrected_Gong_human6_mean']])
 joined.to_csv(O/'tables/HLCA_Gong_joint_cell_evidence.tsv.gz',sep='\t',index_label='cell_id')
 candidates=['Adventitial fibroblasts','Alveolar fibroblasts','Peribronchial fibroblasts','Myofibroblasts','Pericytes','Smooth muscle']
 displayrows=[];patientrows=[]
 for layer,lab in [('supported',su),('expanded',ev)]:
  v=joined.loc[lab.index].copy();v['class_id']=lab.class_id
  for key,g in v.groupby(['dataset','donor_record','site','class_id'],dropna=False):
   row=dict(zip(['dataset','donor_record','site','class_id'],key));row.update(layer=layer,cells=len(g),raw_Gong_human6_mean=float(g.raw_Gong_human6_mean.mean()),corrected_Gong_human6_mean=float(g.corrected_Gong_human6_mean.mean()),mean_uncertainty=float(g.ann_finest_level_uncertainty.mean()))
   for name in candidates:row[name]=float(probs.loc[g.index,name].mean()) if name in probs else 0.
   row['Other_reference_labels']=1-sum(row[name] for name in candidates);patientrows.append(row)
  # Main cohort shown without pooling biological N across possibly overlapping studies.
  g=v[v.dataset.eq('SCPCP000017')]
  for cls,sub in g.groupby('class_id'):
   row={'label':layer+' '+cls,'N':sub.donor_record.nunique()}
   for ph in candidates:
    vals=probs.loc[sub.index,ph] if ph in probs else pd.Series(0.,index=sub.index)
    row[ph]=vals.groupby(sub.donor_record).mean().mean()
   row['Other_reference_labels']=1-sum(row[ph] for ph in candidates);displayrows.append(row)
 pd.DataFrame(patientrows).to_csv(O/'tables/HLCA_Gong_patient_site_class_joint_evidence.tsv',sep='\t',index=False)
 show=pd.DataFrame(displayrows);cols=candidates+['Other_reference_labels']
 fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
 im=ax.imshow(show[cols].to_numpy(),aspect='auto',vmin=0,vmax=1,cmap='Blues')
 ax.set_yticks(range(len(show)),[r.label+' (N='+str(r.N)+')' for r in show.itertuples()]);ax.set_xticks(range(len(cols)),['AdvF','AlvF','Peribronchial','Myofibroblast','Pericyte','Smooth muscle','Other'],rotation=35,ha='right')
 for i in range(len(show)):
  for j in range(len(cols)):
   value=show[cols].iloc[i,j];ax.text(j,i,f'{value:.2f}',ha='center',va='center',fontsize=8,color='white' if value>.55 else 'black')
 ax.set_title('HLCA phenotype weights: SCPCP000017\nEqual patient means within class; no confidence cutoff')
 fig.colorbar(im,ax=ax,label='Mean reference-neighbor label weight');fig.savefig(O/'figures/HLCA_main_cohort_patient_weights.png',dpi=150);plt.close(fig)
 pairs=read(P/'phase31_customer_completion_20260926/resources/bone_lung_pairs_frozen.tsv')
 pairrows=[]
 for layer,lab in [('supported',su),('expanded',ev)]:
  v=joined.loc[lab.index]
  for row in pairs.to_dict('records'):
   for site,col in [('Bone','bone_samples'),('Lung','lung_samples')]:
    sub=v[v['sample'].isin(str(row[col]).split(';'))]
    for ph in candidates:
     pairrows.append(dict(layer=layer,donor_record=row['donor_record'],site=site,phenotype=ph,candidate_cells=len(sub),raw_label_cells=int(sub.ann_finest_level_pred.eq(ph).sum()),raw_label_fraction=float(sub.ann_finest_level_pred.eq(ph).mean()) if len(sub) else np.nan,mean_reference_weight=float(probs.loc[sub.index,ph].mean()) if len(sub) and ph in probs else np.nan,pair_type=row['pair_type'],confidence_cutoff_applied=False))
 pd.DataFrame(pairrows).to_csv(O/'tables/HLCA_paired_bone_lung_descriptive.tsv',sep='\t',index=False)
 (O/'HLCA_JOINT_EVIDENCE_COMPLETE').write_text(now())
 import subprocess
 subprocess.run([sys.executable,str(O/'code/report_delivery.py')],check=True)

print(json.dumps(acceptance,ensure_ascii=False))
