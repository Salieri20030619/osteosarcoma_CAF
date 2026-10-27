from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_cnv_20260918');O=P/'result_review_20260919';T=O/'tables'
samples=pd.read_csv(P/'review/input_summary.tsv',sep='\t');allcells=[];checks=[];params=[];hashes=[]
def read(path):
 hashes.append({'path':str(path.relative_to(P)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 return pd.read_csv(path,sep='\t',dtype={'cluster':str})
for s in samples[samples.cnv_eligible].to_dict('records'):
 sm=s['sample'];raw=read(P/'inputs'/sm/'cells.tsv');a=read(P/'infercnv'/sm/'cell_results.tsv');b=read(P/'copykat'/sm/'cell_results.tsv')
 assert len(raw)==len(a)==len(b)==s['run_cells']
 for c in ['tool_id','cell_id','sample','patient_id','cluster','lineage_phase3','cnv_role','n_counts','n_genes']:
  assert raw[c].tolist()==a[c].tolist()==b[c].tolist(),(sm,c)
 assert not raw.tool_id.duplicated().any() and not raw.cell_id.duplicated().any()
 assert set(a.call)<=set(['low_signal','elevated_signal','ambiguous_signal','not_evaluable'])
 assert set(b.call)<=set(['diploid','aneuploid','not_evaluable'])
 for m,z in [('infercnv',a),('copykat',b)]:
  d=P/m/sm;pa=json.loads((d/'parameters.json').read_text());params.append(dict(pa))
  assert (d/'COMPLETE').exists() and (d/'job.exit_code').read_text().strip()=='0'
  assert pa['input_cells']==len(raw) and pa['reference_cells']==raw.cnv_role.str.startswith('ref_').sum() and pa['heldout_cells']==raw.cnv_role.str.startswith('holdout_').sum()
  assert np.isfinite(z.loc[z.call!='not_evaluable','signal_strength']).all()
 hold=a[a.cnv_role.str.startswith('holdout_')].signal_strength.to_numpy()
 # Reproduce R quantile(type=7) and mad(constant=1); parameters.json rounds numeric output.
 cutoff=max(.05,np.quantile(hold,.99),np.median(hold)+5*max(np.median(abs(hold-np.median(hold))),.005))
 assert np.allclose(a.signal_cutoff,cutoff,rtol=1e-10,atol=1e-12)
 ch=read(P/'infercnv'/sm/'chromosome_log2_medians.tsv').set_index('chromosome')
 ch=ch[a.tool_id].to_numpy(dtype=float);broad=(abs(ch)>=.1).sum(0)>=2; broad|=abs(ch).max(0)>=.25
 reproduced=np.where((a.signal_strength>cutoff)&broad,'elevated_signal',np.where(a.signal_strength<=cutoff,'low_signal','ambiguous_signal'))
 assert (reproduced==a.call).all()
 assert ((abs(ch)>=.1).sum(0)==a.altered_chromosomes).all()
 pred=read(P/'copykat'/sm/(sm+'_copykat_prediction.txt'))
 assert not pred['cell.names'].duplicated().any()
 original=pred.set_index('cell.names')['copykat.pred'].reindex(b.tool_id).fillna('missing')
 mapped=original.where(original.isin(['diploid','aneuploid']),'not_evaluable')
 assert mapped.tolist()==b.call.tolist()
 z=raw[['tool_id','cell_id','sample','patient_id','cluster','lineage_phase3','cnv_role','n_counts','n_genes','pctUMI_Fib','pctUMI_Osteo','pctUMI_Mural','pctUMI_Mineral','pctUMI_CTA','det_Fib_review','det_Osteo_review','det_Mural_review','CTA_context']].copy()
 z['dataset']=s['dataset'];z['site']=s['site'];z['infer_call']=a.call;z['copy_call']=b.call;z['copy_original_call']=original.to_numpy();z['infer_score']=a.signal_strength;z['copy_score']=b.signal_strength
 z['infer_cutoff']=cutoff;z['altered_chromosomes']=a.altered_chromosomes
 z['role_group']=np.where(z.cnv_role.str.startswith('ref_'),'training',np.where(z.cnv_role.str.startswith('holdout_'),'holdout','target'))
 z['joint']=np.select([(z.infer_call=='low_signal')&(z.copy_call=='diploid'),(z.infer_call=='elevated_signal')&(z.copy_call=='aneuploid'),(z.infer_call=='low_signal')&(z.copy_call=='aneuploid'),(z.infer_call=='elevated_signal')&(z.copy_call=='diploid'),z.copy_call=='not_evaluable',z.infer_call=='ambiguous_signal'],['both_low','both_high','infer_low_copy_high','infer_high_copy_low','copy_not_evaluable','infer_ambiguous'],default='other')
 allcells.append(z);checks.append(dict(sample=sm,cells=len(z),identity_match=True,threshold_reproduced=True,infer_calls_reproduced=True,copy_calls_reproduced=True,threshold=cutoff))
 print('CHECKED',sm,flush=True)
z=pd.concat(allcells,ignore_index=True);assert not z.cell_id.duplicated().any()
z.to_csv(T/'cell_review.tsv.gz',sep='\t',index=False)
pd.DataFrame(checks).to_csv(T/'integrity_checks.tsv',sep='\t',index=False)
pd.DataFrame(params).to_csv(T/'parameters.tsv',sep='\t',index=False)
pd.DataFrame(hashes).to_csv(T/'input_sha256.tsv',sep='\t',index=False)
for name,cols in [('sample_role_joint',['sample','dataset','site','role_group','joint']),('lineage_joint',['sample','dataset','site','lineage_phase3','joint']),('cluster_joint',['sample','dataset','site','cluster','lineage_phase3','joint']),('full_confusion',['sample','role_group','infer_call','copy_call']),('original_copy_calls',['sample','role_group','copy_original_call'])]:
 z.groupby(cols,dropna=False).size().rename('cells').reset_index().to_csv(T/(name+'.tsv'),sep='\t',index=False)
rows=[]
for (sm,role),g in z.groupby(['sample','cnv_role']):
 for m,lo,hi in [('infer','low_signal','elevated_signal'),('copy','diploid','aneuploid')]:
  c=g[m+'_call'];rows.append(dict(sample=sm,role=role,method=m,n=len(g),low=int((c==lo).sum()),high=int((c==hi).sum()),not_evaluable=int((c=='not_evaluable').sum()),ambiguous=int((c=='ambiguous_signal').sum()),high_fraction_all=float((c==hi).mean()),median_score=g[m+'_score'].median()))
pd.DataFrame(rows).to_csv(T/'calibration_by_role.tsv',sep='\t',index=False)
rows=[]
for (sm,lin,j),g in z.groupby(['sample','lineage_phase3','joint']):
 r=dict(sample=sm,lineage=lin,joint=j,n=len(g))
 for c in ['n_counts','n_genes','pctUMI_Fib','pctUMI_Osteo','pctUMI_Mural','pctUMI_Mineral','pctUMI_CTA','det_Fib_review','det_Osteo_review','det_Mural_review','infer_score','copy_score']:r[c+'_median']=g[c].median()
 rows.append(r)
pd.DataFrame(rows).to_csv(T/'expression_by_joint.tsv',sep='\t',index=False)
fib=z[z.lineage_phase3=='FibLike'];summary=[]
for sm,g in z.groupby('sample'):
 f=g[g.lineage_phase3=='FibLike'];h=g[g.role_group=='holdout'];t=g[g.role_group=='target'];binary=t[t.infer_call.isin(['low_signal','elevated_signal'])&t.copy_call.isin(['diploid','aneuploid'])]
 r=dict(sample=sm,dataset=g.dataset.iloc[0],site=g.site.iloc[0],cells=len(g),targets=len(t),holdout=len(h),holdout_infer_high=int((h.infer_call=='elevated_signal').sum()),holdout_copy_high=int((h.copy_call=='aneuploid').sum()),holdout_copy_ne=int((h.copy_call=='not_evaluable').sum()),target_binary_n=len(binary),target_binary_agreement=float(binary.joint.isin(['both_low','both_high']).mean()),target_copy_ne=int((t.copy_call=='not_evaluable').sum()),fib_cells=len(f))
 for j in ['both_low','both_high','infer_low_copy_high','infer_high_copy_low','copy_not_evaluable','infer_ambiguous','other']:r['fib_'+j]=int((f.joint==j).sum())
 summary.append(r)
su=pd.DataFrame(summary);su.to_csv(T/'sample_summary.tsv',sep='\t',index=False)
print(su.to_string(index=False));print('ALL_JOINT',z.joint.value_counts().to_dict());print('FIB_JOINT',fib.joint.value_counts().to_dict());print('CALIBRATION',z[z.role_group=='holdout'].groupby(['infer_call','copy_call']).size().to_dict())
(O/'summary.json').write_text(json.dumps(dict(reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),samples=len(su),cells=len(z),fib_cells=len(fib),joint=z.joint.value_counts().to_dict(),fib_joint=fib.joint.value_counts().to_dict(),all_integrity_checks_passed=True,inference_tests_performed=False,labels_changed=False),indent=2))
print('RESULT_TABLE_REVIEW_COMPLETE')
