from pathlib import Path
import json,csv
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');out=r/'phase5_nonpdx_qc_20260916';p4=r/'phase4_download_20260914'
m=json.loads((p4/'download_manifest.json').read_text());c=list(csv.DictReader((p4/'clinical_source_audit.tsv').open(encoding='utf-8-sig'),delimiter='\t'));clinical={(x['dataset'],x['sample']):x for x in c}
validated=list(csv.DictReader((p4/'source_and_matrix_validation.tsv').open(),delimiter='\t'));v={(x['dataset'],x['sample']):x for x in validated}
rows=[]
for x in m:
 cr=clinical[x['dataset'],x['sample']];key=x['dataset']+'__'+x['sample']
 rows.append({'sample':key,'short':x['sample'],'dataset':x['dataset'],'patient_id':x['patient_id'],'matrix_dir':v[x['dataset'],x['sample']]['matrix_dir'],'site':cr['sampling_site_verified'],'group':cr['clinical_group_verified'],'disease_stage_raw':cr['disease_stage_raw'],'treatment_raw':cr['preoperative_treatment_raw'],'technology':'scRNA_author_matrix','library_mapping_status':'tissue_matrix; capture mapping unresolved','source_evidence':x['source_evidence_url'],'input_cells':v[x['dataset'],x['sample']]['input_barcodes']})
with (out/'manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
s=(r/'OS_CAF_phase2/code/phase2.py').read_text()
s=s.replace('from scipy.io import mmread','from scipy.io import mmread\nfrom threadpoolctl import threadpool_limits\nthreadpool_limits(limits=4)')
a=s.index('def read_matrix(');b=s.index('\ndef mad(',a)
s=s[:a]+'''def read_matrix(path):
 records=json.loads((pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase4_download_20260914/matrix_input_resolution_20260916.json')).read_text())
 selected={x['role']:pathlib.Path(x['selected']) for x in records if pathlib.Path(x['selected']).parent==path}
 assert set(selected)=={'matrix','features','barcodes'},path
 def readlines(p):
  op=gzip.open if p.suffix=='.gz' else open
  with op(p,'rt') as f:return [x.rstrip('\\r\\n').split('\\t') for x in f if x.strip()]
 features=readlines(selected['features']);bc=[x[0] for x in readlines(selected['barcodes'])]
 matrix=mmread(str(selected['matrix'])).tocsr()
 assert matrix.shape==(len(features),len(bc))
 assert np.isfinite(matrix.data).all() and (matrix.data>=0).all() and np.equal(matrix.data,np.floor(matrix.data)).all()
 assert len(bc)==len(set(bc))
 keep=np.array([len(f)<3 or f[2]=='Gene Expression' for f in features])
 f=[f for f,k in zip(features,keep) if k];ids=[x[0] for x in f];assert len(ids)==len(set(ids))
 x=matrix[keep].T.tocsr().astype(np.float32);x.eliminate_zeros()
 a=ad.AnnData(x,obs=pd.DataFrame(index=bc),var=pd.DataFrame({'symbol':[x[1] for x in f]},index=ids))
 return a,int((~keep).sum())
''' +s[b:]
a=s.index(" url=row['source_url']");b=s.index(" log('SAMPLE_START'",a)
s=s[:a]+" path=pathlib.Path(row['matrix_dir']);assert path.is_dir()\n"+s[b:]
s=s.replace("a.obs['site']='Lung';a.obs['group']='B2'","a.obs['site']=row['site'];a.obs['group']=row['group']\n for key in ['dataset','disease_stage_raw','treatment_raw','library_mapping_status','source_evidence']:a.obs[key]=row[key]")
s=s.replace("assert q.n_obs>=100,'Too few cells after QC'","assert q.n_obs>=30,('Too few cells after QC for clustering; preserve and review',gsm,q.n_obs)")
a=s.index(' sc.pp.scrublet(q,');b=s.index(" allobs['doublet_score']",a)
s=s[:a]+''' threshold=float('nan');doublet_status='not_evaluated_low_cell_count'
 q.obs['doublet_score']=np.nan;q.obs['predicted_doublet']=pd.array([pd.NA]*q.n_obs,dtype='boolean')
 q.obs['doublet_exclude']=False
 if q.n_obs>=50:
  try:
   sc.pp.scrublet(q,expected_doublet_rate=args.doublet_rate,random_state=SEED,n_prin_comps=min(30,q.n_obs-2),use_approx_neighbors=False,verbose=False)
   threshold=float(q.uns['scrublet'].get('threshold',float('nan')))
   rate=float(q.obs.predicted_doublet.mean())
   valid=np.isfinite(threshold) and not q.obs.predicted_doublet.isna().any()
   doublet_status='automatic_threshold_pass' if valid and q.n_obs>=100 and 0.001<=rate<=0.25 else 'automatic_threshold_review_required'
   if doublet_status=='automatic_threshold_pass':q.obs['doublet_exclude']=q.obs.predicted_doublet.astype(bool)
  except (ValueError,RuntimeError,KeyError) as exc:
   doublet_status='scrublet_not_resolved_'+type(exc).__name__
   (out/'scrublet_issue.txt').write_text(repr(exc))
 q.obs['doublet_assessment_status']=doublet_status
 allobs['doublet_assessment_status']=doublet_status
 allobs['doublet_exclude']=q.obs.doublet_exclude.reindex(allobs.index).fillna(False).astype(bool)
 thresholds.update({'scrublet_threshold':threshold,'doublet_assessment_status':doublet_status,'threshold_review_policy':'review if <100 QC cells, missing threshold, predicted fraction <0.001 or >0.25; no fixed-fraction deletion','capture_mapping':'author tissue matrix; individual capture mapping unresolved'})
 dump(out/'qc_thresholds.json',thresholds)
''' +s[b:]
s=s.replace('a=q[~q.obs.predicted_doublet.values].copy();assert a.n_obs>=100','a=q[~q.obs.doublet_exclude.values].copy();assert a.n_obs>=30')
s=s.replace("axes[2].axvline(threshold,color='#B22222')","axes[2].axvline(threshold if np.isfinite(threshold) else 0,color='#B22222')")
s=s.replace("'doublets':int(q.obs.predicted_doublet.sum())","'doublets':int(q.obs.doublet_exclude.sum()),'predicted_doublets':int(q.obs.predicted_doublet.sum()),'doublet_status':doublet_status")
s=s.replace("'input_cells':n0","'dataset':row['dataset'],'patient_id':row['patient_id'],'site':row['site'],'clinical_group':row['group'],'input_cells':n0")
s=s.replace('PHASE2_COMPLETE','PHASE5_QC_CLUSTERING_COMPLETE')
s=s.replace('"""Per-sample RNA QC, Scrublet, Leiden and provisional lineage evidence."""','"""SYNTAX_CHECK: python -m py_compile phase5.py\nRUN_COMMAND: bash run_qc.sh\nEXPECTED_SUCCESS_MARKER: PHASE5_QC_CLUSTERING_COMPLETE\nAdapted from existing phase2.py; original outputs are preserved.\n"""')
(out/'code/phase5.py').write_text(s)
print('PREPARED',len(rows),sum(int(x['input_cells']) for x in rows))
