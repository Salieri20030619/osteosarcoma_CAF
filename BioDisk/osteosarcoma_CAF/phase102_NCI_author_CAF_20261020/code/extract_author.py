import pathlib,json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
P=R/'phase102_NCI_author_CAF_20261020';Q=R/'phase95_NCI_expanded_reanalysis_20261016'
for s in ['code','tables','resources','logs','report']: (P/s).mkdir(parents=True,exist_ok=True)
genes=json.loads((R/'phase96_all14_literature_20261016/project_registry.json').read_text())
aliases=json.loads((R/'phase97_normal_background_20261020/resources/GENE_ALIAS_MAP.json').read_text())['aliases']
run=json.loads((Q/'COMPUTATION_FINISHED.json').read_text())['samples']
rows=[];summary=[];cells=[];hashes=[]
for rr in run:
 sm=rr['sample'];file=Q/'tables'/(sm+'_all_cell_evidence.tsv.gz')
 c=pd.read_csv(file,sep='\t',encoding='utf-8')
 mask=c.Ann_Level0.eq('Host') & c.Ann_Level2.eq('CAFs')
 assert np.array_equal(mask,c.author_CAF)
 z=c[mask].copy();pid='NCI_PATIENT_'+sm.rsplit('_',1)[-1].removesuffix('b');z['patient_id']=pid
 cells.append(z[['cell_id','sample','patient_id','Ann_Level0','Ann_Level2','full_scDblFinder_class','native_copykat_call','mechanical_expanded_candidate']])
 summary.append(dict(sample=sm,patient_id=pid,author_CAF=int(mask.sum()),author_CAF_singlets=int(z.full_scDblFinder_class.eq('singlet').sum()),old_reference_pass=rr['reference_pass']))
 hashes.append(dict(path=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
 if z.empty:continue
 a=ad.read_h5ad(Q/'inputs'/sm/'available_counts.h5ad',backed='r');ix=a.obs_names.get_indexer(z.cell_id);assert (ix>=0).all()
 b=a[ix,:].to_memory();a.file.close()
 den=np.asarray(b.X.sum(1)).ravel();assert np.array_equal(den,z.n_counts.to_numpy())
 for g in genes:
  hit=np.flatnonzero((b.var_names==g['feature_ids']) | b.var.symbol.isin(aliases[g['gene']]) | b.var_names.isin(aliases[g['gene']]))
  assert len(hit)==1,(sm,g['gene'],len(hit))
  for st,m in [('author_all',np.ones(len(z),bool)),('author_singlets',z.full_scDblFinder_class.eq('singlet').to_numpy())]:
   rows.append(dict(sample=sm,patient_id=pid,stratum=st,gene=g['gene'],feature=g['feature_ids'],cells=int(m.sum()),raw_count=float(b.X[m][:,hit].sum()),library_UMIs=float(den[m].sum())))
 print(sm,len(z),flush=True)
lib=pd.DataFrame(rows)
agg=lib.groupby(['patient_id','stratum','gene','feature'],as_index=False)[['cells','raw_count','library_UMIs']].sum()
agg=agg[agg.cells.gt(0)].copy();agg['expression']=np.log2(agg.raw_count/agg.library_UMIs*1e4+.1)
lib.to_csv(P/'tables/author_CAF_library_counts.tsv',sep='\t',index=False,encoding='utf-8')
agg.to_csv(P/'tables/author_CAF_patient_expression.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(summary).to_csv(P/'tables/author_CAF_library_coverage.tsv',sep='\t',index=False,encoding='utf-8')
pd.concat(cells).to_csv(P/'tables/author_CAF_cell_registry.tsv.gz',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(hashes).to_csv(P/'resources/input_hashes.tsv',sep='\t',index=False,encoding='utf-8')
d=dict(libraries=19,source_patients=17,author_CAF_total=int(sum(x['author_CAF'] for x in summary)),author_singlets_total=int(sum(x['author_CAF_singlets'] for x in summary)),patients_by_stratum=agg.groupby('stratum').patient_id.nunique().to_dict(),selection='Author Host + CAFs only; no CopyKAT, reference-pass or fibroblast-score gate. Singlets separate sensitivity.',patient_grouping='User confirmed paired suffixes 16178/b and 16183/b',normal_reference_threshold=20)
(P/'resources/EXTRACTION_COMPLETE.json').write_text(json.dumps(d,indent=2),encoding='utf-8');print(json.dumps(d),flush=True)
