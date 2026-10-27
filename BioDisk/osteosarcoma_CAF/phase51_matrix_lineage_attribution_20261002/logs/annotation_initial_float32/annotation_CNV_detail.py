from pathlib import Path
import os,sys,json,datetime,hashlib
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';C=P/'phase43_full_transcriptome_discovery_20261001';O=P/'phase51_matrix_lineage_attribution_20261002'
m=pd.read_csv(A/'tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False);m=m[m.dataset.eq('SCPCP000017')]
cw=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);symbols=dict(zip(cw.feature,cw.canonical_symbol));quality=dict(zip(cw.feature,~cw.mapping_method.eq('HGNC_Ensembl_approved')))
a=ad.read_h5ad(A/'objects/all_candidate_full_genes.h5ad',backed='r');ob=a.obs.copy();a.file.close();allfrozen=set(ob.index)
annotation=(O/'resources/receiver_annotation_frozen.py').read_text();flag=lambda x:x.astype(str).str.lower().isin(['true','1'])
genes=['COL1A1','LOX','MMP14'];partition=['expanded__F'+str(i) for i in range(1,6)]+['T_NK','Myeloid','Endothelial','Mural','Fibroblast_other','CNV_altered_unresolved','Other_nonfrozen_unresolved']
rows=[];checks=[]
for row in m.itertuples(index=False):
 f=P/'phase26_client_reanalysis_20260922/inputs'/row.sample/'full_called_counts.h5ad';a=ad.read_h5ad(f);x=a.X.tocsr();meta=a.obs.copy()
 cells=pd.read_csv(f.parent/'cells.tsv',sep='\t',keep_default_na=False).set_index('cell_id');nonfrozen=~meta.index.isin(allfrozen);exec(annotation)
 here=ob[ob['sample'].astype(str).eq(row.sample)]
 for cl in ['F1','F2','F3','F4','F5']:mask['expanded__'+cl]=meta.index.isin(here.index[here.expanded_class_id.astype(str).eq(cl)])
 labels=np.full(a.n_obs,'',dtype=object)
 for group in partition:labels[mask[group]]=group
 assert np.all(np.column_stack([mask[g] for g in partition]).sum(1)==1)
 pref={}
 for j,en in enumerate(a.var_names.astype(str)):
  g=symbols.get(en,'');k=(quality.get(en,1),en)
  if g in genes and (g not in pref or k<pref[g][0]):pref[g]=(k,j)
 assert set(pref)==set(genes)
 raw=x[:,[pref[g][1] for g in genes]].toarray()
 bucket=provider.astype(str).copy();unknown=bucket.str.strip().str.lower().isin(['unknown','unannotated','','nan','na'])
 hknown=~hist.str.strip().str.lower().isin(['unknown','unannotated','unresolved','','nan','na'])
 bucket.loc[unknown&hknown]='historical:'+hist.loc[unknown&hknown]
 bucket.loc[unknown&~hknown]='unassigned_annotation'
 detail=pd.DataFrame({'group':labels,'annotation':bucket.to_numpy(),'CNV_evidence':evidence.CNV_evidence.reindex(meta.index).astype(str).to_numpy(),'reference':reference,'cells':1,'library':pd.to_numeric(cells.loc[a.obs_names,'n_counts']).to_numpy()})
 for j,g in enumerate(genes):detail[g]=raw[:,j]
 detail=detail.groupby(['group','annotation','CNV_evidence','reference'],dropna=False)[['cells','library']+genes].sum().reset_index()
 detail['sample']=row.sample;detail['dataset']=row.dataset;detail['patient_id']=row.donor_record;detail['site']=row.site;rows.append(detail)
 checks.append({'sample':row.sample,'cells':int(detail.cells.sum()),'groups':len(detail),'all_counts_match':bool(np.array_equal(detail[genes].sum().to_numpy(),raw.sum(0)))})
 assert checks[-1]['cells']==a.n_obs and checks[-1]['all_counts_match']
 print('DETAILED_SAMPLE',row.sample,flush=True)
d=pd.concat(rows,ignore_index=True);d.to_csv(O/'tables/main_annotation_CNV_sample_counts.tsv',sep='\t',index=False)
base=pd.read_csv(O/'tables/sample_group_expression.tsv.gz',sep='\t')
b=base[base.dataset.eq('SCPCP000017')&base.group_scheme.eq('exhaustive_partition')]
wide=b.pivot(index=['sample','group'],columns='gene',values='counts').sort_index()
got=d.groupby(['sample','group'])[genes].sum().reindex(wide.index,fill_value=0)
assert np.array_equal(got[genes].to_numpy(),wide[genes].to_numpy())
site=d.groupby(['patient_id','site','group','annotation','CNV_evidence','reference'],dropna=False)[['cells','library']+genes].sum().reset_index()
site.to_csv(O/'tables/main_annotation_CNV_patient_site_counts.tsv',sep='\t',index=False)
# Every observed category gets a true zero contribution for patient-sites where it has no captured cells.
cats=site[['group','annotation','CNV_evidence','reference']].drop_duplicates()
sites=site[['patient_id','site']].drop_duplicates()
grid=sites.merge(cats,how='cross').merge(site,on=['patient_id','site','group','annotation','CNV_evidence','reference'],how='left')
grid[['cells','library']+genes]=grid[['cells','library']+genes].fillna(0)
total=site.groupby(['patient_id','site'])[genes].sum().reset_index().rename(columns={g:g+'_total' for g in genes})
grid=grid.merge(total,on=['patient_id','site'])
long=[]
for g in genes:
 t=grid[['patient_id','site','group','annotation','CNV_evidence','reference']].copy();t['gene']=g
 t['share']=np.divide(grid[g],grid[g+'_total'],out=np.full(len(grid),np.nan),where=grid[g+'_total'].gt(0));long.append(t)
long=pd.concat(long,ignore_index=True)
patient=long.groupby(['patient_id','group','annotation','CNV_evidence','reference','gene'],dropna=False).share.mean().reset_index()
patient.to_csv(O/'tables/main_annotation_CNV_patient_shares.tsv',sep='\t',index=False)
summ=patient.groupby(['group','annotation','CNV_evidence','reference','gene'],dropna=False).agg(N_patients=('patient_id','nunique'),mean_total_gene_share=('share','mean')).reset_index()
summ.to_csv(O/'tables/main_annotation_CNV_summary.tsv',sep='\t',index=False)
assert np.allclose(patient.groupby(['patient_id','gene']).share.sum(),1)
# Category sums must match the original independently validated per-patient partition shares.
expected=pd.read_csv(O/'tables/patient_expression_and_shares.tsv',sep='\t',float_precision='round_trip')
expected=expected[expected.dataset.eq('SCPCP000017')&expected.group_scheme.eq('exhaustive_partition')]
joined=patient.groupby(['patient_id','group','gene']).share.sum().reset_index().merge(expected,on=['patient_id','group','gene'])
assert len(joined)==20*12*3
err=float(np.max(np.abs(joined.share-joined.gene_capture_share)));assert err<1e-12
qa={'samples':len(checks),'cells':sum(c['cells'] for c in checks),'patient_N':20,'all_sample_group_gene_counts_exact':True,'parent_patient_share_max_error':err,'annotation_and_CNV_are_existing_metadata_only':True,'new_genebased_labels':False,'new_tests':0,'after_primary_descriptive_followup':True}
(O/'delivery/ANNOTATION_CNV_QA.json').write_text(json.dumps(qa,indent=2));(O/'ANNOTATION_CNV_COMPLETE').write_text('SUCCESS');print(json.dumps(qa),flush=True)
