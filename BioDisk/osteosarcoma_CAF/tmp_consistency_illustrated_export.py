from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json,shutil,hashlib,zipfile
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'report_consistency_illustrated_20261005';(O/'data').mkdir(parents=True,exist_ok=True)
m=json.loads((B/'CAF_CURRENT_CELL_EVIDENCE.json').read_text());cells=pd.read_csv(m['cell_table'],sep='\t',low_memory=False).set_index('cell_id')
cols=['sample','dataset','donor_record','site','class_id','subgroup_id','lineage_stratum','lineage_confidence','F4_major_branch','working_name','expression_confidence','functional_interpretation_confidence','doublet_positive_any','UMI','mitochondrial_fraction']
t=cells[cols].copy();t.index.name='cell_id';t.to_csv(O/'data/cells.tsv.gz',sep='\t')
paths={'global':'phase30_two_layer_classification_20260926/objects/expanded_integrated.h5ad','F2':'phase30_two_layer_classification_20260926/objects/expanded_myCAF_integrated.h5ad','F3':'phase42_F3F4_identity_states_20261001/objects/expanded_F3_internal_shared_embedding.h5ad','F4':'phase42_F3F4_identity_states_20261001/objects/expanded_F4_internal_shared_embedding.h5ad'}
manifest=[]
for key,rel in paths.items():
 a=ad.read_h5ad(B/rel,backed='r');u=pd.DataFrame(np.asarray(a.obsm['X_umap']),index=a.obs_names,columns=['UMAP1','UMAP2'])
 assert u.index.is_unique and np.isfinite(u.to_numpy()).all()
 expected=cells.index if key=='global' else cells.index[cells.class_id.eq(key)]
 assert set(u.index)==set(expected),(key,len(u),len(expected))
 u.index.name='cell_id';u.to_csv(O/'data'/('umap_'+key+'.tsv.gz'),sep='\t')
 manifest.append({'asset':key,'source':str(B/rel),'cells':len(u),'coordinates':'existing X_umap unchanged'});a.file.close()
for name in ['F1F5_subgroup_summary.tsv','four_axis_internal_overlay.tsv','pathway_display.tsv','F4_branches_by_subgroup.tsv','subgroup_confidence.tsv']:
 shutil.copy2(B/'phase62_F1F5_internal_20261005/tables'/name,O/'data'/name)
shutil.copy2(B/'phase66_doublet_repair_DSC1_20261005/QA/ACCEPTANCE.json',O/'data/phase66_ACCEPTANCE.json')
for k in ['class_id','subgroup_id','dataset','donor_record','F4_major_branch']:
 print(k,t[k].value_counts().to_dict())
print('donor_records_by_class',t.groupby('class_id').donor_record.nunique().to_dict())
print('source_donor_counts',t.groupby('dataset').donor_record.nunique().to_dict())
(O/'data/SOURCE_MANIFEST.json').write_text(json.dumps({'assets':manifest,'metadata_manifest':m,'files':[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((O/'data').iterdir()) if p.is_file()]},ensure_ascii=False,indent=2))
with zipfile.ZipFile(O/'data_export.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in (O/'data').iterdir():z.write(p,p.name)
print('EXPORT',O/'data_export.zip')
