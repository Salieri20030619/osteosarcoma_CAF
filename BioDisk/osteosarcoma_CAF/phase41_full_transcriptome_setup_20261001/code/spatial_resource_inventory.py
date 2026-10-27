from pathlib import Path
import json,datetime,numpy as np,pandas as pd,anndata as ad,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');A=P/'phase41_full_transcriptome_setup_20261001';V=P/'phase35_external_validation_prep_20260927';D=P/'phase36_external_validation_20260927'
m=pd.read_csv(V/'tables/GSE293065_ready_objects.tsv',sep='\t',keep_default_na=False);clinical=pd.read_csv(D/'resources/sample_patient_manifest_corrected.tsv',sep='\t',keep_default_na=False)
rows=[]
for row in m.itertuples(index=False):
 a=ad.read_h5ad(row.object_path,backed='r');has='spatial' in a.obsm;xy=np.asarray(a.obsm['spatial']) if has else np.empty((0,2))
 rows.append({'sample':row.sample,'patient_id':row.patient_id,'spots':a.n_obs,'features':a.n_vars,'object_path':row.object_path,'object_bytes':Path(row.object_path).stat().st_size,'coordinate_array_present':has,'coordinate_shape':str(xy.shape),'coordinates_finite':bool(np.isfinite(xy).all()) if has else False,'spatial_uns_present':'spatial' in a.uns,'uns_keys':';'.join(a.uns.keys()),'obs_columns':';'.join(a.obs.columns),'obsm_keys':';'.join(a.obsm.keys()),'histology_image_or_tumor_boundary_accepted':False,'image_boundary_note':'inventory only; inspect actual source image/mask before boundary endpoint','clinical_source':str(D/'resources/sample_patient_manifest_corrected.tsv')})
 a.file.close()
frame=pd.DataFrame(rows);frame.to_csv(A/'tables/spatial_resource_inventory.tsv',sep='\t',index=False)
assert len(frame)==26 and frame.patient_id.nunique()==19 and frame.spots.sum()==53309
(A/'delivery/SPATIAL_RESOURCE_INVENTORY_QA.json').write_text(json.dumps({'at':datetime.datetime.now().astimezone().isoformat(),'samples':26,'patients':19,'spots':53309,'objects_checked':True,'coordinates_present':int(frame.coordinate_array_present.sum()),'coordinates_finite':int(frame.coordinates_finite.sum()),'geometry_or_immune_exclusion_not_yet_tested':True,'clinical_source':str(D/'resources/sample_patient_manifest_corrected.tsv')},indent=2))
# Formalize successful repaired stages, leaving their failed attempts archived.
B=P/'phase42_F3F4_identity_states_20261001';qa=json.loads((B/'delivery/IDENTITY_STATE_QA.json').read_text());assert (B/'identity_states.exit_code').read_text().strip()=='0' and (B/'COMPLETE').read_text().strip()=='SUCCESS'
qa.update(visual_review='all three internal UMAP PNGs inspected; axes/legends readable',source_and_donor_dominance='retained; state separation does not establish independent biological subtype',reviewed_figures=['supported_F3_internal_UMAP.png','expanded_F3_internal_UMAP.png','expanded_F4_internal_UMAP.png'])
(B/'delivery/IDENTITY_STATE_ACCEPTANCE.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print('SPATIAL_RESOURCES_INVENTORIED',frame[['sample','spatial_uns_present','coordinate_array_present']].to_string(index=False));print('STAGE_B_ACCEPTED')