from pathlib import Path
import json,datetime,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase25_complete_nonPDX_assets_20260921'
state=json.loads((O/'STATUS.json').read_text());assert state['status']=='DOWNLOAD_COMPLETE_VERIFIED'
new=json.loads((O/'matrix_inventory.json').read_text());assert len(new)==165
old=pd.read_csv(P/'phase24_asset_inventory_20260921/tables/ScPCA_90_library_assets.tsv',sep='\t',keep_default_na=False)
expected=old[old.is_xenograft.astype(str).str.lower().eq('false')&old.is_cell_line.astype(str).str.lower().eq('false')];assert len(expected)==67
rows=[]
for r in expected.to_dict('records'):
 sm=r['scpca_sample_id'];out={k:r[k] for k in ['scpca_project_id','scpca_sample_id','scpca_library_id','participant_id','tissue_location','disease_timing','treatment_status','is_xenograft','is_cell_line']}
 out['nonPDX_rechecked']=True
 for kind in ['filtered','processed','unfiltered']:
  matches=[v for v in new if v['sample']==sm and v['file'].endswith('_'+kind+'_rna.h5ad')]
  if matches:
   assert len(matches)==1;p=O/'files'/r['scpca_project_id']/matches[0]['file'];out[kind+'_cells']=matches[0]['cells'];out[kind+'_asset_source']='new_official_download_phase25'
  else:
   parts=r[kind+'_paths'].split(';');assert len(parts)==1 and parts[0];p=Path(parts[0]);out[kind+'_asset_source']='existing_phase15_previously_rechecked_phase24'
  assert p.is_file() and p.stat().st_size>0;out[kind]=str(p)
 rows.append(out)
a=pd.DataFrame(rows);assert a.scpca_sample_id.is_unique and len(a)==67
a.to_csv(O/'all67_nonPDX_asset_manifest.tsv',sep='\t',index=False)
q={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_sample_records':55,'new_expression_objects':165,'existing_sample_records':12,'total_nonPDX_sample_records':67,'total_expression_objects':201,'PDX_in_manifest':0,'cell_lines_in_manifest':0,'data_integrity':'downloaded range lengths, ZIP CRC32 and SHA256 verified by downloader; existing objects reused with previous evidence','source_recheck':'all new objects require is_xenograft=False, is_cell_line=False and matching sample_id','formal_analysis_started':False,'sample_counts_are_not_patient_N':True,'remaining_other_assets':['GSE237070 matrices','GSE217792 reference','HLCA reference implementation','Budhathoki source deduplication','GSE270231 full-called recovery']}
(O/'DOWNLOAD_QA.json').write_text(json.dumps(q,indent=2));(O/'ALL67_ASSETS_VERIFIED').write_text('ALL_67_NON_PDX_SAMPLE_ASSETS_VERIFIED\n');print(json.dumps(q,indent=2))
x=json.loads((P/'PROJECT_STATE.json').read_text());x['phase25_complete_nonPDX_assets_20260921']={'status':'ALL67_ASSETS_VERIFIED','manifest':str(O/'all67_nonPDX_asset_manifest.tsv'),'QA':q};(P/'PROJECT_STATE.json').write_text(json.dumps(x,ensure_ascii=False,indent=2))
