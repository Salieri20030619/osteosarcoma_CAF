from pathlib import Path
import json,csv,hashlib
import pandas as pd,anndata as ad
from range_zip import ROOT,TARGETS
assert int((ROOT/'exit_code').read_text())==0 and (ROOT/'DOWNLOAD_COMPLETE').exists()
state=json.loads((ROOT/'STATUS.json').read_text());files=json.loads((ROOT/'verified_files.json').read_text());inv=json.loads((ROOT/'matrix_inventory.json').read_text())
assert state['status']=='DOWNLOAD_COMPLETE_VERIFIED' and len(inv)==27
rows=[]
for project,samples in TARGETS.items():
    meta_path=ROOT/'files'/project/(project+'_single-cell')/'single-cell_metadata.tsv'
    meta=pd.read_csv(meta_path,sep='\t',keep_default_na=False)
    for sample in samples:
        selected=meta[meta.scpca_sample_id==sample];assert len(selected)==1
        m=selected.iloc[0].to_dict();assert str(m['is_xenograft']).lower()=='false'
        subset={}
        for kind in ['unfiltered','filtered','processed']:
            f=next(r for r in files if r['project']==project and sample in r['archive_member'] and '_'+kind+'_rna.h5ad' in r['archive_member'])
            a=ad.read_h5ad(f['path'],backed='r');subset[kind]=(set(a.obs_names),set(a.var_names),a.n_obs,a.n_vars);a.file.close()
        assert subset['processed'][0]<=subset['filtered'][0]<=subset['unfiltered'][0]
        assert subset['processed'][1]<=subset['filtered'][1]<=subset['unfiltered'][1]
        assert subset['processed'][2]==int(m['processed_cells'])
        assert subset['filtered'][2]==int(m['filtered_cell_count'])
        rows.append({**{k:m.get(k,'FIELD_ABSENT') for k in ['scpca_project_id','scpca_sample_id','scpca_library_id','participant_id','tissue_location','primary_or_metastasis','disease_timing','treatment_status','is_xenograft']},'unfiltered_barcodes':subset['unfiltered'][2],'filtered_cells':subset['filtered'][2],'processed_cells':subset['processed'][2],'genes_processed':subset['processed'][3],'barcode_subset_check':'PASS','official_metadata_count_match':'PASS'})
pd.DataFrame(rows).to_csv(ROOT/'sample_download_acceptance.tsv',sep='\t',index=False)
qa={'all_checks_passed':True,'samples':9,'expression_objects':27,'selected_files':len(files),'ZIP_CRC_verified_by_downloader':True,'SHA256_recorded':True,'AnnData_readable':True,'processed_and_filtered_match_official_metadata':True,'processed_subset_filtered_subset_unfiltered':True,'xenografts':0,'biological_analysis_started':False}
(ROOT/'DOWNLOAD_QA.json').write_text(json.dumps(qa,indent=2));print(pd.DataFrame(rows).to_string(index=False));print(json.dumps(qa))
