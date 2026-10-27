from pathlib import Path
import json,hashlib
import pandas as pd
import anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012'
for sub in ['tables','resources','QA','delivery','figures']:(O/sub).mkdir(exist_ok=True)
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
manifest=read(P/'phase26_client_reanalysis_20260922/tables/analysis_source_manifest.tsv',keep_default_na=False)
final=read(P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz',keep_default_na=False)
cohorts=['GSE152048','SCPCP000018','GSE318840'];rows=[];reasons=[];refs=[];later=[]
for r in manifest[manifest.dataset.isin(cohorts)].itertuples():
    sm=r.sample;f=P/'phase27_identity_review_20260926/tables'/f'{sm}_review.tsv.gz';d=read(f,keep_default_na=False)
    fib=d.provisional_lineage.eq('Fibroblast');low=d.CNV_evidence.isin(['single_method_no_alteration','concordant_no_alteration']);db=truth(d.known_doublet)
    reconstructed=fib&low&~db;original=truth(d.CAF_candidate)
    retained=final[final['sample'].eq(sm)]
    row=dict(dataset=r.dataset,sample=sm,site=r.site,input_cells=len(d),provisional_fibroblast_cells=int(fib.sum()),fibroblast_known_doublet=int((fib&db).sum()),fibroblast_concordant_alteration=int((fib&~db&d.CNV_evidence.eq('concordant_alteration')).sum()),fibroblast_single_alteration=int((fib&~db&d.CNV_evidence.eq('single_method_alteration')).sum()),fibroblast_method_conflict=int((fib&~db&d.CNV_evidence.eq('method_conflict')).sum()),fibroblast_unassessed=int((fib&~db&d.CNV_evidence.eq('unassessed')).sum()),stage27_CAF_candidates=int(original.sum()),final_candidate_cells=len(retained),after_stage27_removed=int(original.sum())-len(retained),gate_mismatches=int((reconstructed!=original).sum()),cell_ids_unique=d.cell_id.is_unique,final_ids_subset=bool(set(retained.cell_id)<=set(d.loc[original,'cell_id'])))
    rows.append(row)
    c=d.groupby(['provisional_lineage','copykat_call','infercnv_call','CNV_evidence','review_class']).size().reset_index(name='cells');c.insert(0,'sample',sm);c.insert(0,'dataset',r.dataset);reasons.append(c)
    ref=d[truth(d.reference)]
    rc=ref.groupby(['copykat_call','infercnv_call','CNV_evidence']).size().reset_index(name='reference_cells');rc.insert(0,'sample',sm);refs.append(rc)
    for phase in ['phase26_client_reanalysis_20260922/retry_once_20260924/tables','phase26_client_reanalysis_20260922/resource_recovery_SCPCS000525_20260925']:
        q=P/phase
        if sm in ['SCPCS000525','GSE152048__BC17','SCPCS000531','SCPCS000532']:
            for ef in q.glob(f'**/*{sm}*identity*'):
                if ef.is_file():later.append(dict(sample=sm,path=str(ef),size=ef.stat().st_size))
save(pd.DataFrame(rows),'sample_retention_waterfall.tsv');save(pd.concat(reasons,ignore_index=True),'lineage_CNV_exclusion_reasons.tsv');save(pd.concat(refs,ignore_index=True),'normal_reference_CNV_calls.tsv');save(pd.DataFrame(later),'later_evidence_files.tsv')
print('RETENTION\n'+pd.DataFrame(rows).to_string(index=False))
print('SAMPLE REFERENCES\n'+pd.concat(refs,ignore_index=True).to_string(index=False))
print('LATER_EVIDENCE',later)
source=P/'phase25_complete_nonPDX_assets_20260921/files/SCPCP000018/SCPCP000018_single-cell/single-cell_metadata.tsv'
md=read(source,keep_default_na=False);print('\nPROVIDER METADATA COLUMNS',md.columns.tolist());print(md[md.astype(str).apply(lambda c:c.str.contains('SARC0105|SARC0103|SARC0129')).any(axis=1)].to_string(index=False));save(md,'official_SCPCP000018_metadata.tsv')
for sm,lib in [('SCPCS000528','SCPCL000860'),('SCPCS000530','SCPCL000862')]:
    f=P/f'phase25_complete_nonPDX_assets_20260921/files/SCPCP000018/SCPCP000018_single-cell/{sm}/{lib}_processed_rna.h5ad'
    a=ad.read_h5ad(f,backed='r');print('\nOFFICIAL OBJECT',sm,'obs_columns',a.obs.columns.tolist(),'uns_keys',list(a.uns.keys()))
    for col in a.obs.columns:
        if any(k in col.lower() for k in ['tissue','site','diagnos','timing','participant']):print(col,a.obs[col].astype(str).unique().tolist()[:8])
    for k,v in a.uns.items():
        if isinstance(v,dict) and any(s in k.lower() for s in ['sample','metadata','project']):
            print('metadata_dict_keys',k,list(v))
    a.file.close()
print('\nREFERENCE_RERUN_FILES')
for sub in ['tables','resources','delivery']:
    p=P/'phase34_reference_rerun_20260927'/sub
    if p.exists():print(sub,[f.name for f in p.iterdir()][:35])
dump=dict(gate_mismatch_cells=sum(r['gate_mismatches'] for r in rows),all_final_ids_subset=all(r['final_ids_subset'] for r in rows),audited_samples=len(rows),all_ids_unique=all(r['cell_ids_unique'] for r in rows))
(O/'QA/RETENTION_CHECKS.json').write_text(json.dumps(dump,indent=2)+'\n')
