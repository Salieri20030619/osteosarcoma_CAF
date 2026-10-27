from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,traceback
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase83_upstream_and_direction_20261012';M=P/'phase82_allCAF_crossdataset_markers_20261012'
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False,na_rep='NA')
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
def digest(f):return hashlib.sha256(f.read_bytes()).hexdigest()

def main():
    main=read(M/'tables/primary_same_layer_crossdataset_markers.tsv');allstats=read(M/'tables/all_candidate_source_layer_design_statistics.tsv')
    evidence=allstats[(allstats.sensitivity=='all_cells')&(allstats.design=='independent')&allstats.symbol.isin(main.symbol)].copy()
    clinical=read(P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv')
    sources=sorted(clinical.dataset.unique());rows=[]
    for layer in ['supported','expanded']:
        for r in main.itertuples():
            expected=1 if r.direction=='lung_higher' else -1
            for ds in sources:
                z=evidence[(evidence.layer==layer)&(evidence.dataset==ds)&(evidence.symbol==r.symbol)]
                if len(z):
                    assert len(z)==1;q=z.iloc[0];effect=q.mean_lung_minus_bone
                    direction='same' if effect*expected>0 else 'opposite' if effect*expected<0 else 'zero' if pd.notna(effect) else 'not_evaluable'
                    rows.append(dict(layer=layer,symbol=r.symbol,dataset=ds,expected_direction=r.direction,bone_N=q.bone_N,lung_N=q.lung_N,effect=effect,P=q.P,FDR=q.FDR,direction=direction,statistical_evidence='nominal_same_direction' if q.P<.05 and direction=='same' else 'nominal_opposite_direction' if q.P<.05 and direction=='opposite' else 'tested_not_nominal' if pd.notna(q.P) else 'direction_only_insufficient_replication',missing_reason=''))
                else:
                    reg=read(P/'phase81_allCAF_bone_lung_20261012/tables/comparison_registry.tsv')
                    rr=reg[(reg.layer==layer)&(reg.dataset==ds)&(reg.sensitivity=='all_cells')&(reg.design=='independent')]
                    b=int(rr.bone_N.iloc[0]) if len(rr) else 0;l=int(rr.lung_N.iloc[0]) if len(rr) else 0
                    rows.append(dict(layer=layer,symbol=r.symbol,dataset=ds,expected_direction=r.direction,bone_N=b,lung_N=l,effect=np.nan,P=np.nan,FDR=np.nan,direction='not_evaluable',statistical_evidence='not_evaluable',missing_reason='one_site_has_no_candidate_patients' if not b or not l else 'feature_not_in_source_canonical_measured_mapping'))
    directions=pd.DataFrame(rows);save(directions,'all14_markers_all_sources_direction_evidence.tsv')
    summary=[]
    for r in main.itertuples():
        z=directions[(directions.layer=='expanded')&(directions.symbol==r.symbol)]
        row=dict(symbol=r.symbol,expected_direction=r.direction,same_sources=int(z.direction.eq('same').sum()),opposite_sources=int(z.direction.eq('opposite').sum()),zero_sources=int(z.direction.eq('zero').sum()),not_evaluable_sources=int(z.direction.eq('not_evaluable').sum()),same_source_list=';'.join(z.loc[z.direction.eq('same'),'dataset']),opposite_source_list=';'.join(z.loc[z.direction.eq('opposite'),'dataset']),zero_source_list=';'.join(z.loc[z.direction.eq('zero'),'dataset']),nominal_same_sources=int(z.statistical_evidence.eq('nominal_same_direction').sum()))
        for ds in ['SCPCP000017','SCPCP000023','SCPCP000018','GSE152048','GSE318840']:
            q=z[z.dataset.eq(ds)].iloc[0]
            for k in ['effect','P','FDR','bone_N','lung_N','direction']:row[ds+'__'+k]=q[k]
        summary.append(row)
    summary=pd.DataFrame(summary);save(summary,'MARKER_DIRECTION_SUMMARY.tsv')
    # Update late source evidence in a separately versioned provisional recovery view.
    base=P/'phase26_client_reanalysis_20260922/resource_recovery_SCPCS000525_20260925'
    qa=json.loads((base/'delivery/QA.json').read_text());assert qa['native_success'] and qa['native_exit_code']==0
    old=read(P/'phase27_identity_review_20260926/tables/SCPCS000525_review.tsv.gz',keep_default_na=False).set_index('cell_id')
    new=read(base/'tables/SCPCS000525_identity_evidence.tsv.gz',keep_default_na=False).set_index('cell_id').loc[old.index]
    keep=old.provisional_lineage.eq('Fibroblast')&~truth(old.known_doublet)&new.CNV_evidence.isin(['concordant_no_alteration','single_method_no_alteration'])
    recovered=old.loc[keep].copy();recovered['updated_CNV_evidence']=new.loc[keep,'CNV_evidence'];recovered['updated_copykat_call']=new.loc[keep,'copykat_call'];recovered['updated_infercnv_call']=new.loc[keep,'infercnv_call'];recovered['status']='recovered_by_original_entry_rule_pending_final_cluster_lineage_review';recovered['evidence_source']=str(base/'tables/SCPCS000525_identity_evidence.tsv.gz')
    assert len(recovered)==384
    save(recovered.reset_index(),'SCPCS000525_provisional_recovered_cells.tsv.gz')
    final=read(P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz');assert not set(recovered.index)&set(final.cell_id)
    source=P/'phase26_client_reanalysis_20260922/inputs/SCPCS000525/full_called_counts.h5ad';a=ad.read_h5ad(source,backed='r');ix=a.obs_names.get_indexer(recovered.index);assert (ix>=0).all()
    x=a.X[ix,:].tocsr();library=float(x.sum());assert np.allclose(np.asarray(x.sum(1)).ravel(),recovered.n_counts.astype(float))
    manifest=read(M/'tables/all_candidate_source_layer_design_statistics.tsv');s18=manifest[(manifest.layer=='expanded')&(manifest.dataset=='SCPCP000018')&(manifest.design=='independent')&(manifest.sensitivity=='all_cells')&manifest.symbol.isin(main.symbol)]
    pat=read(M/'tables/candidate_patient_expression.tsv');pat=pat[(pat.layer=='expanded')&(pat.dataset=='SCPCP000018')&(pat.design=='independent')&pat.symbol.isin(main.symbol)]
    det=read(M/'tables/candidate_detection_and_Fibro_patient_expression.tsv');det=det[(det.layer=='expanded')&(det.dataset=='SCPCP000018')&(det.subset=='all_candidates')]
    donor='SCPCP000018|SARC0088';sensitivity=[]
    for r in s18.itertuples():
        j=a.var_names.get_indexer([r.feature])[0];assert j>=0
        newcount=float(x[:,j].sum());oldrow=det[(det.symbol==r.symbol)&(det.donor_record==donor)&(det.site=='Lung')].iloc[0]
        assert oldrow.measured and oldrow.cells==447
        updated=np.log2((newcount+oldrow.raw_gene_UMIs)/(library+oldrow.library_UMIs)*1e4+.1)
        z=pat[pat.symbol.eq(r.symbol)].copy();z.loc[z.donor_record.eq(donor)&z.site.eq('Lung'),'pseudobulk_log2_CP10k_plus01']=updated
        effect=z.loc[z.site.eq('Lung'),'pseudobulk_log2_CP10k_plus01'].mean()-z.loc[z.site.eq('Bone'),'pseudobulk_log2_CP10k_plus01'].mean()
        sign=1 if main.set_index('symbol').loc[r.symbol,'direction']=='lung_higher' else -1
        sensitivity.append(dict(symbol=r.symbol,bone_N=1,lung_N=4,original_effect=r.mean_lung_minus_bone,recovery_effect=effect,original_direction='same' if r.mean_lung_minus_bone*sign>0 else 'opposite' if r.mean_lung_minus_bone*sign<0 else 'zero',recovery_direction='same' if effect*sign>0 else 'opposite' if effect*sign<0 else 'zero',new_patient_added=False,added_provisional_cells=384,patient_original_cells=447,patient_augmented_cells=831,added_gene_UMIs=newcount,added_library_UMIs=library,updated_patient_expression=updated,scope='provisional recovery sensitivity only; no inferred P or frozen-label change'))
    a.file.close();sensitivity=pd.DataFrame(sensitivity);save(sensitivity,'SCPCP000018_late_recovery_direction_sensitivity.tsv')
    # Quantify calibration warnings rather than reclassifying normal-reference labels as ground truth.
    refs=read(O/'tables/normal_reference_CNV_calls.tsv');refsum=[]
    for sm,g in refs.groupby('sample'):
        n=int(g.reference_cells.sum());alter=int(g.loc[g.infercnv_call.eq('alteration_called'),'reference_cells'].sum())
        refsum.append(dict(sample=sm,designated_reference_cells=n,infercnv_alteration_called=alter,infercnv_alteration_fraction=alter/n if n else np.nan,copykat_aneuploid=int(g.loc[g.copykat_call.eq('aneuploid'),'reference_cells'].sum()),interpretation='Reference calibration/annotation warning; not proof every excluded cell is malignant or normal'))
    save(pd.DataFrame(refsum),'designated_reference_CNV_calibration_summary.tsv')
    result=dict(status='analysis_complete_review_pending',completed_at=datetime.now(timezone.utc).isoformat(),retention_audit_samples=30,site_labels_changed=False,reason_no_site_change='Official metadata and processed-object field remain Arm/Leg/Thigh; primary osteosarcoma timing is confirmed but sampled tissue compartment not resolved',provisional_cells_recovered=384,independent_patients_added=0,main_candidate_pool_unchanged=True,new_inferential_P_values=0,late_recovery_marker_direction_changes=int((sensitivity.original_direction!=sensitivity.recovery_direction).sum()),main_marker_count=14,direction_evaluable_sources=['SCPCP000017','SCPCP000023','SCPCP000018','GSE152048','GSE318840'])
    (O/'ANALYSIS_COMPLETE.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2));print(summary[['symbol','same_sources','opposite_sources','zero_sources','same_source_list','opposite_source_list']].to_string(index=False));print('LATE_RECOVERY\n'+sensitivity[['symbol','original_effect','recovery_effect','original_direction','recovery_direction']].to_string(index=False))

if __name__=='__main__':
    try:main()
    except BaseException:
        (O/'logs').mkdir(exist_ok=True);(O/'logs/analysis.exit_code').write_text('1\n');traceback.print_exc();raise
    else:
        (O/'logs').mkdir(exist_ok=True);(O/'logs/analysis.exit_code').write_text('0\n')
