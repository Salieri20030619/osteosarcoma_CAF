from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');R=P/'phase53_de_novo_reannotation_20261002'
curdoc=json.loads((R/'resources/CURATION_LOCKED.json').read_text());cur=pd.DataFrame(curdoc['clusters']);assert len(cur)==307 and not cur.old_labels_consulted.any()
cur['curation_id']=cur.dataset+'__K'+cur.cluster.astype(str)
cur.to_csv(R/'tables/proposed_cluster_annotations.tsv',sep='\t',index=False)
maps=json.loads((R/'unblinding/identifier_maps.json').read_text())
manifest=pd.read_csv(R/'resources/COPYKAT_PROFILE_MANIFEST.tsv',sep='\t',keep_default_na=False).set_index('sample')
clin=pd.read_csv(P/'phase41_full_transcriptome_setup_20261001/tables/clinical_sample_registry.tsv',sep='\t',keep_default_na=False)
snap={x['sample']:x for x in json.loads((P/'phase27_identity_review_20260926/resources/SOURCE_SNAPSHOT.json').read_text())}
cnv=[];evledger=[]
for row in clin.itertuples(index=False):
 ef=Path(snap[row.sample]['evidence']);e=pd.read_csv(ef,sep='\t',keep_default_na=False)
 selected=manifest.loc[row.sample,'copykat_rds']
 if selected:
  calls=pd.read_csv(Path(selected).parent/'cell_results.tsv',sep='\t',keep_default_na=False).set_index('cell_id')['call']
  assert set(calls.index)==set(e.cell_id)
  e['copykat_call']=e.cell_id.map(calls)
 ca=e.copykat_call.eq('aneuploid');cn=e.copykat_call.eq('diploid');ia=e.infercnv_call.eq('alteration_called');inn=e.infercnv_call.eq('no_alteration_called')
 state=np.select([ca&ia,cn&inn,(ca&inn)|(cn&ia),ca|ia,cn|inn],['concordant_alteration','concordant_no_alteration','method_conflict','single_method_alteration','single_method_no_alteration'],default='unassessed')
 out=e[['cell_id','copykat_call','infercnv_call','infercnv_source']].copy();out['CNV_evidence_current']=state;out['sample']=row.sample;cnv.append(out)
 evledger.append({'sample':row.sample,'identity_evidence_file':str(ef),'identity_SHA256':hashlib.sha256(ef.read_bytes()).hexdigest(),'current_copykat_profile_source':selected})
cnv=pd.concat(cnv,ignore_index=True);assert cnv.cell_id.is_unique and len(cnv)==876648
parts=[];allrefs=[]
for ds in sorted(cur.dataset.unique()):
 d=R/'tables'/ds
 obs=pd.read_csv(d/'cell_structure.tsv.gz',sep='\t',dtype={'r05':str,'r10':str,'r20':str},keep_default_na=False)
 pm=pd.read_csv(d/'reference_profile_metadata.tsv',sep='\t',dtype={'cluster':str},keep_default_na=False)
 ra=pd.read_csv(d/'reference_annotations.tsv',sep='\t',keep_default_na=False);joined=ra.merge(pm,on='profile_id');allrefs.append(joined.assign(dataset=ds))
 for res in ['r05','r10','r20']:
  for ref in ['HPCA','BlueprintEncode']:
   tab=joined[(joined.kind=='global')&(joined.resolution==res)&(joined.reference==ref)].set_index('cluster')
   obs[res+'_'+ref+'_label']=obs[res].map(tab.label)
   obs[res+'_'+ref+'_pruned_label']=obs[res].map(tab.pruned_label)
 for ref in ['HPCA','BlueprintEncode']:
  tab=joined[(joined.kind=='patient')&(joined.resolution=='r10')&(joined.reference==ref)].copy()
  tab['join_id']=tab.patient_id+'|||'+tab.cluster;tab=tab.set_index('join_id')
  obs['patient_'+ref+'_label']=(obs.patient_id+'|||'+obs.r10).map(tab.label)
 review=cur[cur.dataset.eq(ds)].rename(columns={'cluster':'r10'})
 obs=obs.merge(review[['r10','curation_id','proposed_lineage','state','status']],on='r10',how='left',validate='many_to_one')
 assert obs.proposed_lineage.notna().all()
 obs['reference_resolution_disagreement']=obs.r10_HPCA_label.ne(obs.r20_HPCA_label)|obs.r10_BlueprintEncode_label.ne(obs.r20_BlueprintEncode_label)
 obs['reference_patient_disagreement']=obs.r10_HPCA_label.ne(obs.patient_HPCA_label)|obs.r10_BlueprintEncode_label.ne(obs.patient_BlueprintEncode_label)
 obs['native_reference_pruning_flag']=obs.r10_HPCA_pruned_label.isin(['','NA'])|obs.r10_BlueprintEncode_pruned_label.isin(['','NA'])
 noinfo=obs.nonpositive_analysis_library.astype(str).str.lower().isin(['true','1'])
 obs.loc[noinfo,'proposed_lineage']='Unresolved';obs.loc[noinfo,'status']='unresolved_no_nonheldout_library'
 obs=obs.merge(cnv.drop(columns='sample'),on='cell_id',how='left',validate='one_to_one')
 assert obs.CNV_evidence_current.notna().all()
 parts.append(obs)
allcells=pd.concat(parts,ignore_index=True);assert len(allcells)==876648 and allcells.cell_id.is_unique
allcells['identity_interpretation']=allcells.CNV_evidence_current.map({'concordant_alteration':'RNA_CNV_alteration_supported_not_diagnostic_of_malignancy','concordant_no_alteration':'no_alteration_in_two_callers_not_proof_of_benign_identity','method_conflict':'CNV_conflict_identity_unresolved','single_method_alteration':'single_method_alteration_identity_unresolved','single_method_no_alteration':'single_method_no_alteration_identity_unresolved','unassessed':'CNV_unassessed_identity_unresolved'})
allcells.to_csv(R/'tables/all_cell_reannotation_v1.tsv.gz',sep='\t',index=False)
pd.concat(allrefs,ignore_index=True).to_csv(R/'tables/all_reference_profile_annotations.tsv.gz',sep='\t',index=False)
pd.DataFrame(evledger).to_csv(R/'resources/CNV_EVIDENCE_PROVENANCE.tsv',sep='\t',index=False)
counts=allcells.groupby(['dataset','proposed_lineage','status']).agg(cells=('cell_id','size'),patients=('patient_id','nunique'),samples=('sample','nunique')).reset_index()
counts.to_csv(R/'tables/new_annotation_counts.tsv',sep='\t',index=False)
allcells.groupby(['proposed_lineage','CNV_evidence_current']).size().rename('cells').reset_index().to_csv(R/'tables/lineage_CNV_cross_table.tsv',sep='\t',index=False)
allcells[allcells.proposed_lineage.eq('Fibroblast_like')].to_csv(R/'tables/fibroblast_like_review_candidates.tsv.gz',sep='\t',index=False)
# Old-state comparison is performed only after the complete proposed cluster annotation lock.
old=ad.read_h5ad(P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',backed='r')
oldobs=old.obs[['expanded_class_id','supported_class_id','in_supported']].copy();old.file.close()
oldobs=oldobs.reset_index(names='cell_id');cmp=oldobs.merge(allcells[['cell_id','dataset','sample','patient_id','proposed_lineage','status','CNV_evidence_current','curation_id']],on='cell_id',validate='one_to_one')
assert len(cmp)==8174
cmp.to_csv(R/'unblinding/old_frozen_to_new_annotation.tsv.gz',sep='\t',index=False)
for layer in ['expanded','supported']:
 t=cmp if layer=='expanded' else cmp[cmp.in_supported.astype(bool)]
 t.groupby([layer+'_class_id','proposed_lineage']).size().rename('cells').reset_index().to_csv(R/'unblinding'/f'old_{layer}_new_lineage_counts.tsv',sep='\t',index=False)
summary={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cells':len(allcells),'samples':allcells['sample'].nunique(),'sources':allcells.dataset.nunique(),'primary_clusters':len(cur),'all_clusters_reviewed':True,'labels_status':'proposed_v1_pending_external_blind_review','old_labels_not_used_before_curation_lock':True,'curation_lock_SHA256':hashlib.sha256((R/'resources/CURATION_LOCKED.json').read_bytes()).hexdigest(),'old_frozen_comparison_rows':len(cmp),'no_cell_removed':True,'lineage_and_CNV_identity_separate':True,'counts_by_lineage':allcells.proposed_lineage.value_counts().to_dict(),'fibroblast_like_review_cells':int(allcells.proposed_lineage.eq('Fibroblast_like').sum())}
(R/'delivery/ANNOTATION_INTEGRATION_QA.json').write_text(json.dumps(summary,indent=2));(R/'ANNOTATION_INTEGRATION_COMPLETE').write_text('SUCCESS');print(json.dumps(summary),flush=True)
