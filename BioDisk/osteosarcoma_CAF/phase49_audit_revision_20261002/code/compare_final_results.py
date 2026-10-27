from pathlib import Path
import pandas as pd,numpy as np,json
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');C=P/'phase43_full_transcriptome_discovery_20261001';R=P/'phase49_audit_revision_20261002'
assert (R/'COMPUTATION_VALIDATED').read_text().strip()=='SUCCESS'
ex=pd.read_csv(R/'tables/pathway_execution_summary.tsv',sep='	');rows=[];selected=[]
for r in ex.itertuples():
 fn=f'{r.context_id}__{r.library}__GSEA.tsv.gz';newfile=R/'tables/pathways'/fn;oldfile=C/'tables/pathways'/fn
 if not newfile.exists():continue
 new=pd.read_csv(newfile,sep='	');old=pd.read_csv(oldfile,sep='	');assert not new.Term.duplicated().any();d=new.merge(old,on='Term',suffixes=('_new','_old'),how='outer',indicator=True);both=d[d._merge.eq('both')]
 rows.append(dict(context_id=r.context_id,library=r.library,state=r.state,old_terms=len(old),new_terms=len(new),common_terms=len(both),newly_unevaluable=int(d._merge.eq('right_only').sum()),newly_evaluable=int(d._merge.eq('left_only').sum()),model_P_nominal_changes=int((both.GSEA_model_P_new.lt(.05)!=both.GSEA_model_P_old.lt(.05)).sum()),model_BH_FDR_changes=int((both.model_BH_FDR_new.lt(.05)!=both.model_BH_FDR_old.lt(.05)).sum()),max_NES_change=float(np.max(np.abs(pd.to_numeric(both.NES_new)-pd.to_numeric(both.NES_old))))))
 if r.context_id in ['supported__SCPCP000017__F3_vs_rest','expanded__SCPCP000017__F4_vs_rest']:
  for term in ['NABA_CORE_MATRISOME','NABA_COLLAGENS','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_GLYCOLYSIS','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_G2M_CHECKPOINT','HALLMARK_PROTEIN_SECRETION']:
   t=d[d.Term.eq(term)]
   if len(t):
    v=t.iloc[0];selected.append(dict(context_id=r.context_id,library=r.library,term=term,NES_old=v.NES_old,NES_new=v.NES_new,P_old=v.GSEA_model_P_old,P_new=v.GSEA_model_P_new,BH_old=v.model_BH_FDR_old,BH_new=v.model_BH_FDR_new))
s=pd.DataFrame(rows);s.to_csv(R/'tables/GSEA_before_after_summary.tsv',sep='	',index=False);pd.DataFrame(selected).to_csv(R/'tables/main_GSEA_before_after.tsv',sep='	',index=False)
ora=pd.read_csv(R/'tables/ORA_before_after_summary.tsv',sep='	');surface=pd.read_csv(R/'surface/tables/surface_old7_before_after.tsv',sep='	');qa=json.loads((R/'delivery/PATHWAY_INDEPENDENT_QA.json').read_text());nes=next(t['NES_new'] for t in selected if t['context_id']=='supported__SCPCP000017__F3_vs_rest' and t['term']=='NABA_CORE_MATRISOME')
out={'computations_accepted':True,'ORA_nominal_status_changes':int(ora.nominal_status_changed.sum()),'ORA_FDR_status_changes':int(ora.FDR_status_changed.sum()),'ORA_rows':int(ora.rows.sum()),'ORA_library_contexts':len(ora),'GSEA_recomputed':int(ex.state.eq('recomputed').sum()),'GSEA_reused':int(ex.state.eq('reused_identical_inputs').sum()),'GSEA_model_P_nominal_changes_common_terms':int(s.model_P_nominal_changes.sum()),'GSEA_model_BH_changes_common_terms':int(s.model_BH_FDR_changes.sum()),'GSEA_newly_unevaluable':int(s.newly_unevaluable.sum()),'F3_core_matrisome_NES':float(nes),'old7_FDR_status_changes':int((surface.FDR_new.lt(.05)!=surface.FDR_old.lt(.05)).sum()),'surface_candidates':10,'surface_pairs':45,'pathway_independent_QA':qa,'main_GSEA_before_after':selected};(R/'delivery/FINAL_REVISION_FACTS.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
