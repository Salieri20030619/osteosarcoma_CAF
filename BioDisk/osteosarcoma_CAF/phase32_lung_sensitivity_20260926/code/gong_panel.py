from pathlib import Path
import json,datetime
import pandas as pd,numpy as np,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase32_lung_sensitivity_20260926';M=P/'phase30_two_layer_classification_20260926'
genes=['PTGS2','IL6','CCL2','CXCL1','PTX3','TNFAIP6']
protocol={'at':datetime.datetime.now().astimezone().isoformat(),'name':'Gong2022 human PTGS2 coexpression panel','genes':genes,'source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC9830653/','doi':'10.1016/j.immuni.2022.07.001','location':'Results, Figure3H human normal-tissue coexpression paragraph','same_as_previous_frozen_human_panel':True,'not_claimed_complete_mouse_scRNA_signature':True,'mouse_top10_record':['Ptgs2','Has1','Cxcl1','Il6','Tnfaip6','Mt1','Zfp36','Ccl7','Ptx3','Ccl2'],'mouse_Mt1_not_arbitrarily_replaced_by_one_human_MT1_paralog':True,'score':'continuous mean log1p(CP10k) across exact six human genes, original full-library normalization; DecontX uses corrected full-library denominator','positivity_or_cell_entry_cutoff':None}
(O/'resources/Gong2022_human_panel_verified.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2))
a=ad.read_h5ad(M/'objects/expanded_integrated.h5ad');c=ad.read_h5ad(O/'objects/expanded_decontx_counts.h5ad')
assert a.obs_names.equals(c.obs_names) and a.var_names.equals(c.var_names)
ev=pd.read_csv(M/'tables/expanded_classified_cells.tsv.gz',sep='\t').set_index('cell_id').loc[a.obs_names]
su=pd.read_csv(M/'tables/supported_classified_cells.tsv.gz',sep='\t').set_index('cell_id')
ix=[np.flatnonzero(a.var.symbol.eq(g).to_numpy()) for g in genes];assert all(len(z)==1 for z in ix)
raw=np.column_stack([np.asarray(a.X[:,z].sum(1)).ravel() for z in ix])
den=c.obs.corrected_full_library_counts.to_numpy(float)
corrected=np.column_stack([np.log1p(np.asarray(c.X[:,z].sum(1)).ravel()/np.maximum(den,1)*1e4) for z in ix])
out=ev[['sample','dataset','donor_record','site','class_id','in_supported']].rename(columns={'class_id':'expanded_class_id'}).copy()
out['supported_class_id']=su.class_id.reindex(out.index)
for j,g in enumerate(genes):out['raw_'+g]=raw[:,j];out['corrected_'+g]=corrected[:,j]
out['raw_Gong_human6_mean']=raw.mean(1);out['corrected_Gong_human6_mean']=corrected.mean(1)
out.to_csv(O/'tables/Gong_human6_cell_scores.tsv.gz',sep='\t',index_label='cell_id')
for layer,lab in [('supported',su),('expanded',ev)]:
 v=out.loc[lab.index].copy();v['class_id']=lab.class_id
 v.groupby(['dataset','donor_record','site','class_id'])[['raw_Gong_human6_mean','corrected_Gong_human6_mean']].mean().to_csv(O/'tables'/f'{layer}_Gong_human6_patient_scores.tsv',sep='\t')
qa={'cells':len(out),'supported_cells':len(su),'all_six_genes_in_shared_measured_space':True,'raw_and_corrected_scores_finite':bool(np.isfinite(raw).all() and np.isfinite(corrected).all()),'no_signature_positivity_threshold_or_cell_filter':True,'score_not_fraction_of_COX2_positive_cells':True}
(O/'delivery/GONG_PANEL_QA.json').write_text(json.dumps(qa,indent=2));print('GONG_PANEL_COMPLETE',qa)
