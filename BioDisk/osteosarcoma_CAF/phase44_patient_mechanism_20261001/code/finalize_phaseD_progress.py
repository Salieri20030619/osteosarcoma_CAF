from pathlib import Path
import json,datetime,hashlib,pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');C=P/'phase43_full_transcriptome_discovery_20261001';D=P/'phase44_patient_mechanism_20261001';L=P/'phase45_LR_expression_20261001'
assert (C/'controller.exit_code').read_text().strip()=='0' and (C/'PIPELINE_COMPLETE').read_text().strip()=='SUCCESS' and (C/'C2_COMPLETE').read_text().strip()=='SUCCESS'
assert (D/'qa_resume.exit_code').read_text().strip()=='0' and (D/'CORE_COMPLETE').read_text().strip()=='SUCCESS'
assert (L/'qa_resume.exit_code').read_text().strip()=='0' and (L/'PIPELINE_COMPLETE').read_text().strip()=='SUCCESS'
assert (D/'clinical_site_review/INDEPENDENT_R_QA_COMPLETE').read_text().strip()=='SUCCESS' and (D/'clinical_site_review/model_sensitivity/sensitivity.exit_code').read_text().strip()=='0'
rows=[]
for layer,cls in [('supported','F3'),('expanded','F3'),('expanded','F4')]:
 ctx=f'{layer}__SCPCP000017__{cls}_vs_rest'
 for lib in ['hallmark','reactome','go_bp','matrisome']:
  f=C/'tables/pathways'/f'{ctx}__{lib}__GSEA.tsv.gz'
  if f.exists():
   x=pd.read_csv(f,sep='\t');x['NOM p-val']=pd.to_numeric(x['NOM p-val']);x['NES']=pd.to_numeric(x.NES);x['FDR q-val']=pd.to_numeric(x['FDR q-val']);x.to_csv(C/'delivery'/f'{ctx}__{lib}__full_GSEA.tsv.gz',sep='\t',index=False)
   select=x[x['NOM p-val'].lt(.05)].sort_values(['NOM p-val','NES'],ascending=[True,False]);select.to_csv(C/'delivery'/f'{ctx}__{lib}__nominal_GSEA.tsv',sep='\t',index=False)
   rows.append({'layer':layer,'class_id':cls,'library':lib,'terms':len(x),'model_nominal':int(x['NOM p-val'].lt(.05).sum()),'native_model_FDR':int(x['FDR q-val'].lt(.05).sum()),'gene_set_model_P_not_patient_P':True})
pd.DataFrame(rows).to_csv(C/'delivery/main_F3F4_pathway_result_summary.tsv',sep='\t',index=False)
# Produce a quantitative tie-order diagnostic for the primary class programs, without inventing biological P values.
import numpy as np
cross=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='\t',keep_default_na=False);quality=cross.set_index('feature').mapping_method.eq('HGNC_Ensembl_approved').map({True:0,False:1}).to_dict();diagnostic=[]
for layer,cls in [('supported','F3'),('expanded','F3'),('expanded','F4')]:
 frame=pd.read_csv(C/'tables'/f'{layer}__SCPCP000017__full_gene_statistics.tsv.gz',sep='\t');g=frame[frame.contrast.eq(cls+'_vs_rest')&frame.P.notna()&frame.N.ge(2)&frame.symbol.notna()&frame.symbol.ne('')].copy();g['quality']=g.feature.map(quality);g=g.sort_values(['symbol','quality','feature']).drop_duplicates('symbol');frac=1-g.mean_delta.nunique()/len(g)
 for lib in ['hallmark','reactome','matrisome']:
  native=pd.read_csv(C/'tables/pathways'/f'{layer}__SCPCP000017__{cls}_vs_rest__{lib}__GSEA.tsv.gz',sep='\t')
  # Native output leading-edge genes are retained; this diagnostic only checks observed ES against opposite lexical ordering inside equal scores.
  sets={}
  fname='curated_cgp.gmt' if lib=='matrisome' else lib+'.gmt'
  for line in (P/'phase41_full_transcriptome_setup_20261001/resources'/fname).read_text().splitlines():
   parts=line.split('\t')
   if lib=='matrisome' and not parts[0].startswith('NABA_'):continue
   sets[parts[0]]=set(parts[2:])
  order1=g.sort_values(['mean_delta','symbol'],ascending=[False,True]);order2=g.sort_values(['mean_delta','symbol'],ascending=[False,False])
  for row in native.itertuples(index=False):
   term=row.Term
   if lib=='reactome' and not any(s in term for s in ['COLLAGEN','EXTRACELLULAR','TGF','IMMUNE','INTERLEUKIN','MATRIX']):continue
   members=sets.get(term,set());scores=[]
   for ordered in [order1,order2]:
    hit=ordered.symbol.isin(members).to_numpy();w=np.abs(ordered.mean_delta.to_numpy());den=(w*hit).sum();miss=(~hit).sum()
    if not hit.any() or den<=0 or miss<=0:scores.append(np.nan);continue
    run=np.cumsum(w*hit/den-(~hit)/miss);scores.append(float(run[np.argmax(np.abs(run))]))
   diagnostic.append({'layer':layer,'class_id':cls,'library':lib,'term':term,'background_tie_fraction':frac,'ES_lexical_ascending':scores[0],'ES_lexical_descending':scores[1],'ES_absolute_change':abs(scores[0]-scores[1]),'ES_sign_changed':bool(np.sign(scores[0])!=np.sign(scores[1])) if np.isfinite(scores).all() else False,'scope':'observed_ES_order_diagnostic_not_new_P; current_HGNC_aliases_in_primary_analysis_remain_authoritative'})
pd.DataFrame(diagnostic).to_csv(C/'delivery/main_pathway_tie_order_ES_diagnostic.tsv',sep='\t',index=False)
status={'at':datetime.datetime.now().astimezone().isoformat(),'current_phase':'phase44/45 core mechanism and clinical work accepted; C2 pathway computation completed','C2':json.loads((C/'delivery/OPEN_PATHWAY_QA.json').read_text()),'D_core':json.loads((D/'delivery/ROBUSTNESS_TF_QA.json').read_text()),'clinical_preferred':str(D/'clinical_site_review'),'clinical_QA':json.loads((D/'clinical_site_review/delivery/CLINICAL_FRACTION_QA.json').read_text()),'clinical_reporting_extension':json.loads((D/'clinical_reporting_extension/delivery/CLINICAL_FRACTION_QA.json').read_text()),'clinical_sensitivity':json.loads((D/'clinical_site_review/model_sensitivity/delivery/SENSITIVITY_QA.json').read_text()),'LR':json.loads((L/'delivery/LR_EXPRESSION_QA.json').read_text()),'exit_codes':{'C2':0,'D_core_QA_resume':0,'LR_QA_resume':0,'clinical_site_review':0,'clinical_sensitivity':0},'proxy_used':False,'key_findings':['F3 COL1A1/LOX/MMP14 mean direction stays positive in every one-patient deletion; N9 nominal, not FDR significant','main F4 CBFB/TGIF1 signed target response FDR0.042969, 3/11 targets; inferred regulon response, not direct activity','strict-site unadjusted fraction tests no raw P<0.05; primary adjusted F3/F4 signals fail t-reference sensitivity; all original model P/FDR retained','RNA LR compatibility nominal candidates retained; no communication or immune inhibition experimentally established'],'remaining':['final functional synthesis after review of all pathway output including tie diagnostic and small-N model statistics','independent source/heldout phenotype migration not claimed completed','spatial tumor boundary/image acceptance and immune-exclusion testing','cross-source global patient deduplication','final client report'],'preferred_status_files':{'D':str(D/'delivery/CORE_ACCEPTED_STATUS.json'),'LR':str(L/'delivery/QA_ACCEPTED_STATUS.json'),'C2':str(C/'PIPELINE_STATUS.json')},'historical_failures':'initial CORE_PIPELINE_STATUS and LR PIPELINE_STATUS failures retained; successful QA resume records supersede them'}
(P/'CAF_NEXT_STAGE_CURRENT.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));(D/'delivery/CURRENT_PROJECT_PROGRESS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
files=[D/'delivery/main_selected_gene_patient_stability.tsv',D/'delivery/main_all_nominal_TF_response.tsv',D/'delivery/TF_regulon_coverage_review.tsv',D/'clinical_site_review/tables/all_clinical_fraction_statistics.tsv',D/'clinical_site_review/model_sensitivity/tables/all_adjusted_HC3_t_reference_sensitivity.tsv',D/'figures/selected_gene_patient_leaveout.png',D/'figures/selected_gene_patient_leaveout.pdf',C/'delivery/main_F3F4_pathway_result_summary.tsv',C/'delivery/main_pathway_tie_order_ES_diagnostic.tsv',D/'delivery/CURRENT_PROJECT_PROGRESS.json']
(D/'delivery/PROGRESS_ARTIFACT_MANIFEST.json').write_text(json.dumps([{'path':str(f),'file':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files],indent=2))
print('ALL_CURRENT_PHASE_EXIT_AND_MARKER_CHECKS_PASSED');print('PATHWAYS',pd.DataFrame(rows).to_string(index=False));print('TIE_DIAGNOSTIC',len(diagnostic),'rows',sum(x['ES_sign_changed'] for x in diagnostic),'sign changes; exploratory model P retained')