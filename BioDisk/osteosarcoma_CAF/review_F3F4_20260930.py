from pathlib import Path
import json, hashlib, datetime
import pandas as pd
import anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
T=P/'phase30_two_layer_classification_20260926/tables'
O=P/'phase40_F3F4_matrix_review_20260930'
O.mkdir(exist_ok=True)
for sub in ['tables','code','delivery']: (O/sub).mkdir(exist_ok=True)
log=[]; summaries=[]; top=[]; focused=[]
axes={
'collagen_synthesis':['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','COL6A2','COL6A3','COL11A1','COL12A1','COL14A1','P4HA1','P4HA2','P4HB','PLOD1','PLOD2','PLOD3','SERPINH1'],
'matrix_assembly_crosslinking':['FN1','POSTN','SPARC','TNC','THBS1','THBS2','VCAN','DCN','LUM','FBN1','ELN','LOX','LOXL1','LOXL2','LOXL3','LOXL4','LAMA4','LAMB1','LAMC1'],
'matrix_remodelling':['MMP2','MMP9','MMP11','MMP14','TIMP1','TIMP2','TIMP3','ADAM12','ADAMTS2','ADAMTS4','CTHRC1','FAP'],
'signalling_candidates':['TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','SMAD2','SMAD3','SMAD7','CTGF','CCN2','CXCL12','CXCL14','CCL2','IL6','PTGS2','CD274','PDCD1LG2','LGALS1','LGALS3','NT5E','ENTPD1','IDO1','DDR1','DDR2','ITGA11','ITGB1','YAP1','WWTR1','TEAD1'],
'other_direction_examples':['VEGFA','ANGPT1','ANGPT2','HIF1A','CA9','SLC2A1','LDHA','HGF','IGF1','WNT5A','RSPO3','RUNX2','SP7','ALPL','IBSP','BGLAP','MGP']}
reverse={}
for axis,genes in axes.items():
 for gene in genes: reverse.setdefault(gene,[]).append(axis)
files=sorted(T.glob('*_all_gene_statistics.tsv.gz'))
all_rows=[]
for f in files:
 d=pd.read_csv(f,sep='\t'); g=d[d.class_id.isin(['F3','F4'])].copy()
 if g.empty: continue
 assert not g.duplicated(['layer','pool','class_id','feature']).any(), f.name
 g['source_table']=f.name
 g['direction']=g.median_delta.map(lambda v:'up' if v>0 else ('down' if v<0 else 'zero_or_NA'))
 g['raw_P_lt_0_05']=g.P.lt(.05)
 g['FDR_lt_0_05']=g.FDR.lt(.05)
 g['illustrative_axis_annotation']=g.symbol.map(lambda x:';'.join(reverse.get(x,[])))
 all_rows.append(g)
 for (layer,pool,cls),h in g.groupby(['layer','pool','class_id']):
  n=h.n_independent_donors.unique(); assert len(n)==1
  up=h[h.P.lt(.05)&h.median_delta.gt(0)]
  down=h[h.P.lt(.05)&h.median_delta.lt(0)]
  summaries.append(dict(layer=layer,pool=pool,class_id=cls,N=int(n[0]),genes=len(h),valid_P=int(h.P.notna().sum()),nominal_up=len(up),nominal_down=len(down),FDR_up=int((up.FDR<.05).sum()),FDR_down=int((down.FDR<.05).sum()),min_P=h.P.min(),positive_effect_genes=int(h.median_delta.gt(0).sum()),comparison='same patient and site: target class vs all other candidate classes'))
  r=up.sort_values(['P','median_delta'],ascending=[True,False]).head(30)
  top.append(r)
 focused.append(g[g.symbol.isin(reverse)])
full=pd.concat(all_rows,ignore_index=True)
full.to_csv(O/'tables/F3F4_all_existing_gene_statistics.tsv.gz',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(O/'tables/F3F4_analysis_availability.tsv',sep='\t',index=False)
pd.concat(top,ignore_index=True).to_csv(O/'tables/F3F4_top30_nominal_up_by_pool.tsv',sep='\t',index=False)
pd.concat(focused,ignore_index=True).to_csv(O/'tables/F3F4_mechanistic_annotation_examples.tsv.gz',sep='\t',index=False)
programs=pd.read_csv(T/'class_program_statistics.tsv',sep='\t')
programs=programs[programs.class_id.isin(['F3','F4'])]
programs.to_csv(O/'tables/F3F4_existing_program_statistics.tsv',sep='\t',index=False)
profiles=pd.read_csv(T/'class_profiles.tsv',sep='\t')
profiles[profiles.class_id.isin(['F3','F4'])].to_csv(O/'tables/F3F4_class_profiles.tsv',sep='\t',index=False)
coverage=[]
for layer in ['supported','expanded']:
 a=ad.read_h5ad(T.parent/'objects'/f'{layer}_integrated.h5ad',backed='r')
 coverage.append(dict(layer=layer,cells=a.n_obs,genes=a.n_vars,symbol_column='symbol' in a.var.columns,absent_examples=[x for x in reverse if x not in set(a.var['symbol'].astype(str))],var_columns=list(a.var.columns)))
 a.file.close()
qa={'at':datetime.datetime.now().astimezone().isoformat(),'scope':'existing full-gene statistics review; no new tests or enrichment executed','input_tables':len(files),'F3F4_rows':len(full),'availability_rows':len(summaries),'model_coverage':coverage,'P_FDR_families_reused':True,'new_gene_filter':False,'mechanistic_examples_are_annotations_not_gene_universe':True,'layers_not_independent_cohorts':True,'cohort_pools_overlap_do_not_combine_N':True,'F4_identity_conflict_retained':True,'all_genes_claim_boundary':'existing integrated gene universe; next stage must audit per-source full gene coverage','protocol':'within patient and site class versus rest; not F3 versus F4, not metastasis contrast','immune_exclusion_or_shell_geometry_tested':False,'required_next':'per-source full genes; unbiased pathway analysis; direct F3/F4 contrasts where evaluable; spatial tumor-edge and CD8 neighborhood validation'}
(O/'delivery/REVIEW_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print('AVAILABILITY_MAIN',pd.DataFrame(summaries).query("pool == 'SCPCP000017'").to_string(index=False))
print('COVERAGE',json.dumps(coverage,ensure_ascii=False))
cols=['layer','class_id','symbol','n_independent_donors','median_delta','positive_donors','P','FDR']
for layer,cls in [('supported','F3'),('expanded','F3'),('expanded','F4')]:
 h=full[(full.layer==layer)&(full.pool=='SCPCP000017')&(full.class_id==cls)]
 print('TOP',layer,cls,h[h.P.lt(.05)&h.median_delta.gt(0)].sort_values(['P','median_delta'],ascending=[True,False])[cols].head(15).to_string(index=False))
 print('MATRIX_EXAMPLES',h[h.symbol.isin(reverse)&h.P.lt(.05)&h.median_delta.gt(0)][cols].sort_values('P').head(20).to_string(index=False))
print('PROGRAMS_MAIN',programs[programs.pool.eq('SCPCP000017')].to_string(index=False))
manifest=[{'file':str(f.relative_to(O)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted((O/'tables').iterdir())]
(O/'delivery/EVIDENCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
(O/'COMPLETE').write_text('SUCCESS\n')
print('F3F4_EXISTING_GENOMEWIDE_REVIEW_COMPLETE')