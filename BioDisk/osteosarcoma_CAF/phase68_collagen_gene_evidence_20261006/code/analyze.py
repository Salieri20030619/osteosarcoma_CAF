from pathlib import Path
import os,sys,json,datetime,hashlib,shutil,traceback,zipfile
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase68_collagen_gene_evidence_20261006'
for p in ['code','tables','resources','QA','delivery','logs']:(O/p).mkdir(parents=True,exist_ok=True)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def run():
 manifest=json.loads((B/'CAF_CURRENT_CELL_EVIDENCE.json').read_text())
 ob=pd.read_csv(manifest['cell_table'],sep='\t',low_memory=False).set_index('cell_id')
 s=pd.read_csv(B/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv',sep='\t').set_index('subgroup_id')
 z=np.load(B/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz');assert z['cells'].tolist()==ob.index.tolist()
 genes=z['genes'];idx={g:i for i,g in enumerate(genes)}
 old=pd.read_csv(B/'phase66_doublet_repair_DSC1_20261005/tables/gene_before_after_comparison.tsv.gz',sep='\t')
 q=pd.read_csv(B/'priority_collagen_matrix_20261005/collagen_matrix_upregulated_primary.tsv',sep='\t')
 watch=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL5A3','LOX','LOXL1','LOXL2','LOXL3','LOXL4','PCOLCE','BMP1','PLOD2','P4HA1','SERPINH1','FN1','DCN','LUM','VCAN','COL4A1','COL4A2','PXDN']
 display=list(dict.fromkeys(['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL5A3','PCOLCE','LOX','LOXL2','FN1','VCAN','COL4A1']))
 selected=sorted(set(q.gene)|set(watch));jj=[idx[g] for g in selected]
 cfg=json.loads((B/'phase66_doublet_repair_DSC1_20261005/resources/PROGRAMS.json').read_text())
 protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'request':'Proceed with collagen fibril assembly/crosslinking and matrix-deposition-related gene upregulation; defer resistance and inhibition mechanisms.','frozen_screening':'Phase66 corrected expression, source defined by largest original subgroup cell count, mean_effect_excluded>0 and raw P_excluded<.05; retain original full-family FDR and all baseline comparisons; 45 gene-state rows,35genes.','primary_stats':'Reuse accepted Phase66 statistics without new fitting or FDR rescoping. Comparison is target subgroup versus other candidate cells within source/patient/site, not treatment versus control.','new_work':'Reconstruct patient and site expression contrasts; produce patient/UMAP-linked subgroup expression displays; add descriptive gene-level sensitivity for min5 and Fibro-only comparators, no new P/FDR for these descriptive arms.','sensitivity_arms':['raw_excluded_all_min1','corrected_excluded_all_min1','corrected_excluded_all_min5','corrected_excluded_Fibro_min1'],'screening_not_gated_by_sensitivity':'Do not silently drop candidates for unavailable min5/Fibro comparisons; report effectiveN and direction separately.','coverage':'Use assayed mask, missing not zero; effective independent donor count per source; source records not pooled into independent external validation.','selected_genes':selected,'expression_display_genes':display,'deferred':['drug_resistance','inhibition_related_mechanisms'],'gene_membership':'Inherited curated PROGRAMS.json; membership can include both structural and degradation/regulatory components; no claim every member drives deposition.','metadata_source':manifest}
 dump(O/'resources/PROTOCOL.json',protocol);dump(O/'resources/PROGRAMS.json',cfg)
 for n in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
  p=B/n;backup=O/'resources'/('PREVIOUS_'+n)
  if not backup.exists():shutil.copy2(p,backup)
  d=json.loads(p.read_text());d.update({'active_task':{'title':'Collagen and matrix gene-level evidence, patient effects and figures','phase':O.name,'status':'running'},'task_status':'collagen_gene_evidence_running'});dump(p,d)
 q.to_csv(O/'tables/frozen_nominal_upregulated_45_rows.tsv',sep='\t',index=False)
 old.to_csv(O/'tables/inherited_all_658_gene_statistics.tsv.gz',sep='\t',index=False)
 programs=pd.read_csv(B/'phase66_doublet_repair_DSC1_20261005/tables/program_before_after_comparison.tsv',sep='\t')
 programs.to_csv(O/'tables/inherited_program_statistics.tsv',sep='\t',index=False)
 s.to_csv(O/'tables/subgroup_summary.tsv',sep='\t')
 patient=[];sites=[];sens=[];expr=[];maxerr=0;baseline_rows=0
 keep=~ob.doublet_positive_any.to_numpy(bool)
 def avg(x):
  n=np.isfinite(x).sum(0);return np.divide(np.nansum(x,0),n,out=np.full(x.shape[1],np.nan),where=n>0)
 for st,sr in s.iterrows():
  ds=sr.primary_source;source=ob.dataset.eq(ds).to_numpy()
  for mode,scope,minimum in [('raw','all',1),('corrected','all',1),('corrected','all',5),('corrected','Fibro',1)]:
   arm=f'{mode}_excluded_{scope}_min{minimum}'
   target=source&ob.subgroup_id.eq(st).to_numpy()&keep;rest=source&ob.subgroup_id.ne(st).to_numpy()&keep
   if scope=='Fibro':
    fibro=ob.lineage_stratum.eq('Fibro').to_numpy();target&=fibro;rest&=fibro
   x=z[mode];lib=z[mode+'_library'];donoreffects=[];names=[]
   for donor,dg in ob.loc[source].groupby('donor_record',sort=True):
    sd=[];tc=[];rc=[];nt=nr=0
    for site,sg in dg.groupby('site',dropna=False):
     ix=ob.index.get_indexer(sg.index);il=ix[target[ix]];ir=ix[rest[ix]]
     if min(len(il),len(ir))<minimum:continue
     covered=z['assayed'][np.r_[il,ir]].all(0)
     t=x[il].sum(0)/lib[il].sum()*1e4;r=x[ir].sum(0)/lib[ir].sum()*1e4
     effect=np.log2(t+.1)-np.log2(r+.1);effect[~covered]=np.nan;t[~covered]=np.nan;r[~covered]=np.nan
     sd.append(effect);tc.append(t);rc.append(r);nt+=len(il);nr+=len(ir)
     for j in jj:sites.append({'subgroup_id':st,'dataset':ds,'arm':arm,'donor_record':donor,'site':site,'gene':genes[j],'target_cells':len(il),'reference_cells':len(ir),'target_CP10k':t[j],'reference_CP10k':r[j],'effect':effect[j]})
    if sd:
     dv=avg(np.array(sd));donoreffects.append(dv);names.append(donor);tv=avg(np.array(tc));rv=avg(np.array(rc))
     for j in jj:patient.append({'subgroup_id':st,'dataset':ds,'arm':arm,'donor_record':donor,'gene':genes[j],'sites':len(sd),'target_cells':nt,'reference_cells':nr,'target_CP10k_site_mean':tv[j],'reference_CP10k_site_mean':rv[j],'effect':dv[j]})
   gd=np.asarray(donoreffects).reshape((-1,len(genes)))
   means=avg(gd);nn=np.isfinite(gd).sum(0);npos=(gd>0).sum(0)
   if minimum==1 and scope=='all':
    ref=old[(old.subgroup_id==st)&(old.dataset==ds)&(old.matrix==mode)].set_index('gene').reindex(genes)
    assert ref.index.is_unique and len(ref)==658
    assert np.array_equal(nn,ref.N_excluded.to_numpy())
    assert np.array_equal(npos,ref.positive_patients_excluded.to_numpy())
    assert np.allclose(means,ref.mean_effect_excluded.to_numpy(),equal_nan=True,atol=1e-9)
    diff=np.abs(means-ref.mean_effect_excluded.to_numpy());maxerr=max(maxerr,float(np.nanmax(diff)) if np.isfinite(diff).any() else 0);baseline_rows+=len(ref)
   for j in jj:sens.append({'subgroup_id':st,'dataset':ds,'gene':genes[j],'arm':arm,'N':int(nn[j]),'mean_effect':means[j],'positive_patients':int(npos[j]),'direction_agreement':npos[j]/nn[j] if nn[j] else np.nan,'new_P_computed':False})
  # Descriptive per-cell expression within the existing primary source; no cell-level inference.
  cells=np.flatnonzero(source&ob.subgroup_id.eq(st).to_numpy()&keep)
  for mode in ['raw','corrected']:
   for gene in display:
    j=idx[gene];ids=cells[z['assayed'][cells,j]];v=z[mode][ids,j];lib=z[mode+'_library'][ids]
    expr.append({'subgroup_id':st,'primary_source':ds,'matrix':mode,'gene':gene,'cells':len(cells),'assayed_cells':len(ids),'detection_ge1_pct':float((v>=1).mean()*100) if len(v) else np.nan,'detection_gt0_pct':float((v>0).mean()*100) if len(v) else np.nan,'mean_log1p_CP10k':float(np.log1p(v/lib*1e4).mean()) if len(v) else np.nan,'mean_CP10k':float((v/lib*1e4).mean()) if len(v) else np.nan})
  print('SUBGROUP_COMPLETE',st,ds,flush=True)
 pat=pd.DataFrame(patient);se=pd.DataFrame(sens)
 pat.to_csv(O/'tables/patient_gene_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(sites).to_csv(O/'tables/patient_site_gene_effects.tsv.gz',sep='\t',index=False)
 se.to_csv(O/'tables/gene_direction_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(expr).to_csv(O/'tables/subgroup_expression_display.tsv',sep='\t',index=False)
 rank=q.copy()
 for arm,prefix in [('raw_excluded_all_min1','raw'),('corrected_excluded_all_min5','min5'),('corrected_excluded_Fibro_min1','fibro')]:
  t=se[se.arm.eq(arm)][['dataset','subgroup_id','gene','N','mean_effect','positive_patients']].rename(columns={k:prefix+'_'+k for k in ['N','mean_effect','positive_patients']})
  rank=rank.merge(t,on=['dataset','subgroup_id','gene'],validate='one_to_one',how='left')
 rank['within_patient_direction_fraction']=rank.positive_patients_excluded/rank.N_excluded
 # LOPO direction is descriptive robustness, not independent validation.
 lopo=[]
 for row in rank.itertuples(index=False):
  vals=pat[(pat.subgroup_id==row.subgroup_id)&(pat.dataset==row.dataset)&(pat.gene==row.gene)&(pat.arm=='corrected_excluded_all_min1')].effect.dropna().to_numpy()
  leave=(vals.sum()-vals)/(len(vals)-1) if len(vals)>1 else np.array([])
  lopo.append({'dataset':row.dataset,'subgroup_id':row.subgroup_id,'gene':row.gene,'LOPO_min_mean':float(leave.min()) if len(leave) else np.nan,'LOPO_all_means_positive':bool((leave>0).all()) if len(leave) else None})
 rank=rank.merge(pd.DataFrame(lopo),on=['dataset','subgroup_id','gene'],validate='one_to_one')
 rank.to_csv(O/'tables/candidate_gene_evidence_with_sensitivity.tsv',sep='\t',index=False)
 checks={'passed':True,'baseline_gene_rows_reproduced':baseline_rows,'max_effect_error':maxerr,'frozen_upregulated_rows':len(q),'unique_frozen_genes':q.gene.nunique(),'subgroups':len(s),'metadata_cells':len(ob),'model_positive_excluded':int((~keep).sum()),'all_45_candidates_retained':len(rank)==45,'inherited_P_FDR_unchanged':True,'new_sensitivity_is_descriptive':True}
 assert len(rank)==45 and q.gene.nunique()==35
 dump(O/'QA/ACCEPTANCE.json',checks)
 summary={'candidate_counts':q.groupby('subgroup_id').size().to_dict(),'nominal_genes':q.gene.nunique(),'all45_LOPO_positive':int(rank.LOPO_all_means_positive.sum()),'fibro_negative_or_zero':rank[(rank.fibro_N>0)&(rank.fibro_mean_effect<=0)][['subgroup_id','gene','fibro_N','fibro_mean_effect']].to_dict('records'),'min5_unavailable':int(rank.min5_N.eq(0).sum()),'F2_C3_key_genes':rank[(rank.subgroup_id=='F2_C3')&rank.gene.isin(['COL3A1','COL5A1','PCOLCE','LOXL2','FN1','VCAN'])][['gene','N_excluded','mean_effect_excluded','positive_patients_excluded','P_excluded','FDR_excluded','raw_mean_effect','min5_N','min5_mean_effect','fibro_N','fibro_mean_effect','LOPO_min_mean']].to_dict('records')}
 dump(O/'resources/SUMMARY.json',summary)
 (O/'ANALYSIS_COMPLETE').write_text('SUCCESS\n')
 print(json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=='__main__':
 status=1
 try:run();status=0
 except Exception:traceback.print_exc()
 finally:dump(O/'QA/EXECUTION_STATUS.json',{'exit_code':status,'at':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 sys.exit(status)
