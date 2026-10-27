
from pathlib import Path
import os,sys,json,datetime,shutil,traceback,zipfile
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase69_all_cluster_matrix_20261006'
for d in ['code','tables','resources','QA','delivery','logs']:(O/d).mkdir(parents=True,exist_ok=True)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def run():
 manifest=json.loads((B/'CAF_CURRENT_CELL_EVIDENCE.json').read_text())
 ob=pd.read_csv(manifest['cell_table'],sep='\t',low_memory=False).set_index('cell_id')
 s=pd.read_csv(B/'phase62_F1F5_internal_20261005/tables/F1F5_subgroup_summary.tsv',sep='\t').set_index('subgroup_id')
 z=dict(np.load(B/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz'));assert z['cells'].tolist()==ob.index.tolist()
 genes=z['genes'];idx={g:i for i,g in enumerate(genes)}
 old=pd.read_csv(B/'phase66_doublet_repair_DSC1_20261005/tables/gene_before_after_comparison.tsv.gz',sep='\t')
 prog=pd.read_csv(B/'phase66_doublet_repair_DSC1_20261005/tables/program_before_after_comparison.tsv',sep='\t')
 cfg=json.loads((B/'phase66_doublet_repair_DSC1_20261005/resources/PROGRAMS.json').read_text())
 roles=['collagen_biosynthesis','collagen_assembly','collagen_crosslinking','collagen_components','ECM_organization']
 related=set().union(*(set(cfg['sets'][k]) for k in roles))
 q=old[(old.matrix=='corrected')&old.gene.isin(related)&(old.mean_effect_excluded>0)&(old.P_excluded<.05)].copy()
 q['primary_source']=q.subgroup_id.map(s.primary_source);q['is_primary_source']=q.dataset.eq(q.primary_source)
 display=['COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A3','COL4A1','COL4A2','ADAMTS14','PCOLCE','PLOD2','SERPINH1','LOX','LOXL2','PXDN','FN1','SPARC','BGN','DCN','LUM','SDC2','TNC','VCAN','ACAN','ASPN','ADAMTSL1','MMP2','MMP11','PPIB','ACTB','ACTG1']
 display=[g for g in display if g in idx]
 selected=sorted(set(q.gene)|set(display));jj=[idx[g] for g in selected]
 dump(O/'resources/PROTOCOL.json',{'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task':'All 33 subgroups, every observed source: collagen assembly/crosslinking and matrix-related expression','metadata':manifest,'comparison':'Within-source/donor/site target subgroup vs all other candidate cells; equal weight sites then donors. No pooling sources into independent patients.','screening':'All sources, corrected, after excluding 26 model-positive doublets; positive effect and inherited P<0.05; full original family FDR unchanged. Broad pathway membership does not establish collagen-specific function.','sensitivity':'Reconstruct raw/corrected baseline all 658 genes; add descriptive corrected minimum5 and Fibro-only minimum1 arms for selected genes; no new P/FDR.','gene_effect':'Mean over patients of log2(target pseudobulk CP10k+0.1)-log2(reference pseudobulk CP10k+0.1)','unassayed':'Missing, never zero','expression':'Within-source descriptive per-cell display; not inferential evidence or physical contribution fractions.','selected_genes':selected,'display_genes':display,'deferred':['drug_resistance','inhibition'],'scope_change':'Removes previous primary-source-only filter; retains all-source contradictions and low-replication descriptive evidence.'})
 dump(O/'resources/PROGRAMS.json',cfg)
 for n in ['CAF_NEXT_STAGE_CURRENT.json','PROJECT_STATE.json','CAF_RESUME_STATUS_CURRENT.json']:
  p=B/n;backup=O/'resources'/('PREVIOUS_'+n)
  if not backup.exists():shutil.copy2(p,backup)
  d=json.loads(p.read_text());d.update({'active_task':{'title':'All-cluster matrix evidence expansion','phase':O.name,'status':'running'},'task_status':'all_cluster_matrix_running'});dump(p,d)
 old.to_csv(O/'tables/inherited_all_gene_statistics.tsv.gz',sep='\t',index=False)
 old[old.gene.isin(related)].to_csv(O/'tables/all_source_related_gene_statistics.tsv.gz',sep='\t',index=False)
 prog.to_csv(O/'tables/inherited_program_statistics.tsv.gz',sep='\t',index=False)
 ps=prog[prog.program.isin(roles+['fibrillar_I_III_V','collagen_IV','cartilage_matrix','LOX_family'])]
 ps.to_csv(O/'tables/all_source_program_statistics.tsv',sep='\t',index=False)
 s.to_csv(O/'tables/subgroup_summary.tsv',sep='\t')
 patient=[];sens=[];expr=[];composition=[];maxerr=0;baseline_rows=0;pairs=0
 keep=~ob.doublet_positive_any.to_numpy(bool)
 def avg(x):
  n=np.isfinite(x).sum(0);return np.divide(np.nansum(x,0),n,out=np.full(x.shape[1],np.nan),where=n>0)
 for st,sr in s.iterrows():
  for ds in sorted(ob.loc[ob.subgroup_id.eq(st),'dataset'].unique()):
   pairs+=1;source=ob.dataset.eq(ds).to_numpy()
   base=source&ob.subgroup_id.eq(st).to_numpy()&keep
   for lineage,n in ob.loc[base,'lineage_stratum'].value_counts().items():
    composition.append({'subgroup_id':st,'dataset':ds,'lineage':lineage,'cells':n,'total_cells':int(base.sum())})
   for mode,scope,minimum in [('raw','all',1),('corrected','all',1),('corrected','all',5),('corrected','Fibro',1)]:
    arm=f'{mode}_excluded_{scope}_min{minimum}'
    target=base.copy();rest=source&ob.subgroup_id.ne(st).to_numpy()&keep
    if scope=='Fibro':
     fibro=ob.lineage_stratum.eq('Fibro').to_numpy();target&=fibro;rest&=fibro
    x=z[mode];lib=z[mode+'_library'];effects=[]
    for donor,dg in ob.loc[source].groupby('donor_record',sort=True):
     sd=[];nt=nr=0
     for site,sg in dg.groupby('site',dropna=False):
      ix=ob.index.get_indexer(sg.index);il=ix[target[ix]];ir=ix[rest[ix]]
      if min(len(il),len(ir))<minimum:continue
      covered=z['assayed'][np.r_[il,ir]].all(0)
      t=x[il].sum(0)/lib[il].sum()*1e4;r=x[ir].sum(0)/lib[ir].sum()*1e4
      effect=np.log2(t+.1)-np.log2(r+.1);effect[~covered]=np.nan
      sd.append(effect);nt+=len(il);nr+=len(ir)
     if sd:
      dv=avg(np.array(sd));effects.append(dv)
      for j in jj:patient.append({'subgroup_id':st,'dataset':ds,'arm':arm,'donor_record':donor,'gene':genes[j],'sites':len(sd),'target_cells':nt,'reference_cells':nr,'effect':dv[j]})
    gd=np.asarray(effects).reshape((-1,len(genes)));means=avg(gd);nn=np.isfinite(gd).sum(0);npos=(gd>0).sum(0)
    if minimum==1 and scope=='all':
     ref=old[(old.subgroup_id==st)&(old.dataset==ds)&(old.matrix==mode)].set_index('gene').reindex(genes)
     assert ref.index.is_unique and len(ref)==658
     assert np.array_equal(nn,ref.N_excluded.fillna(0).to_numpy()),(st,ds,mode,'N')
     assert np.array_equal(npos,ref.positive_patients_excluded.fillna(0).to_numpy()),(st,ds,mode,'positiveN')
     assert np.allclose(means,ref.mean_effect_excluded.to_numpy(),equal_nan=True,atol=1e-9),(st,ds,mode,'effect')
     diff=np.abs(means-ref.mean_effect_excluded.to_numpy());maxerr=max(maxerr,float(np.nanmax(diff)) if np.isfinite(diff).any() else 0);baseline_rows+=int(ref.N_excluded.notna().sum())
    for j in jj:sens.append({'subgroup_id':st,'dataset':ds,'gene':genes[j],'arm':arm,'N':int(nn[j]),'mean_effect':means[j],'positive_patients':int(npos[j]),'direction_agreement':npos[j]/nn[j] if nn[j] else np.nan,'new_P_computed':False})
   cells=np.flatnonzero(base)
   for mode in ['raw','corrected']:
    for gene in display:
     j=idx[gene];ids=cells[z['assayed'][cells,j]];v=z[mode][ids,j];lib=z[mode+'_library'][ids]
     expr.append({'subgroup_id':st,'dataset':ds,'primary_source':sr.primary_source,'is_primary_source':ds==sr.primary_source,'matrix':mode,'gene':gene,'cells':len(cells),'assayed_cells':len(ids),'detection_ge1_pct':float((v>=1).mean()*100) if len(v) else np.nan,'mean_log1p_CP10k':float(np.log1p(v/lib*1e4).mean()) if len(v) else np.nan,'mean_CP10k':float((v/lib*1e4).mean()) if len(v) else np.nan})
  print('DONE',st,flush=True)
 se=pd.DataFrame(sens);pat=pd.DataFrame(patient)
 se.to_csv(O/'tables/all_source_gene_sensitivity.tsv.gz',sep='\t',index=False);pat.to_csv(O/'tables/patient_gene_effects.tsv.gz',sep='\t',index=False)
 pd.DataFrame(expr).to_csv(O/'tables/all_source_expression.tsv',sep='\t',index=False)
 pd.DataFrame(composition).to_csv(O/'tables/all_source_lineage_composition.tsv',sep='\t',index=False)
 for arm,prefix in [('raw_excluded_all_min1','raw'),('corrected_excluded_all_min5','min5'),('corrected_excluded_Fibro_min1','fibro')]:
  t=se[se.arm.eq(arm)][['dataset','subgroup_id','gene','N','mean_effect','positive_patients']].rename(columns={k:prefix+'_'+k for k in ['N','mean_effect','positive_patients']})
  q=q.merge(t,on=['dataset','subgroup_id','gene'],validate='one_to_one',how='left')
 lo=[]
 for r in q.itertuples():
  v=pat[(pat.subgroup_id==r.subgroup_id)&(pat.dataset==r.dataset)&(pat.gene==r.gene)&(pat.arm=='corrected_excluded_all_min1')].effect.dropna().to_numpy()
  x=(v.sum()-v)/(len(v)-1) if len(v)>1 else np.array([])
  lo.append(float(x.min()) if len(x) else np.nan)
 q['LOPO_min_mean']=lo
 q.to_csv(O/'tables/nominal_upregulated_all_sources_with_sensitivity.tsv',sep='\t',index=False)
 extra=q[~q.is_primary_source];extra.to_csv(O/'tables/additional_nonprimary_source_signals.tsv',sep='\t',index=False)
 primary=q[q.is_primary_source]
 assert len(primary)==45 and len(q)==54 and len(extra)==9 and len(s)==33
 pd.DataFrame([{'subgroup_id':st,'primary_source':sr.primary_source,'working_name':sr.working_name,'all_cells':sr.cells,'all_source_nominal_gene_rows':len(q[q.subgroup_id==st]),'nominal_sources':';'.join(sorted(q.loc[q.subgroup_id==st,'dataset'].unique())),'nominal_genes':';'.join(sorted(q.loc[q.subgroup_id==st,'gene'].unique()))} for st,sr in s.iterrows()]).to_csv(O/'tables/all_33_subgroup_audit.tsv',sep='\t',index=False)
 qa={'passed':True,'subgroups':len(s),'source_subgroup_pairs':pairs,'baseline_gene_rows_reproduced':baseline_rows,'max_effect_error':maxerr,'all_source_nominal_rows':len(q),'new_nonprimary_rows':len(extra),'inherited_P_FDR_unchanged':True,'new_sensitivity_no_P':True,'metadata_cells':len(ob),'excluded_model_positive':int((~keep).sum())}
 dump(O/'QA/ACCEPTANCE.json',qa);(O/'ANALYSIS_COMPLETE').write_text('SUCCESS\n')
 with zipfile.ZipFile(O/'analysis_for_report.zip','w',zipfile.ZIP_DEFLATED) as a:
  for folder in ['tables','resources','QA']:
   for f in (O/folder).iterdir():
    if f.is_file():a.write(f,f.relative_to(O).as_posix())
 print(json.dumps(qa),flush=True)
if __name__=='__main__':
 status=1
 try:run();status=0
 except Exception:traceback.print_exc()
 finally:dump(O/'QA/EXECUTION_STATUS.json',{'exit_code':status,'at':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 sys.exit(status)
