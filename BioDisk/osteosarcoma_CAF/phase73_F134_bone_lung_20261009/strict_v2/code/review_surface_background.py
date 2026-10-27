from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,datetime,importlib.util,traceback
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase73_F134_bone_lung_20261009/strict_v2'
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def save(d,name):d.to_csv(O/'tables'/name,sep='\t',index=False,na_rep='NA')
def dump(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def main():
 assert (O/'core.exit_code').read_text().strip()=='0'
 spec=importlib.util.spec_from_file_location('analysis',P/'phase73_F134_bone_lung_20261009/code/analyze_bone_lung_v2.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 candidates=read(O/'tables/nominal_surface_marker_candidates.tsv')
 assert len(candidates)>0
 pops=candidates[['layer','dataset','target','level','parent']].drop_duplicates().to_dict('records')
 cm=read(O/'resources/clinical_sample_manifest.tsv',keep_default_na=False).set_index('sample')
 annotation=read(P/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz',usecols=['original_cell_id','sample','review_lineage_code']).set_index('original_cell_id')
 frozen=read(P/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',usecols=['cell_id']).cell_id
 layers={layer:read(O/'tables'/f'{layer}_fixed_cells_with_sites.tsv.gz').set_index('cell_id') for layer in ['supported','expanded']}
 groups={'Fibroblast_like_other':['Fibro'],'Osteogenic_chondrogenic':['Osteo','Chondro','Osteochondro'],'Mural':['Mural','SMC'],'Endothelial':['Endo'],'T_NK':['T','T_NK','NK','Treg'],'B_Plasma':['B','Plasma'],'Myeloid':['Macro','Myeloid','MonoDC','Mono','DC','pDC'],'Epithelial':['Epi','Renal'],'Mast':['Mast']}
 mapping={v:k for k,vs in groups.items() for v in vs};bg_names=list(groups)+['Other_unresolved','Other_frozen_CAF']
 protocol={'frozen_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'candidate_source':'All nominal surface candidates from strict_v2, no top-N filter','backgrounds':groups,'remaining_annotations':'Other_unresolved','background_scope':'All source-called cells outside original8174 grouped by frozen Phase55 broad lineage; other cells within same frozen layer are a separate CAF competitor','target_scope':'Exact source/layer/class or fixed internal state that generated site association','matching':'Same source, patient and anatomical site, restricted to samples containing target cells; repeated specimens pooled; both-site summary averages site effects within patient','normalization':'Raw full-library pseudobulk log2(CP10k+0.1), target minus competitor','test':'Paired patient signflip; N<2/all-zero effect no P; same target-source-site candidate x background family BH','RNA_detection':'Raw count>0, descriptive; not measured flow positivity or purity','identity':'Background lineage and targets remain expression/CNV candidates, not histological ground truth'}
 dump(O/'resources/SURFACE_BACKGROUND_PROTOCOL.json',protocol)
 source_records=[];data={};checks=[];target_detection=[]
 for ds in sorted(candidates.dataset.unique()):
  ss=candidates[candidates.dataset.eq(ds)][['feature','symbol']].drop_duplicates('feature').sort_values('feature');features=ss.feature.to_numpy(str);symbols=ss.symbol.to_numpy(str)
  ds_pops=[r for r in pops if r['dataset']==ds]
  for sm,c in cm[cm.dataset.eq(ds)&cm.strict_site_group.isin(['Bone','Lung'])].iterrows():
   # Samples without any requested target cannot contribute matched background evidence.
   target_ids={}
   for pop in ds_pops:
    ob=layers[pop['layer']];z=ob[ob['sample'].eq(sm)]
    q=z[z.analysis_class.eq(pop['target'])] if pop['level']=='main' else z[z.analysis_internal.eq(pop['target'])]
    if len(q):target_ids[(pop['layer'],pop['target'])]=q.index
   if not target_ids:continue
   src=ad.read_h5ad(P/'phase26_client_reanalysis_20260922/inputs'/sm/'full_called_counts.h5ad',backed='r')
   obs_ids=src.obs_names;vi=src.var_names.get_indexer(features);valid=vi>=0
   panel=np.zeros((src.n_obs,len(features)),dtype=float)
   if valid.any():panel[:,valid]=src.X[:,vi[valid]].toarray()
   lib=read(P/'phase26_client_reanalysis_20260922/inputs'/sm/'cells.tsv',usecols=['cell_id','n_counts']).set_index('cell_id').loc[obs_ids,'n_counts'].to_numpy(float)
   assert len(lib)==src.n_obs and np.all(lib>0)
   lineage=annotation.loc[obs_ids,'review_lineage_code'].map(mapping).fillna('Other_unresolved').to_numpy()
   outside=~obs_ids.isin(frozen)
   def aggregate(ix):
    ix=np.asarray(ix,dtype=int);return dict(counts=panel[ix].sum(0),positive=(panel[ix]>0).sum(0),library=float(lib[ix].sum()),cells=len(ix),valid=valid.copy())
   bgs={name:aggregate(np.flatnonzero(outside&(lineage==name))) for name in bg_names if name!='Other_frozen_CAF' and np.any(outside&(lineage==name))}
   for (layer,target),ids in target_ids.items():
    ix=obs_ids.get_indexer(ids);assert (ix>=0).all();t=aggregate(ix)
    wanted=set(candidates.loc[candidates.dataset.eq(ds)&candidates.layer.eq(layer)&candidates.target.eq(target),'feature'])
    for j,f in enumerate(features):
     if f in wanted and valid[j]:target_detection.append(dict(dataset=ds,layer=layer,target=target,site=c.strict_site_group,donor_record=c.donor_record,feature=f,symbol=symbols[j],sample=sm,cells=t['cells'],positive_cells=t['positive'][j]))
    ob=layers[layer];other=ob.index[ob['sample'].eq(sm)&~ob.index.isin(ids)];oi=obs_ids.get_indexer(other);assert (oi>=0).all()
    backgrounds=dict(bgs)
    if len(oi):backgrounds['Other_frozen_CAF']=aggregate(oi)
    for bg,b in backgrounds.items():
     key=(ds,layer,target,c.donor_record,c.strict_site_group,bg)
     if key not in data:data[key]=dict(target=t.copy(),background=b.copy(),samples=[sm])
     else:
      dest=data[key]
      for role,extra in [('target',t),('background',b)]:
       for k in ['counts','positive','library','cells']:dest[role][k]=dest[role][k]+extra[k]
       dest[role]['valid']=dest[role]['valid']&extra['valid']
      dest['samples'].append(sm)
    checks.append(dict(dataset=ds,sample=sm,layer=layer,target=target,target_cells=len(ix),panel_features=len(features),source_cells=src.n_obs))
   src.file.close();print('SURFACE_RAW',ds,sm,len(features),len(target_ids),flush=True)
  for (sds,layer,target,donor,site,bg),d in data.items():
   if sds!=ds:continue
   t=d['target'];b=d['background'];valid=t['valid']&b['valid']
   delta=np.log2(t['counts']/t['library']*1e4+.1)-np.log2(b['counts']/b['library']*1e4+.1);delta[~valid]=np.nan
   wanted=set(candidates.loc[candidates.dataset.eq(ds)&candidates.layer.eq(layer)&candidates.target.eq(target),'feature'])
   for j,f in enumerate(features):
    if f not in wanted:continue
    source_records.append(dict(dataset=ds,layer=layer,target=target,donor_record=donor,site=site,background=bg,feature=f,symbol=symbols[j],effect=delta[j],target_cells=t['cells'],background_cells=b['cells'],target_positive_fraction=t['positive'][j]/t['cells'],background_positive_fraction=b['positive'][j]/b['cells'],samples=';'.join(d['samples'])))
 effects=pd.DataFrame(source_records);save(effects,'surface_background_patient_site_effects.tsv.gz')
 combined=effects.groupby(['dataset','layer','target','donor_record','background','feature','symbol'],as_index=False).agg(effect=('effect','mean'),target_positive_fraction=('target_positive_fraction','mean'),background_positive_fraction=('background_positive_fraction','mean'))
 combined['site']='Both_patient_equal_sites'
 tables=[]
 for (ds,layer,target,site,bg),z in pd.concat([effects,combined],ignore_index=True).groupby(['dataset','layer','target','site','background']):
  mat=z.pivot(index='donor_record',columns='feature',values='effect');arr=mat.to_numpy();st=m.stats(np.zeros_like(arr),arr,True);st['feature']=mat.columns
  symbols=z.drop_duplicates('feature').set_index('feature').symbol;st['symbol']=st.feature.map(symbols)
  for k,v in dict(dataset=ds,layer=layer,target=target,site=site,background=bg).items():st[k]=v
  st['family']='|'.join([ds,layer,target,site,'all_surface_candidates_all_backgrounds'])
  st['target_positive_fraction_median']=st.feature.map(z.groupby('feature').target_positive_fraction.median())
  st['background_positive_fraction_median']=st.feature.map(z.groupby('feature').background_positive_fraction.median())
  tables.append(st)
 allstats=m.fdr(pd.concat(tables,ignore_index=True));allstats=allstats.rename(columns={'mean_lung_minus_bone':'mean_target_minus_background','median_lung_minus_bone':'median_target_minus_background','bone_N':'N','lung_N':'N_check'})
 save(allstats,'surface_background_statistics.tsv')
 detection=pd.DataFrame(target_detection).groupby(['dataset','layer','target','site','donor_record','feature','symbol'],as_index=False).agg(cells=('cells','sum'),positive_cells=('positive_cells','sum'))
 detection['RNA_positive_fraction']=detection.positive_cells/detection.cells
 save(detection,'surface_target_patient_RNA_detection.tsv')
 reviews=[]
 for r in candidates.to_dict('records'):
  site='Lung' if r['direction']=='lung_higher' else 'Bone'
  z=allstats[allstats.dataset.eq(r['dataset'])&allstats.layer.eq(r['layer'])&allstats.target.eq(r['target'])&allstats.feature.eq(r['feature'])&allstats.site.eq(site)]
  det=detection[detection.dataset.eq(r['dataset'])&detection.layer.eq(r['layer'])&detection.target.eq(r['target'])&detection.feature.eq(r['feature'])&detection.site.eq(site)]
  r.update(background_review_site=site,background_types_evaluable=int(z.mean_target_minus_background.notna().sum()),background_types_target_higher=int(z.mean_target_minus_background.gt(0).sum()),background_types_nominally_higher=int((z.mean_target_minus_background.gt(0)&z.P.lt(.05)).sum()),background_types_target_lower=';'.join(z.loc[z.mean_target_minus_background.lt(0),'background']),RNA_target_detection_median=float(det.RNA_positive_fraction.median()) if len(det) else np.nan,flow_status='RNA candidate only; extracellular epitope/clone and protein background require confirmation')
  reviews.append(r)
 save(pd.DataFrame(reviews),'surface_marker_evidence_review.tsv')
 save(pd.DataFrame(checks),'surface_raw_sample_checks.tsv')
 qa=dict(status='PASS',candidates=len(candidates),unique_genes=candidates.symbol.nunique(),raw_sample_target_checks=len(checks),patient_site_effect_rows=len(effects),background_statistics=len(allstats),all_candidates_retained=True,all_source_cell_IDs_matched=True,no_flow_purity_claim=True)
 dump(O/'QA/SURFACE_ACCEPTANCE.json',qa);dump(O/'SURFACE_COMPLETE.json',qa);print(json.dumps(qa),flush=True)
if __name__=='__main__':
 try:main()
 except BaseException:
  traceback.print_exc();(O/'surface.exit_code').write_text('1\n');raise
 else:(O/'surface.exit_code').write_text('0\n')
