from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='2'
import json,hashlib,datetime,traceback,itertools,math,warnings,shutil,importlib.util
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
O=P/'phase73_F134_bone_lung_20261009/strict_v2'
A=P/'phase41_full_transcriptome_setup_20261001'
CLASSES=['F1','F3','F4']; SEED=20261009
def read(f,**kw):return pd.read_csv(f,sep='\t',low_memory=False,**kw)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def dump(f,x):f.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def save(d,f):d.to_csv(f,sep='\t',index=False,na_rep='NA')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as s:
  for b in iter(lambda:s.read(8388608),b''):h.update(b)
 return h.hexdigest()
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
def mean(x):
 out=np.full(x.shape[1],np.nan);n=np.isfinite(x).sum(0)
 np.divide(np.nansum(x,axis=0),n,out=out,where=n>0);return out
weights_cache={}
def test(x,y,paired):
 """Mean effect, exact or registered Monte Carlo patient test; featurewise missingness."""
 nfeat=x.shape[1];pv=np.full(nfeat,np.nan);status=np.full(nfeat,'insufficient_replication',object)
 eff=mean(y)-mean(x);n1=np.isfinite(x).sum(0);n2=np.isfinite(y).sum(0)
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',RuntimeWarning)
  med=(np.nanmedian(y-x,axis=0) if paired and len(x) else np.nanmedian(y,axis=0)-np.nanmedian(x,axis=0)) if len(x) and len(y) else np.full(nfeat,np.nan)
 pos=np.sum((y-x)>0,axis=0) if paired else np.full(nfeat,np.nan)
 neg=np.sum((y-x)<0,axis=0) if paired else np.full(nfeat,np.nan)
 if not len(x) or not len(y):return eff,med,pv,n1,n2,pos,neg,status
 valid=np.isfinite(x)&np.isfinite(y) if paired else np.vstack([np.isfinite(x),np.isfinite(y)])
 patterns,inverse=np.unique(valid.T,axis=0,return_inverse=True)
 for j,pat in enumerate(patterns):
  gi=np.flatnonzero(inverse==j)
  if paired:
   nx=int(pat.sum());ny=nx;v=(y-x)[pat][:,gi];obs=np.abs(v.mean(0)) if nx else np.zeros(len(gi))
   n1[gi]=nx;n2[gi]=ny
  else:
   p1=pat[:len(x)];p2=pat[len(x):];nx=int(p1.sum());ny=int(p2.sum());v=np.vstack([x[p1][:,gi],y[p2][:,gi]])
   obs=np.abs(v[nx:].mean(0)-v[:nx].mean(0)) if nx and ny else np.zeros(len(gi))
  if nx<2 or ny<2:continue
  nonzero=np.any(v!=0,axis=0) if paired else np.any(v!=v[0],axis=0)
  status[gi[~nonzero]]='all_zero_effect_or_no_variation'
  gi=gi[nonzero];v=v[:,nonzero];obs=obs[nonzero]
  if not len(gi):continue
  key=(paired,nx,ny)
  if key not in weights_cache:
   if paired:
    exact=nx<=15;w=np.array(list(itertools.product([-1.,1.],repeat=nx))) if exact else np.random.default_rng(SEED).choice([-1.,1.],size=(65536,nx));w=w/nx
   else:
    exact=math.comb(nx+ny,ny)<=65536
    combinations=list(itertools.combinations(range(nx+ny),ny)) if exact else [np.random.default_rng(SEED+i).choice(nx+ny,ny,replace=False) for i in range(65536)]
    w=np.full((len(combinations),nx+ny),-1/nx)
    for ii,ix in enumerate(combinations):w[ii,list(ix)]=1/ny
   weights_cache[key]=(w,exact)
  w,exact=weights_cache[key]
  for st in range(0,len(gi),256):
   null=np.abs(w@v[:,st:st+256]);hits=(null>=obs[None,st:st+256]-1e-12).sum(0)
   pv[gi[st:st+256]]=hits/len(w) if exact else (hits+1)/(len(w)+1)
  status[gi]=('exact_' if exact else 'MC65536_plus_observed_')+('patient_signflip' if paired else 'patient_label_permutation')
 return eff,med,pv,n1,n2,pos,neg,status
def stats(x,y,paired):
 e,m,p,n1,n2,pos,neg,s=test(x,y,paired)
 return pd.DataFrame(dict(bone_N=n1,lung_N=n2,mean_lung_minus_bone=e,median_lung_minus_bone=m,P=p,positive_patients=pos,negative_patients=neg,test=s))
def fdr(d):
 d=d.copy();d['FDR']=np.nan
 for fam,g in d.groupby('family',sort=False):
  q=multipletests(g.P.fillna(1).to_numpy(),method='fdr_bh')[1];d.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 d['nominal']=d.P.lt(.05);d['FDR_significant']=d.FDR.lt(.05);return d

def main():
 assert not (O/'CORE_COMPLETE.json').exists()
 for sub in ['code','tables','objects','resources','QA','figures','delivery','logs']:(O/sub).mkdir(parents=True,exist_ok=True)
 paths={
  'metadata':P/'phase66_doublet_repair_DSC1_20261005/tables/cells_repaired_with_confidence.tsv.gz',
  'clinical':P/'phase71_allcell_AB1B2_grouping_20261009/tables/sample_group_manifest.tsv',
  'supported':P/'phase30_two_layer_classification_20260926/tables/supported_classified_cells.tsv.gz',
  'expanded':P/'phase30_two_layer_classification_20260926/tables/expanded_classified_cells.tsv.gz',
  'states':P/'phase42_F3F4_identity_states_20261001/tables/all_F3F4_internal_cell_states.tsv.gz',
  'counts':A/'objects/all_candidate_full_genes.h5ad',
  'presence':A/'objects/sample_feature_presence.npz',
  'canonical':P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',
  'programs':P/'phase58_F1F5_round1_20261004/resources/programs_frozen.json',
  'libraries':P/'phase58_F1F5_round1_20261004/resources/libraries_canonical.json',
  'surface':P/'phase44_patient_mechanism_20261001/resources/intercell/receptor_roles.tsv',
  'legacy_surface':P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json',
  'adjudicated_sites':P/'phase44_patient_mechanism_20261001/clinical_site_review/tables/clinical_manifest_with_strict_site_group.tsv',
 }
 dump(O/'resources/INPUT_MANIFEST.json',{k:dict(path=str(f),sha256=sha(f),bytes=f.stat().st_size) for k,f in paths.items()})
 spec=importlib.util.spec_from_file_location('current_metadata',P/'phase66_doublet_repair_DSC1_20261005/code/current_metadata.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 cells=mod.load_current_cells().set_index('cell_id',drop=False)
 cm=read(paths['clinical'],keep_default_na=False).set_index('sample')
 adjudicated=read(paths['adjudicated_sites'],keep_default_na=False).set_index('sample')
 assert set(cm.index)==set(adjudicated.index)
 assert (cm.site==adjudicated.loc[cm.index,'site']).all()
 corrections=cm.loc[cm.strict_site_group.ne(adjudicated.loc[cm.index,'strict_site_group']),['site','strict_site_group']].copy()
 corrections['adjudicated_site']=adjudicated.loc[corrections.index,'strict_site_group']
 save(corrections.reset_index(),O/'QA/site_adjudication_corrections.tsv')
 cm['strict_site_group']=adjudicated.loc[cm.index,'strict_site_group']
 cells['strict_site']=cells['sample'].map(cm.strict_site_group)
 layers={k:read(paths[k],keep_default_na=False).set_index('cell_id') for k in ['supported','expanded']}
 assert len(layers['supported'])==1224 and len(layers['expanded'])==8174
 assert layers['expanded'].class_id.to_dict()==cells.class_id.to_dict()
 frozen_states=read(paths['states'],keep_default_na=False)
 protocol={
  'frozen_at':now(),'reference_report':'User supplied 20260928 integrated CAF/myCAF PDF, pages 2,4,9,20',
  'labels':'Original supported1224 and expanded8174 class memberships exactly preserved; F4 absent from supported layer',
  'two_comparison_levels':['patient-level captured proportions','within-class/state expression and function'],
  'targets':'F1 F3 F4 and fixed internal subdivisions',
  'internal_labels':'Expanded Phase66 subgroup_id; supported F3 Phase42 fixed states; supported F1 original within-layer cluster labels (not new subclustering)',
  'site':'Accepted Phase44 clinical_site_review anatomical adjudication; Paraspinal remains Unknown_or_mixed; original clinical group labels unchanged',
  'design':'Source-specific; same donor with target/denominator observed at both sites analyzed as paired, removed from independent contrast for that endpoint; no cross-source patient pooling',
  'proportion_denominators':['full source-called cells for parent classes','all cells in frozen layer for parent classes','parent class cells for internal states'],
  'proportion_policy':'Sum numerators and denominators within patient-site; numerator zero retained; zero denominator excluded as undefined',
  'expression':'Raw counts patient-site pseudobulk, full original library denominator; log2(CP10k+0.1); no Harmony expression values',
  'programs':'Frozen 1867 sets; mean per-cell log1p(CP10k), then mean within patient-site and measured canonical set members; no extra coverage gate',
  'missing_features':'Sample-feature presence; absent measurement remains NA, not zero',
  'testing':'Exact paired sign-flip N<=15; independent label permutation exact if <=65536 assignments otherwise fixed65536 MC+observed; fewer than2 patients per independent arm or paired N<2 => no P',
  'gene_families':'layer|source|design|main_classes genes; internal genes separately per parent; all source feature hypotheses retained',
  'program_families':'Same class partitions; programs separate from genes',
  'proportion_families':'layer|source|design|denominator|main or internal-parent',
  'ORA':'Canonical full evaluable gene background; raw P<0.05 up/down lists; each library and direction separate BH family; no ORA if list empty',
  'surface':'Nominal site-differential genes intersect full frozen receptor-role plasma membrane transmembrane/peripheral annotation OR original prespecified16 surface panel; no filtering by previous F-class marker discovery; protein topology and flow feasibility separate',
  'amendment':'v1 used a superseded anatomical field and a derivative108-gene surface evidence subset. v2 restores accepted site adjudication and full annotation coverage. v1 preserved as superseded; no outcome-driven numerical gate added.',
  'sensitivity':'Remove repaired doublet-positive cells with same labels, contrasts and families; evaluate changes',
  'scope_limits':'No new clustering or relabeling; no new numerical cell-entry cutoff; source/clinical confounding not erased by grouping; overlapping supported/expanded layers not independent validation',
 }
 dump(O/'resources/PROTOCOL_FROZEN.json',protocol)
 save(cm.reset_index(),O/'resources/clinical_sample_manifest.tsv')
 save(cm.loc[~cm.strict_site_group.isin(['Bone','Lung'])].reset_index(),O/'tables/excluded_unknown_mixed_samples.tsv')
 a=ad.read_h5ad(paths['counts']);assert a.shape==(8174,63398)
 presence=sparse.load_npz(paths['presence']).tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict()
 features=a.var_names.astype(str).to_numpy();symbols=a.var.symbol.astype(str).to_numpy()
 assert set(a.obs_names)==set(cells.index)
 assert np.allclose(np.asarray(a.X.sum(1)).ravel(),cells.loc[a.obs_names,'full_library_counts'].to_numpy(float))
 canonical=read(paths['canonical'],keep_default_na=False)
 programs=json.loads(paths['programs'].read_text())['sets'];program_names=list(programs)
 libs=json.loads(paths['libraries'].read_text())
 surface=read(paths['surface'],keep_default_na=False);surface_symbols=set(surface.loc[truth(surface.plasma_membrane_transmembrane)|truth(surface.plasma_membrane_peripheral),'genesymbol'])
 legacy_surface=set(json.loads(paths['legacy_surface'].read_text())['surface_candidates']);surface_symbols|=legacy_surface
 save(surface[surface.genesymbol.isin(surface_symbols)],O/'resources/full_surface_role_evidence.tsv')
 dump(O/'resources/surface_screen_universe.json',{'genes':sorted(surface_symbols),'legacy_prespecified':sorted(legacy_surface),'count':len(surface_symbols)})
 allmeta=[];targets_by_layer={};registry=[]
 for layer,old in layers.items():
  ob=cells.loc[old.index].copy();ob['analysis_class']=old.class_id
  ob['analysis_internal']='';ob['layer']=layer
  if layer=='expanded':ob['analysis_internal']=ob.subgroup_id
  else:
   sm=frozen_states[frozen_states.layer.eq('supported')&frozen_states.parent_class.eq('F3')].set_index('cell_id').state_id
   ob.loc[ob.analysis_class.eq('F3'),'analysis_internal']=sm.reindex(ob.index[ob.analysis_class.eq('F3')]).to_numpy()
   ob.loc[ob.analysis_class.eq('F1'),'analysis_internal']='supported_F1_originalC'+old.loc[ob.analysis_class.eq('F1'),'cluster'].astype(str)
  assert ob.loc[ob.analysis_class.isin(['F1','F3']),'analysis_internal'].notna().all()
  targets=[]
  for parent in CLASSES:
   targets.append(dict(target=parent,parent=parent,level='main'))
   for state in sorted(ob.loc[ob.analysis_class.eq(parent),'analysis_internal'].unique()):
    if state:targets.append(dict(target=state,parent=parent,level='internal'))
  targets_by_layer[layer]=targets
  for t in targets:
   mask=ob.analysis_class.eq(t['target']) if t['level']=='main' else ob.analysis_internal.eq(t['target'])
   registry.append(dict(layer=layer,**t,cells=int(mask.sum()),original_label_preserved=True))
  allmeta.append(ob)
  save(ob.reset_index(drop=True),O/'tables'/f'{layer}_fixed_cells_with_sites.tsv.gz')
 save(pd.DataFrame(registry),O/'tables/frozen_target_registry.tsv')
 # Patient-site fractions, including zero target counts in every eligible denominator.
 prop_values=[];prop_statistics=[];prop_designs=[]
 for sensitivity in ['all_cells','exclude_doublets']:
  for layer,obfull in zip(layers,allmeta):
   ob=obfull if sensitivity=='all_cells' else obfull[~truth(obfull.doublet_positive_any)]
   for ds,ms in cm[cm.strict_site_group.isin(['Bone','Lung'])].groupby('dataset'):
    for target in targets_by_layer[layer]:
     denoms=['full_called','layer'] if target['level']=='main' else ['parent']
     for denominator in denoms:
      vals={'Bone':{},'Lung':{}}
      for (donor,site),sg in ms.groupby(['donor_record','strict_site_group']):
       z=ob[ob['sample'].isin(sg.index)]
       numerator=int((z.analysis_class.eq(target['target']) if target['level']=='main' else z.analysis_internal.eq(target['target'])).sum())
       denom=int(sg.full_called_cells.astype(int).sum()) if denominator=='full_called' else (len(z) if denominator=='layer' else int(z.analysis_class.eq(target['parent']).sum()))
       value=numerator/denom if denom else np.nan
       prop_values.append(dict(sensitivity=sensitivity,layer=layer,dataset=ds,**target,denominator=denominator,donor_record=donor,site=site,numerator=numerator,denominator_cells=denom,fraction=value,samples=';'.join(sg.index)))
       if denom:vals[site][donor]=value
      paired=sorted(set(vals['Bone'])&set(vals['Lung']))
      for design in ['paired','independent']:
       ids1=paired if design=='paired' else sorted(set(vals['Bone'])-set(paired));ids2=paired if design=='paired' else sorted(set(vals['Lung'])-set(paired))
       x=np.array([vals['Bone'][d] for d in ids1]).reshape(-1,1);y=np.array([vals['Lung'][d] for d in ids2]).reshape(-1,1)
       st=stats(x,y,design=='paired');st['bone_mean_fraction']=mean(x);st['lung_mean_fraction']=mean(y)
       meta=dict(sensitivity=sensitivity,layer=layer,dataset=ds,**target,denominator=denominator,design=design)
       for k,v in meta.items():st[k]=v
       st['family']='|'.join([sensitivity,layer,ds,design,denominator,'main' if target['level']=='main' else target['parent']+'_internal'])
       prop_statistics.append(st)
       prop_designs.append(dict(**meta,bone_patients=';'.join(ids1),lung_patients=';'.join(ids2)))
 save(pd.DataFrame(prop_values),O/'tables/patient_site_proportions.tsv.gz')
 prop=fdr(pd.concat(prop_statistics,ignore_index=True));save(prop,O/'tables/all_proportion_statistics.tsv');save(prop[prop.nominal],O/'tables/nominal_proportions.tsv')
 save(pd.DataFrame(prop_designs),O/'tables/proportion_test_patient_roster.tsv')
 print('PROPORTIONS_COMPLETE',len(prop),flush=True)
 contexts=[];expr_roster=[];gene_paths=[];program_frames=[];coverage=[];expression_effects=[]
 # Sparse matrix of per-cell log1p full-library normalized raw expression.
 logx=a.X.astype(float).tocsr().copy();den=np.asarray(a.X.sum(1)).ravel();logx=sparse.diags(1e4/den)@logx;logx.data=np.log1p(logx.data)
 source_mats={}
 for ds,sg in cm.groupby('dataset'):
  gi=np.flatnonzero(np.asarray(presence[[sidx[s] for s in sg.index]].sum(0)).ravel()>0)
  localmap={features[j]:k for k,j in enumerate(gi)}
  canon=canonical[canonical.dataset.eq(ds)];canonical_pos={r.symbol:localmap[r.feature] for r in canon.itertuples() if r.feature in localmap}
  rr=[];cc=[]
  for k,name in enumerate(program_names):
   mapped=[canonical_pos[s] for s in programs[name] if s in canonical_pos]
   rr.extend(mapped);cc.extend([k]*len(mapped))
   coverage.append(dict(dataset=ds,program=name,definition_members=len(programs[name]),source_measured_members=len(mapped),missing_members=';'.join(s for s in programs[name] if s not in canonical_pos)))
  mat=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(gi),len(program_names)))
  source_mats[ds]=(gi,mat,canon)
 for layer,obfull in zip(layers,allmeta):
  for ds in cm.dataset.unique():
   gi,progmat,canon=source_mats[ds];fe=features[gi];sy=symbols[gi]
   source_frames=[]
   for sensitivity in ['all_cells','exclude_doublets']:
    ob=obfull if sensitivity=='all_cells' else obfull[~truth(obfull.doublet_positive_any)]
    ob=ob[ob.dataset.eq(ds)&ob.strict_site.isin(['Bone','Lung'])]
    for target in targets_by_layer[layer]:
     z=ob[ob.analysis_class.eq(target['target'])] if target['level']=='main' else ob[ob.analysis_internal.eq(target['target'])]
     vals={'Bone':{},'Lung':{}};pvals={'Bone':{},'Lung':{}}
     for (donor,site),q in z.groupby(['donor_record','strict_site']):
      ix=a.obs_names.get_indexer(q.index);assert (ix>=0).all()
      raw=np.asarray(a.X[ix][:,gi].sum(0)).ravel();library=float(q.full_library_counts.sum())
      samples=q['sample'].unique().tolist();valid=np.asarray(presence[[sidx[s] for s in samples]][:,gi].sum(0)).ravel()==len(samples)
      v=np.log2(raw/library*1e4+0.1);v[~valid]=np.nan;vals[site][donor]=v
      lm=np.asarray(logx[ix][:,gi].mean(0)).ravel();nm=np.asarray(valid.astype(float)@progmat).ravel()
      numer=np.asarray(np.where(valid,lm,0)@progmat).ravel();pv=np.full(len(program_names),np.nan);np.divide(numer,nm,out=pv,where=nm>0);pvals[site][donor]=pv
      expr_roster.append(dict(sensitivity=sensitivity,layer=layer,dataset=ds,**target,donor_record=donor,site=site,cells=len(q),samples=';'.join(samples),library_UMIs=library))
     shared=sorted(set(vals['Bone'])&set(vals['Lung']))
     for design in ['paired','independent']:
      ids1=shared if design=='paired' else sorted(set(vals['Bone'])-set(shared));ids2=shared if design=='paired' else sorted(set(vals['Lung'])-set(shared))
      key='__'.join([sensitivity,layer,ds,target['target'],design])
      meta=dict(sensitivity=sensitivity,layer=layer,dataset=ds,**target,design=design,context_id=key)
      contexts.append(dict(**meta,bone_N=len(ids1),lung_N=len(ids2),bone_patients=';'.join(ids1),lung_patients=';'.join(ids2),status='testable' if min(len(ids1),len(ids2))>=2 else ('direction_only' if ids1 and ids2 else 'one_or_both_groups_absent')))
      if not ids1 or not ids2:continue
      x=np.vstack([vals['Bone'][d] for d in ids1]);y=np.vstack([vals['Lung'][d] for d in ids2])
      frame=stats(x,y,design=='paired');frame['feature']=fe;frame['symbol']=sy
      for k,v in meta.items():frame[k]=v
      fam='|'.join([sensitivity,layer,ds,design,'main' if target['level']=='main' else target['parent']+'_internal'])
      frame['family']=fam+'|genes';source_frames.append(frame)
      ps=stats(np.vstack([pvals['Bone'][d] for d in ids1]),np.vstack([pvals['Lung'][d] for d in ids2]),design=='paired');ps['program']=program_names
      for k,v in meta.items():ps[k]=v
      ps['family']=fam+'|programs';program_frames.append(ps)
      np.savez_compressed(O/'objects'/(key+'.npz'),features=fe.astype(str),symbols=sy.astype(str),bone_patients=np.array(ids1,dtype=str),lung_patients=np.array(ids2,dtype=str),bone_expression=x,lung_expression=y,bone_programs=np.vstack([pvals['Bone'][d] for d in ids1]),lung_programs=np.vstack([pvals['Lung'][d] for d in ids2]),programs=np.array(program_names,dtype=str))
      print('EXPRESSION',key,len(ids1),len(ids2),int(frame.P.notna().sum()),flush=True)
   if source_frames:
    gd=fdr(pd.concat(source_frames,ignore_index=True));path=O/'tables'/f'{layer}__{ds}__all_gene_statistics.tsv.gz';save(gd,path);gene_paths.append(path)
    save(gd[gd.nominal],O/'tables'/f'{layer}__{ds}__nominal_genes.tsv.gz')
    for key,g in gd.groupby('context_id'):
     expression_effects.append(dict(context_id=key,rows=len(g),valid_tests=int(g.P.notna().sum()),nominal_up=int((g.nominal&g.mean_lung_minus_bone.gt(0)).sum()),nominal_down=int((g.nominal&g.mean_lung_minus_bone.lt(0)).sum()),FDR_significant=int(g.FDR_significant.sum())))
 save(pd.DataFrame(contexts),O/'tables/expression_comparison_registry.tsv')
 save(pd.DataFrame(expr_roster),O/'tables/expression_patient_roster.tsv')
 save(pd.DataFrame(coverage),O/'tables/program_coverage.tsv')
 save(pd.DataFrame(expression_effects),O/'tables/gene_result_summary.tsv')
 progs=fdr(pd.concat(program_frames,ignore_index=True));save(progs,O/'tables/all_program_statistics.tsv.gz');save(progs[progs.nominal],O/'tables/nominal_programs.tsv')
 # Canonical gene ORA from nominal patient-level contrasts. Libraries and candidates frozen before tests.
 ora=[];markerframes=[];sensitivityrows=[];clinical_targets=[]
 for path in gene_paths:
  d=read(path);ds=d.dataset.iloc[0];can=canonical[canonical.dataset.eq(ds)][['feature','symbol']].drop_duplicates('feature')
  for key,g in d.groupby('context_id',sort=False):
   meta={k:g[k].iloc[0] for k in ['sensitivity','layer','dataset','target','parent','level','design','context_id']}
   if meta['sensitivity']!='all_cells':continue
   canon=g.drop(columns='symbol').merge(can,on='feature',how='inner',validate='one_to_one');valid=canon[canon.P.notna()];background=set(valid.symbol)
   # Compare to the exact same source/target/design after doublet exclusion.
   senskey=key.replace('all_cells__','exclude_doublets__',1);sens=d[d.context_id.eq(senskey)]
   hit=canon[canon.P.lt(.05)&canon.symbol.isin(surface_symbols)].copy()
   if len(hit):
    ss=sens[['feature','mean_lung_minus_bone','P','FDR']].rename(columns={'mean_lung_minus_bone':'doublet_excluded_effect','P':'doublet_excluded_P','FDR':'doublet_excluded_FDR'})
    hit=hit.merge(ss,on='feature',how='left');hit['direction']=np.where(hit.mean_lung_minus_bone>0,'lung_higher','bone_higher');hit['surface_evidence']='full_frozen_membrane_annotation_or_original_surface_panel';hit['protein_flow_validated']=False;markerframes.append(hit)
   if not background:continue
   for direction in ['up','down']:
    selected=set(valid.loc[valid.P.lt(.05)&(valid.mean_lung_minus_bone.gt(0) if direction=='up' else valid.mean_lung_minus_bone.lt(0)),'symbol'])
    if not selected:continue
    for library,terms in libs.items():
     rows=[]
     for name,members in terms.items():
      measured=set(members)&background
      if not measured or len(measured)==len(background):continue
      overlap=selected&measured
      pp=float(hypergeom.sf(len(overlap)-1,len(background),len(measured),len(selected)))
      rows.append(dict(**meta,direction=direction,library=library,term=name,background_genes=len(background),selected_genes=len(selected),term_measured_genes=len(measured),overlap_genes=len(overlap),overlap_symbols=';'.join(sorted(overlap)),P=pp,family=key+'|ORA|'+library+'|'+direction))
     ora.extend(rows)
 if ora:
  od=fdr(pd.DataFrame(ora));save(od,O/'tables/all_ORA_statistics.tsv.gz');save(od[od.nominal],O/'tables/nominal_ORA.tsv')
 else:save(pd.DataFrame(columns=['context_id','term','P','FDR']),O/'tables/nominal_ORA.tsv')
 if markerframes:
  markers=pd.concat(markerframes,ignore_index=True)
  bg=read(P/'phase58_F1F5_round1_20261004/tables/surface_competition_statistics.tsv')
  evidence=[]
  for row in markers.itertuples():
   rr=bg[bg.dataset.eq(row.dataset)&bg.class_id.eq(row.parent)&bg.gene.eq(row.symbol)] if row.layer=='expanded' and row.level=='main' else bg.iloc[0:0]
   evidence.append(dict(context_id=row.context_id,feature=row.feature,existing_exact_parent_background_comparisons=len(rr),background_negative_or_nonsignificant=int(((rr.mean_effect<=0)|(rr.P>=.05)|rr.P.isna()).sum()) if len(rr) else np.nan,background_evidence_scope='existing expanded parent class only' if len(rr) else 'not_yet_tested_for_this_exact_population'))
  markers=markers.merge(pd.DataFrame(evidence),on=['context_id','feature'],validate='one_to_one')
  save(markers,O/'tables/nominal_surface_marker_candidates.tsv')
 else:save(pd.DataFrame(columns=['context_id','symbol','P','FDR']),O/'tables/nominal_surface_marker_candidates.tsv')
 # Raw source verification for all main-class cells at both sites in all studies with testable main expression.
 ctx=pd.DataFrame(contexts);vctx=ctx[(ctx.sensitivity=='all_cells')&(ctx.level=='main')&(ctx.status=='testable')]
 check_datasets=set(vctx.dataset);rawqa=[]
 for sm in sorted(set(cells.loc[cells.dataset.isin(check_datasets)&cells.strict_site.isin(['Bone','Lung']),'sample'])):
  source=ad.read_h5ad(P/'phase26_client_reanalysis_20260922/inputs'/sm/'full_called_counts.h5ad',backed='r')
  for layer,ob in zip(layers,allmeta):
   for cls in CLASSES:
    ids=ob.index[ob['sample'].eq(sm)&ob.analysis_class.eq(cls)]
    if not len(ids):continue
    ix=source.obs_names.get_indexer(ids);assert (ix>=0).all();raw=source.X[np.sort(ix),:].tocsr()
    ai=a.obs_names.get_indexer(ids);vi=a.var_names.get_indexer(source.var_names);assert (vi>=0).all()
    err=float(np.max(np.abs(np.asarray(raw.sum(0)).ravel()-np.asarray(a.X[ai][:,vi].sum(0)).ravel())));assert err==0
    rawqa.append(dict(layer=layer,class_id=cls,sample=sm,cells=len(ids),max_count_error=err))
  source.file.close()
 save(pd.DataFrame(rawqa),O/'QA/raw_source_count_verification.tsv')
 # Structural acceptance: full frozen memberships, unpaired donor disjointness and sample-size rules.
 assert all(not(set(r['bone_patients'].split(';'))&set(r['lung_patients'].split(';'))- {''}) for r in contexts if r['design']=='independent')
 assert not progs.loc[(progs.bone_N<2)|(progs.lung_N<2),'P'].notna().any()
 assert not prop.loc[(prop.bone_N<2)|(prop.lung_N<2),'P'].notna().any()
 qa=dict(at=now(),status='PASS',reference_class_counts={k:old.class_id.value_counts().to_dict() for k,old in layers.items()},source_cell_sets_preserved=True,strict_sites_only=True,raw_source_checks=len(rawqa),max_raw_count_error=0,contexts=len(contexts),testable_expression_contexts=int((ctx.status=='testable').sum()),proportion_rows=len(prop),program_rows=len(progs),gene_files=[str(f) for f in gene_paths],no_fabricated_P_for_insufficient_N=True)
 dump(O/'QA/CORE_ACCEPTANCE.json',qa);dump(O/'CORE_COMPLETE.json',qa)
 print('CORE_COMPLETE',json.dumps({k:v for k,v in qa.items() if k!='gene_files'},ensure_ascii=False),flush=True)

if __name__=='__main__':
 try:main()
 except BaseException:
  if O.exists():(O/'core.exit_code').write_text('1\n')
  traceback.print_exc();raise
 else:(O/'core.exit_code').write_text('0\n')
