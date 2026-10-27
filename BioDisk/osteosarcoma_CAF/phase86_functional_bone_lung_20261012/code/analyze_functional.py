from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,importlib.util,shutil,traceback
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase86_functional_bone_lung_20261012'
for d in ['resources','tables','objects','figures','logs','QA','delivery','state_backups']:(O/d).mkdir(parents=True,exist_ok=True)
DS=['SCPCP000017','SCPCP000023']
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def save(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
def dump(x,p):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
def adjust(df):
 df=df.copy();df['FDR']=np.nan
 for _,g in df.groupby('family'):
  ix=g.index[g.P.notna()]
  if len(ix):df.loc[ix,'FDR']=multipletests(df.loc[ix,'P'],method='fdr_bh')[1]
 df['nominal']=df.P.lt(.05);return df
spec=importlib.util.spec_from_file_location('core',P/'phase73_F134_bone_lung_20261009/strict_v2/code/analyze_bone_lung_v2.py')
core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.SEED=20261012
def stat(x,y,paired):
 z=core.stats(x,y,paired)
 if paired:z['mean_lung_minus_bone']=core.mean(y-x)
 return z

def main():
 assert not (O/'COMPUTATION_COMPLETE.json').exists()
 paths=dict(metadata=P/'phase73_F134_bone_lung_20261009/strict_v2/tables/expanded_fixed_cells_with_sites.tsv.gz',clinical=P/'phase81_allCAF_bone_lung_20261012/tables/clinical_sample_manifest.tsv',programs=P/'phase66_doublet_repair_DSC1_20261005/resources/PROGRAMS.json',counts=P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',canonical=P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',presence=P/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz',roles=P/'phase44_patient_mechanism_20261001/resources/intercell/receptor_roles.tsv',legacy_panels=P/'phase58_F1F5_round1_20261004/resources/programs_frozen.json',mechanism=P/'phase63_matrix_resistance_20261005/objects/mechanism_counts.npz')
 manifest={k:dict(path=str(v),sha256=sha(v)) for k,v in paths.items()};dump(manifest,O/'resources/INPUT_MANIFEST.json')
 cfg=json.loads(paths['programs'].read_text());panels=cfg['sets'].copy()
 extra=json.loads(paths['legacy_panels'].read_text())['sets']
 for key,label in [('hallmark|HALLMARK_TGF_BETA_SIGNALING','TGF_beta'),('hallmark|HALLMARK_IL6_JAK_STAT3_SIGNALING','IL6_JAK_STAT3'),('customer|iCAF_customer','iCAF_reference')]:panels[label]=extra[key]
 protocol=dict(frozen_at=datetime.now(timezone.utc).isoformat(),clinical_framework='Anatomical Bone versus Lung is the primary comparison; clinical groups are not data-source IDs; A/B1/B2 retained but not used to filter current scheme1',objective='Locate clinically differing functional CAF states and assess surface-marker proxies',primary_sources=DS,record_only_sources=['SCPCP000018','GSE152048','GSE318840'],all14_retained=True,primary_cells='Existing reviewed Fibro lineage, exclude repaired model-positive doublets; no new identity reassignment',sensitivity='All expanded candidates excluding model-positive doublets; raw/corrected scoring where available; minimum5 patient-site cells as recorded sensitivity only',programs=list(panels),scoring='Frozen mean measured-member log1p(CP10k); reuse independently checked existing 11 raw/corrected scores; three extra frozen panels raw only',primary_endpoint='Patient-site continuous program mean; paired and independent designs separate',high_fraction='Secondary operational upper-quartile state. Per source/scope/matrix/program weighted 75th percentile with equal donor weight and equal sites per donor, blind to site labels. Strictly above threshold. Not a newly validated discrete CAF subtype.',independence='All donors recorded at both anatomical sites in clinical manifest excluded from independent comparisons, including orphaned scope-specific observations; paired comparison only complete pairs',statistics='Patient exact permutation / signflip via audited engine; insufficient biological repetition no P; valid-P BH within source/scope/matrix/design/endpoint/minimum-cell family',surface_search='All existing plasma-membrane annotation plus prior panel and all14; within-patient/site top versus bottom score quartiles; minimum5 cells per arm; mean donor effects after equal site weighting; exact signflip over donors; tested genes in defining program skipped as circular; family per source/scope/program across all noncomponent assayed marker genes',surface_effect='Patient high-vs-low pseudo-bulk log2(CP10k+0.1), not site effect; marker association exploratory and separate from functional site test',no_claim='Inflammation/TGF expression is not proof of immunosuppression; ECM transcription is not measured deposition, wrapping or drug resistance; RNA detection is not flow positivity',nextstep='No resolution tuning or outcome-dependent panel selection')
 dump(protocol,O/'resources/PROTOCOL_FROZEN.json');dump(panels,O/'resources/PROGRAMS.json')
 stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
 for n in ['CAF_NEXT_STAGE_CURRENT.json','CAF_RESUME_STATUS_CURRENT.json','CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json']:
  p=P/n;shutil.copy2(p,O/'state_backups'/f'{stamp}_{n}');v=json.loads(p.read_text());v['clinical_comparison_framework']=protocol['clinical_framework'];v['functional_bone_lung_analysis']=dict(phase=O.name,status='running',protocol=str(O/'resources/PROTOCOL_FROZEN.json'))
  if n!='CAF_FUNCTIONAL_STATE_OBJECTIVE_CURRENT.json':v['active_task']=dict(phase=O.name,status='running',title='Bone-Lung functional-state and surface-marker analysis')
  dump(v,p)
 ob=read(paths['metadata'],keep_default_na=False).set_index('cell_id',drop=False)
 loader_spec=importlib.util.spec_from_file_location('metadata_current',P/'phase66_doublet_repair_DSC1_20261005/code/current_metadata.py');loader=importlib.util.module_from_spec(loader_spec);loader_spec.loader.exec_module(loader)
 current=loader.load_current_cells().set_index('cell_id');assert set(current.index)==set(ob.index)
 assert np.array_equal(truth(ob.doublet_positive_any),current.loc[ob.index,'doublet_positive_any'].to_numpy(bool))
 cm=read(paths['clinical'],keep_default_na=False).set_index('sample');assert (ob.strict_site==ob['sample'].map(cm.strict_site_group)).all()
 clinical=cm[cm.dataset.isin(DS)].reset_index();save(clinical,'clinical_manifest.tsv')
 save(clinical.groupby(['dataset','strict_site_group','provider_disease_timing','treatment'],dropna=False).agg(samples=('sample','size'),donors=('donor_record','nunique')).reset_index(),'clinical_context_availability.tsv')
 a=ad.read_h5ad(paths['counts']);assert a.shape==(8174,63398)
 ob=ob.loc[a.obs_names].copy();ob.index.name='cell_index';raw=a.X.tocsr();library=np.asarray(raw.sum(1)).ravel();assert np.allclose(library,ob.full_library_counts.astype(float))
 logx=sparse.diags(1e4/library)@raw;logx.data=np.log1p(logx.data);logx=logx.tocsr()
 can=read(paths['canonical'],keep_default_na=False);fp={f:i for i,f in enumerate(a.var_names)}
 presence=sparse.load_npz(paths['presence']).tocsr();sidx=read(P/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict()
 maps={ds:{r.symbol:fp[r.feature] for r in can[can.dataset.eq(ds)].itertuples() if r.feature in fp} for ds in DS}
 score={};coverage=[]
 for pr in cfg['sets']:
  for matrix in ['raw','corrected']:score[matrix,pr]=ob[matrix+'__'+pr].astype(float).to_numpy()
 for pr in ['TGF_beta','IL6_JAK_STAT3','iCAF_reference']:score['raw',pr]=np.full(len(ob),np.nan)
 for sm,q in ob[ob.dataset.isin(DS)].groupby('sample'):
  ds=q.dataset.iloc[0];ix=ob.index.get_indexer(q.index);present=presence.getrow(sidx[sm]).toarray().ravel().astype(bool)
  for pr,members in panels.items():
   js=[maps[ds][g] for g in members if g in maps[ds] and present[maps[ds][g]]]
   coverage.append(dict(dataset=ds,sample=sm,program=pr,defined=len(members),measured=len(js)))
   if pr not in cfg['sets'] and js:score['raw',pr][ix]=np.asarray(logx[ix][:,js].mean(1)).ravel()
 save(pd.DataFrame(coverage),'program_coverage.tsv')
 # Independent reconstruction of all inherited scores from their frozen mechanism arrays.
 mz=np.load(paths['mechanism']);order=pd.Index(mz['cells'].astype(str)).get_indexer(ob.index);assert (order>=0).all();mg={g:i for i,g in enumerate(mz['genes'].astype(str))};ass=mz['assayed'][order].astype(bool)
 recon_errors=[]
 for matrix in ['raw','corrected']:
  cx=mz[matrix][order];den=mz[matrix+'_library'][order];vals=np.log1p(cx/den[:,None]*1e4);vals[~ass]=np.nan
  for pr,members in cfg['sets'].items():
   js=[mg[g] for g in members];reconstructed=np.nanmean(vals[:,js],axis=1);err=float(np.nanmax(np.abs(reconstructed-score[matrix,pr])));recon_errors.append(dict(matrix=matrix,program=pr,max_error=err));assert err<1e-6,(matrix,pr,err)
 save(pd.DataFrame(recon_errors),'score_reconstruction_QA.tsv')
 keep=ob.dataset.isin(DS)&ob.strict_site.isin(['Bone','Lung'])&~truth(ob.doublet_positive_any)
 scopes={'Fibro_reviewed':keep&ob.review_lineage_code.eq('Fibro'),'expanded_candidates':keep}
 cc=ob[['cell_id','dataset','sample','donor_record','strict_site','analysis_class','subgroup_id','review_lineage_code','doublet_positive_any','full_library_counts','detected_features','mitochondrial_fraction']].copy()
 for (matrix,pr),vals in score.items():cc[matrix+'__'+pr]=vals
 save(cc.loc[keep].reset_index(drop=True),'cell_scores_and_original_labels.tsv.gz')
 patients=[];thresholds=[];statecomp=[];quality=[];highs={}
 for scope,mask in scopes.items():
  for ds in DS:
   q=ob[mask&ob.dataset.eq(ds)];ix=ob.index.get_indexer(q.index)
   donor_sites=q[['donor_record','strict_site']].drop_duplicates().groupby('donor_record').size().to_dict()
   group_counts=q.groupby(['donor_record','strict_site']).size().to_dict()
   weights=np.array([1/(donor_sites[d]*group_counts[d,s]) for d,s in zip(q.donor_record,q.strict_site)])
   for (matrix,pr),vals in score.items():
    v=vals[ix];valid=np.isfinite(v);ss=np.argsort(v[valid]);vv=v[valid][ss];ww=weights[valid][ss];cut=float(vv[np.searchsorted(np.cumsum(ww),.75*ww.sum())]) if len(vv) else np.nan
    highs[scope,ds,matrix,pr]=cut;thresholds.append(dict(scope=scope,dataset=ds,matrix=matrix,program=pr,threshold=cut,rule='patient-balanced pooled 75th percentile; strict >'))
    tmp=q[['donor_record','strict_site','analysis_class','subgroup_id']].copy();tmp['score']=v;tmp['high']=np.where(valid,(v>cut).astype(float),np.nan)
    for (donor,site),g in tmp.groupby(['donor_record','strict_site']):
     gg=q.loc[g.index];patients.append(dict(scope=scope,dataset=ds,matrix=matrix,program=pr,donor_record=donor,site=site,cells=len(g),measured_cells=int(g.score.notna().sum()),mean_score=g.score.mean(),high_fraction=g.high.mean(),high_cells=int(g.high.fillna(0).sum()),detected_features_mean=pd.to_numeric(gg.detected_features).mean(),mitochondrial_fraction_mean=pd.to_numeric(gg.mitochondrial_fraction).mean(),samples=';'.join(sorted(gg['sample'].unique()))))
    for (cl,st),g in tmp.groupby(['analysis_class','subgroup_id']):statecomp.append(dict(scope=scope,dataset=ds,matrix=matrix,program=pr,analysis_class=cl,subgroup_id=st,cells=len(g),patients=g.donor_record.nunique(),mean_score=g.score.mean(),high_cells=int(g.high.fillna(0).sum()),high_fraction=g.high.mean()))
   print('SCORED',scope,ds,len(q),q.donor_record.nunique(),flush=True)
 pat=pd.DataFrame(patients);save(pat,'patient_site_programs.tsv');save(pd.DataFrame(thresholds),'state_thresholds.tsv');save(pd.DataFrame(statecomp),'functional_state_original_cluster_crosswalk.tsv')
 comparisons=[];roster=[]
 for (scope,ds,matrix,pr),g in pat.groupby(['scope','dataset','matrix','program'],sort=False):
  clin=clinical[(clinical.dataset==ds)&clinical.strict_site_group.isin(['Bone','Lung'])];paired_clin=set(clin.groupby('donor_record').strict_site_group.nunique().loc[lambda x:x==2].index)
  for minimum in [1,5]:
   gg=g[g.cells>=minimum];b=gg[gg.site=='Bone'].set_index('donor_record');l=gg[gg.site=='Lung'].set_index('donor_record')
   ids_pair=sorted(set(b.index)&set(l.index)&paired_clin)
   for design in ['independent','paired']:
    ib=ids_pair if design=='paired' else sorted(set(b.index)-paired_clin);il=ids_pair if design=='paired' else sorted(set(l.index)-paired_clin)
    for endpoint in ['mean_score','high_fraction']:
     x=b.loc[ib,endpoint].to_numpy(float).reshape(-1,1);y=l.loc[il,endpoint].to_numpy(float).reshape(-1,1)
     st=stat(x,y,design=='paired').iloc[0].to_dict();effect=st['mean_lung_minus_bone']
     if design=='independent':
      loo=[float(np.mean(y)-np.mean(np.delete(x,i,axis=0))) for i in range(len(x)) if len(x)>1 and len(y)]+[float(np.mean(np.delete(y,i,axis=0))-np.mean(x)) for i in range(len(y)) if len(y)>1 and len(x)]
     else:loo=[float(np.mean(np.delete(y-x,i,axis=0))) for i in range(len(x)) if len(x)>1]
     st.update(scope=scope,dataset=ds,matrix=matrix,program=pr,minimum_cells=minimum,design=design,endpoint=endpoint,bone_mean=float(x.mean()) if x.size else np.nan,lung_mean=float(y.mean()) if y.size else np.nan,bone_cells=int(b.loc[ib,'cells'].sum()),lung_cells=int(l.loc[il,'cells'].sum()),bone_patients=';'.join(ib),lung_patients=';'.join(il),LOPO_min=min(loo) if loo else np.nan,LOPO_max=max(loo) if loo else np.nan,LOPO_direction_retained=bool(all(v*effect>0 for v in loo)) if loo else None,family=f'{scope}|{ds}|{matrix}|{minimum}|{design}|{endpoint}')
     comparisons.append(st)
 stats=adjust(pd.DataFrame(comparisons));save(stats,'all_functional_bone_lung_statistics.tsv');save(stats[stats.nominal],'nominal_functional_bone_lung_statistics.tsv')
 print('FUNCTIONAL_STATS',len(stats),'nominal',int(stats.nominal.sum()),flush=True)
 # Continuous score depth dependence is a diagnostic, not evidence of mechanism.
 for (scope,ds,matrix,pr),g in pat.groupby(['scope','dataset','matrix','program']):
  quality.append(dict(scope=scope,dataset=ds,matrix=matrix,program=pr,patient_site_records=len(g),score_depth_spearman=g.mean_score.corr(g.detected_features_mean,method='spearman'),score_mito_spearman=g.mean_score.corr(g.mitochondrial_fraction_mean,method='spearman')))
 save(pd.DataFrame(quality),'program_technical_covariation.tsv')
 # Prepare open membrane search plus all14, with no source18/other incomplete data used.
 roles=read(paths['roles'],keep_default_na=False);surface=set(roles.loc[truth(roles.plasma_membrane_transmembrane)|truth(roles.plasma_membrane_peripheral),'genesymbol'])
 legacy=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates'];surface.update(legacy)
 all14=json.loads((P/'CAF_ALL14_RESEARCH_CURRENT.json').read_text())['retained_genes'];symbols=sorted(surface|set(all14));dump(dict(surface_symbols=sorted(surface),retained14=all14),O/'resources/MARKER_SEARCH_UNIVERSE.json')
 associations=[];marker_patients=[];detection=[];skip=[]
 for scope,mask in scopes.items():
  for ds in DS:
   q=ob[mask&ob.dataset.eq(ds)];genes=[g for g in symbols if g in maps[ds]];js=[maps[ds][g] for g in genes];rc=raw[:,js].toarray()
   for (donor,site),gg in q.groupby(['donor_record','strict_site']):
    ix=ob.index.get_indexer(gg.index);v=rc[ix];ss=gg['sample'].unique();assayed=np.asarray(presence[[sidx[s] for s in ss]][:,js].sum(0)).ravel()==len(ss);ex=np.log2(v.sum(0)/library[ix].sum()*1e4+.1);det=(v>0).mean(0);ex[~assayed]=np.nan;det[~assayed]=np.nan
    for j,gene in enumerate(genes):detection.append(dict(scope=scope,dataset=ds,donor_record=donor,site=site,gene=gene,cells=len(ix),RNA_detected_fraction=det[j],log2_pseudobulk_CP10k_plus01=ex[j],surface_annotation=gene in surface,retained14=gene in all14))
   for pr in panels:
    donor_effects={};donor_deteffects={}
    for (donor,site),gg in q.groupby(['donor_record','strict_site']):
     ix=ob.index.get_indexer(gg.index);v=score['raw',pr][ix];valid=np.isfinite(v)
     lo=float(np.nanquantile(v,.25)) if valid.any() else np.nan;hi=float(np.nanquantile(v,.75)) if valid.any() else np.nan
     high=ix[valid&(v>hi)];low=ix[valid&(v<lo)]
     if len(high)<5 or len(low)<5:
      skip.append(dict(scope=scope,dataset=ds,program=pr,donor_record=donor,site=site,cells=len(ix),high_cells=len(high),low_cells=len(low),reason='fewer_than5_per_arm_or_score_ties'));continue
     ss=gg['sample'].unique();assayed=np.asarray(presence[[sidx[s] for s in ss]][:,js].sum(0)).ravel()==len(ss)
     hh=np.log2(rc[high].sum(0)/library[high].sum()*1e4+.1);ll=np.log2(rc[low].sum(0)/library[low].sum()*1e4+.1);d=hh-ll;dd=(rc[high]>0).mean(0)-(rc[low]>0).mean(0)
     circular=np.array([g in panels[pr] for g in genes]);d[~assayed|circular]=np.nan;dd[~assayed|circular]=np.nan
     donor_effects.setdefault(donor,[]).append(d);donor_deteffects.setdefault(donor,[]).append(dd)
    ids=sorted(donor_effects)
    if not ids:continue
    d=np.vstack([np.nanmean(donor_effects[i],axis=0) for i in ids]);dd=np.vstack([np.nanmean(donor_deteffects[i],axis=0) for i in ids]);st=stat(np.zeros_like(d),d,True).rename(columns={'mean_lung_minus_bone':'mean_high_minus_low','median_lung_minus_bone':'median_high_minus_low','bone_N':'N','lung_N':'N_repeated'})
    st['gene']=genes;st['mean_detection_high_minus_low']=np.nanmean(dd,axis=0);st['scope']=scope;st['dataset']=ds;st['program']=pr;st['surface_annotation']=[g in surface for g in genes];st['retained14']=[g in all14 for g in genes];st['component_of_score']=[g in panels[pr] for g in genes];st['family']=f'{scope}|{ds}|{pr}|within_patient_marker_association'
    associations.append(st)
    for i,donor in enumerate(ids):
     for j,gene in enumerate(genes):
      if np.isfinite(d[i,j]):marker_patients.append(dict(scope=scope,dataset=ds,program=pr,donor_record=donor,gene=gene,high_minus_low=d[i,j],detection_high_minus_low=dd[i,j],patient_sites=len(donor_effects[donor])))
   print('SURFACE',scope,ds,len(genes),flush=True)
 assocs=adjust(pd.concat(associations,ignore_index=True));save(assocs,'all_marker_functional_associations.tsv.gz');save(assocs[assocs.nominal],'nominal_marker_functional_associations.tsv');save(assocs[assocs.retained14],'ALL14_functional_associations.tsv');save(pd.DataFrame(marker_patients),'marker_patient_effects.tsv.gz');save(pd.DataFrame(detection),'marker_patient_detection.tsv.gz');save(pd.DataFrame(skip),'marker_association_insufficient_cells.tsv')
 # Explicitly retain a row for every accepted candidate, even if no evaluable association.
 roster14=[]
 for gene in all14:
  g=assocs[assocs.gene.eq(gene)];roster14.append(dict(gene=gene,retain=True,status='retained',surface_annotation=gene in surface,evaluable_associations=int(g.P.notna().sum()),nominal_associations=int(g.nominal.sum()),not_evaluable_reason='' if g.P.notna().any() else 'insufficient_repeated_donors_or_no_variation_or_unmapped'))
 save(pd.DataFrame(roster14),'ALL14_retention_registry.tsv')
 scope_counts=pd.concat([ob[m].groupby(['dataset','strict_site']).agg(cells=('cell_id','size'),patients=('donor_record','nunique')).reset_index().assign(scope=k) for k,m in scopes.items()]);save(scope_counts,'scope_patient_cell_counts.tsv')
 unchanged=all(sha(paths[k])==v['sha256'] for k,v in manifest.items());assert unchanged
 done=dict(status='computed_pending_independent_QA',completed_at=datetime.now(timezone.utc).isoformat(),programs=len(panels),functional_tests=len(stats),marker_association_rows=len(assocs),all14_retained=True,original_inputs_unchanged=unchanged,new_patient_level_tests=True,max_score_reconstruction_error=max(x['max_error'] for x in recon_errors))
 dump(done,O/'COMPUTATION_COMPLETE.json');print(json.dumps(done),flush=True)
 mainst=stats[(stats.scope=='Fibro_reviewed')&(stats.matrix=='raw')&(stats.minimum_cells==1)&(stats.design=='independent')&(stats.endpoint=='mean_score')]
 print('MAIN_RESULTS\n'+mainst[['dataset','program','bone_N','lung_N','mean_lung_minus_bone','P','FDR']].to_string(index=False),flush=True)

if __name__=='__main__':
 try:main()
 except BaseException:
  (O/'logs/analysis.exit_code').write_text('1\n');traceback.print_exc();raise
 else:(O/'logs/analysis.exit_code').write_text('0\n')
