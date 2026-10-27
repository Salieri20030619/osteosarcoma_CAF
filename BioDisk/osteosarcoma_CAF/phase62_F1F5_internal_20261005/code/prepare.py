from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,gc,sys,itertools,shutil
import pandas as pd,numpy as np,anndata as ad
from scipy import sparse
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'); B=P/'phase58_F1F5_round1_20261004'; A=P/'phase41_full_transcriptome_setup_20261001'; O=P/'phase62_F1F5_internal_20261005'; Q=P/'phase59_four_functional_axes_20261004'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(8388608),b''):h.update(z)
 return h.hexdigest()
for d in ['resources','objects','tables/genes','tables/gsea','tables/paired','figures','QA','delivery']:(O/d).mkdir(parents=True,exist_ok=True)
ob=read(Q/'tables/frozen_cells_with_identity.tsv.gz').set_index('cell_id');qc=read(B/'tables/cell_QC_and_class.tsv.gz').set_index('cell_id');ob=ob.join(qc[['UMI','detected_features','mitochondrial_fraction','ribosomal_fraction']]);assert len(ob)==8174 and ob.index.is_unique
ob['subgroup_id']=ob.round1_state.str.replace('expanded_','',regex=False);ob['F4_major_branch']='not_F4';f4=read(P/'phase60_F3_review_20261005/tables/F4_two_main_strata_with_residual.tsv.gz').set_index('cell_id');ob.loc[f4.index,'F4_major_branch']=f4.analysis_group
states=ob.groupby(['class_id','subgroup_id']).agg(cells=('sample','size'),donor_records=('donor_record','nunique'),sources=('dataset','nunique'),max_donor_fraction=('donor_record',lambda v:v.value_counts(normalize=True).max()),median_UMI=('UMI','median'),median_MT=('mitochondrial_fraction','median')).reset_index()
primary={s:g.dataset.value_counts().sort_index().sort_values(ascending=False,kind='stable').index[0] for s,g in ob.groupby('subgroup_id')};states['primary_source']=states.subgroup_id.map(primary);states.to_csv(O/'tables/subgroup_roster.tsv',sep='\t',index=False)
ob.groupby(['class_id','subgroup_id','dataset','donor_record','site','lineage_stratum','lineage_confidence'],dropna=False).size().rename('cells').reset_index().to_csv(O/'tables/patient_lineage_distribution.tsv',sep='\t',index=False)
protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'All F1-F5 internal working subgroups, retained confidence. Reuse 33 frozen clusters, F4 two main expression branches plus residual; no new clustering or biological subtype claim.','primary_source':'Per subgroup source with most target cells, selected by counts before expression; all other source contrasts retained where both arms exist. Never pool cross-source gene measurements.','full_genes':'Every canonical measured feature per source, restricted to measured in all source cells used; absent features excluded, never zero. Full-library denominators.','gene_metrics':'Descriptive pseudobulk log2(CP10k+.1) target minus other subgroups of SAME F; detection raw>=1, corrected>=.5; no cell-as-N P. Full tables retained.','sensitivity':'raw/corrected; corrected Fibro restricted when >=5 cells per arm; per-donor/site matched corrected effects with >=5 cells per arm, site-equal and donor-equal. Within-source only.','biological_statistics':'Paired donor signflip exact N<=15, MC19999+observed otherwise; N<2/allzero no P; N<10 exploratory. Gene BH family = source x parentF x all internal contrasts, separately from all pathways. Raw P<.05 nominal screening, FDR reported not required.','GSEA':'Primary source corrected full transcriptome ranks for all states plus evaluable Fibro sensitivity; Hallmark/Reactome/Matrisome; 1000 gene-set permutations seed20261005; model P and native q plus BH by context/library, not biological sample inference; no P-prefilter. F3 phase61 original portraits also preserved as prior evidence.','display':'Top12 descriptive markers ranked positive delta*detection_difference (target detection>=.1); not statistical nomination. Stable display markers delta>=.5, detection difference>=.1, raw same direction.','confidence':'Heuristic ordinal EVIDENCE grades, not calibrated assignment probabilities. Expression high >=3 stable markers and >=80% raw sign preservation among top12, medium >=1, otherwise low. Cross-patient support high >=10 paired source-local donors with top-marker positive direction >=70%; medium >=3 and >=60%; otherwise limited. Each source graded separately; strongest source recorded with its name. Overall = high if expression high and replication high, medium if expression >=medium and replication >=medium, else exploratory. No cluster removed because low confidence.','retained_confidence':'Original cell lineage_confidence, candidate_tier, CNV evidence, mixed/doublet flags unchanged. Appended subgroup-level confidence is NOT cell-level certainty.','external_F_comparisons':'Not needed for defining internal states; not computed.'}
dump(O/'resources/PROTOCOL.json',protocol);dump(O/'resources/PREVIOUS_CURRENT.json',json.loads((P/'CAF_NEXT_STAGE_CURRENT.json').read_text()));libs=json.loads((B/'resources/libraries_canonical.json').read_text());dump(O/'resources/libraries_canonical.json',{k:libs[k] for k in ['hallmark','reactome','matrisome']})
paths=[Q/'tables/frozen_cells_with_identity.tsv.gz',A/'objects/all_candidate_full_genes.h5ad',A/'objects/sample_feature_presence.npz',B/'objects/expanded_corrected_full_genes.h5ad',B/'resources/source_canonical_feature_map.tsv',B/'resources/libraries_canonical.json'];dump(O/'resources/INPUT_MANIFEST.json',[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in paths])
presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict();mp=read(B/'resources/source_canonical_feature_map.tsv');mats={};den={}
for mode,path in [('raw',paths[1]),('corrected',paths[3])]:
 a=ad.read_h5ad(path);assert a.obs_names.tolist()==ob.index.tolist();mats[mode]=a.X.tocsr();den[mode]=np.asarray(a.X.sum(1)).ravel();features=a.var_names.copy();del a
ob['corrected_library']=den['corrected'];ob['removed_fraction']=1-den['corrected']/den['raw'];ob.to_csv(O/'tables/cells_original_with_review.tsv.gz',sep='\t')
jobs=[];contrasts=[];top=[];paired_registry=[];coverage=[];statframes=[];marker_patient=[];signcache={}
def test(d):
 n=len(d);out=np.full(d.shape[1],np.nan)
 if n<2:return out
 active=np.any(d!=0,axis=0);ii=np.flatnonzero(active)
 if n not in signcache:signcache[n]=np.array(list(itertools.product([-1.,1.],repeat=n))) if n<=15 else np.random.default_rng(20261005).choice([-1.,1.],size=(19999,n))
 signs=signcache[n]
 for start in range(0,len(ii),128):
  ix=ii[start:start+128];hits=(np.abs(signs@d[:,ix]/n)>=np.abs(d[:,ix].mean(0))[None,:]-1e-12).sum(0);out[ix]=hits/len(signs) if n<=15 else (hits+1)/(len(signs)+1)
 return out
for ds in sorted(ob.dataset.unique()):
 rr=np.flatnonzero(ob.dataset.eq(ds));meta=ob.iloc[rr];mapping=mp[mp.dataset.eq(ds)].sort_values('symbol');jj=features.get_indexer(mapping.feature);assert (jj>=0).all();common=np.asarray(presence[[sidx[s] for s in meta['sample'].unique()]][:,jj].sum(0)).ravel()==meta['sample'].nunique();mapping=mapping.iloc[np.flatnonzero(common)];jj=jj[common];symbols=mapping.symbol.to_numpy(str);assert len(set(symbols))==len(symbols);data={k:v[rr][:,jj].toarray() for k,v in mats.items()};dd={k:v[rr] for k,v in den.items()};coverage.append({'dataset':ds,'canonical_features':len(common),'measured_all_source_samples':len(jj),'samples':meta['sample'].nunique()})
 for cl in sorted(meta.class_id.unique()):
  parent=meta.class_id.eq(cl).to_numpy()
  for st in sorted(meta.loc[parent,'subgroup_id'].unique()):
   left=meta.subgroup_id.eq(st).to_numpy();right=parent&~left
   if not right.any():continue
   primary_here=primary[st]==ds;cache={}
   for mode,scope in [('raw','all'),('corrected','all')]+([('corrected','Fibro')] if primary_here else []):
    l=left.copy();r=right.copy()
    if scope=='Fibro':l &= meta.lineage_stratum.eq('Fibro').to_numpy();r &= meta.lineage_stratum.eq('Fibro').to_numpy()
    il=np.flatnonzero(l);ir=np.flatnonzero(r)
    if not len(il) or not len(ir) or (scope=='Fibro' and min(len(il),len(ir))<5):continue
    x=data[mode];lib=dd[mode];lc=x[il].sum(0)/lib[il].sum()*1e4;rc=x[ir].sum(0)/lib[ir].sum()*1e4;eff=np.log2(lc+.1)-np.log2(rc+.1);threshold=1 if mode=='raw' else .5;lf=(x[il]>=threshold).mean(0);rf=(x[ir]>=threshold).mean(0);ctx=f'{ds}__{st}__{mode}__{scope}'
    fr=pd.DataFrame({'symbol':symbols,'feature':mapping.feature.to_numpy(),'mean_delta':eff,'target_CP10k':lc,'reference_CP10k':rc,'target_detection':lf,'reference_detection':rf,'detection_difference':lf-rf,'any_expression':(lc+rc)>0});fr.to_csv(O/'tables/genes'/f'{ctx}.tsv.gz',sep='\t',index=False);cache[(mode,scope)]=fr
    contrasts.append({'context':ctx,'dataset':ds,'class_id':cl,'subgroup_id':st,'matrix':mode,'scope':scope,'primary':primary_here,'target_cells':len(il),'reference_cells':len(ir),'target_donors':meta.iloc[il].donor_record.nunique(),'reference_donors':meta.iloc[ir].donor_record.nunique(),'features':len(fr)})
    if primary_here:
     t=fr[(fr.target_detection>=.1)&(fr.mean_delta>0)&(fr.detection_difference>0)].copy();t['display_score']=t.mean_delta*t.detection_difference;t=t.sort_values(['display_score','symbol'],ascending=[False,True]).head(12);t['context']=ctx;t['subgroup_id']=st;t['matrix']=mode;t['scope']=scope;top.append(t)
     if mode=='corrected':
      rank=fr[fr.any_expression][['symbol','mean_delta']].sort_values(['mean_delta','symbol'],ascending=[False,True]);rankpath=O/'resources'/f'{ctx}__rank.tsv';rank.to_csv(rankpath,sep='\t',index=False)
      for library in ['hallmark','reactome','matrisome']:jobs.append({'context':ctx,'matrix':mode,'contrast':st,'scope':scope,'library':library,'rank':str(rankpath),'output':str(O/'tables/gsea'/f'{ctx}__{library}.tsv.gz'),'target_cells':len(il),'reference_cells':len(ir),'target_donors':meta.iloc[il].donor_record.nunique(),'reference_donors':meta.iloc[ir].donor_record.nunique()})
   # Matched patient evidence uses full genes and source-local independent records.
   dv=[];donors=[];nsite=[]
   for donor,dg in meta.groupby('donor_record',sort=True):
    siteeffects=[]
    for site,sg in dg.groupby('site',dropna=False):
     ix=meta.index.get_indexer(sg.index);il=ix[left[ix]];ir=ix[right[ix]]
     if min(len(il),len(ir))<5:continue
     x=data['corrected'];lib=dd['corrected'];e=np.log2(x[il].sum(0)/lib[il].sum()*1e4+.1)-np.log2(x[ir].sum(0)/lib[ir].sum()*1e4+.1);siteeffects.append(e);paired_registry.append({'dataset':ds,'class_id':cl,'subgroup_id':st,'donor_record':donor,'site':site,'target_cells':len(il),'reference_cells':len(ir)})
    if siteeffects:dv.append(np.mean(siteeffects,axis=0));donors.append(donor);nsite.append(len(siteeffects))
   d=np.array(dv).reshape((-1,len(symbols)));np.savez_compressed(O/'objects'/f'{ds}__{st}__paired.npz',effects=d,symbols=symbols,donors=np.array(donors,str),sites=np.array(nsite))
   pf=pd.DataFrame({'symbol':symbols,'dataset':ds,'class_id':cl,'subgroup_id':st,'N':len(d),'mean_patient_effect':d.mean(0) if len(d) else np.nan,'positive_patients':(d>0).sum(0),'P':test(d),'family':ds+'|'+cl+'|corrected_internal_genes','test':('not_tested_N_lt2' if len(d)<2 else 'exact_patient_signflip' if len(d)<=15 else 'patient_signflip_MC19999_plus_observed')});statframes.append(pf)
   print('PREPARED',ds,st,'primary',primary_here,'pairedN',len(d),flush=True)
 del data;gc.collect()
stats=pd.concat(statframes,ignore_index=True);stats['FDR']=np.nan
for fam,g in stats.groupby('family'):stats.loc[g.index,'FDR']=np.where(g.P.notna(),multipletests(g.P.fillna(1),method='fdr_bh')[1],np.nan)
stats['nominal']=stats.P.lt(.05);stats['FDR_significant']=stats.FDR.lt(.05);stats.to_csv(O/'tables/paired_gene_statistics.tsv.gz',sep='\t',index=False);stats[stats.nominal].to_csv(O/'tables/paired_nominal_genes.tsv.gz',sep='\t',index=False)
pd.DataFrame(contrasts).to_csv(O/'tables/contrast_registry.tsv',sep='\t',index=False);pd.concat(top).to_csv(O/'tables/marker_display.tsv',sep='\t',index=False);pd.DataFrame(paired_registry).to_csv(O/'tables/paired_donor_site_roster.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/feature_coverage.tsv',sep='\t',index=False);dump(O/'resources/GSEA_JOBS.json',jobs)
# Preserve accepted F3 evidence and explicit labels rather than silently replacing it.
for name in ['F3_internal_working_labels.tsv.gz','contrast_registry.tsv']:
 shutil.copy2(P/'phase61_F3_internal_portraits_20261005/tables'/name,O/'resources'/('phase61_'+name))
(O/'PREPARE_COMPLETE').write_text('SUCCESS\n');print('PREPARE_COMPLETE',len(contrasts),'contrasts',len(jobs),'GSEA jobs',flush=True)
