from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='4'
import json, hashlib, datetime, itertools, sys, gc
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'); B=P/'phase58_F1F5_round1_20261004'; A=P/'phase41_full_transcriptome_setup_20261001'; O=P/'phase59_four_functional_axes_20261004'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def validmean(x,axis=0):
 n=np.isfinite(x).sum(axis); out=np.full(np.shape(n),np.nan);np.divide(np.nansum(x,axis=axis),n,out=out,where=n>0);return out
def truth(s):return s.astype(str).str.lower().isin(['true','1'])
for d in ['code','resources','tables','objects','logs','figures','delivery','QA']:(O/d).mkdir(parents=True,exist_ok=True)
libs=json.loads((B/'resources/libraries_canonical.json').read_text()); definitions={}; roles={}; origins={}
primary={'inflammation':('hallmark','HALLMARK_INFLAMMATORY_RESPONSE'),'contraction':('reactome','REACTOME_SMOOTH_MUSCLE_CONTRACTION'),'matrix':('matrisome','NABA_CORE_MATRISOME'),'osteogenic':('go_bp','GOBP_OSTEOBLAST_DIFFERENTIATION')}
for axis,(lib,key) in primary.items():definitions[axis]=sorted(set(libs[lib][key]));roles[axis]='primary_four_axes';origins[axis]=lib+'|'+key
secondary={'IL6_JAK_STAT3':('hallmark','HALLMARK_IL6_JAK_STAT3_SIGNALING'),'TGF_beta':('hallmark','HALLMARK_TGF_BETA_SIGNALING'),'collagen_formation':('reactome','REACTOME_COLLAGEN_FORMATION'),'ECM_degradation':('reactome','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX'),'ECM_regulators':('matrisome','NABA_ECM_REGULATORS'),'ossification':('go_bp','GOBP_OSSIFICATION'),'mineralization':('go_bp','GOBP_BONE_MINERALIZATION'),'RUNX2_osteoblast':('reactome','REACTOME_RUNX2_REGULATES_OSTEOBLAST_DIFFERENTIATION')}
for name,(lib,key) in secondary.items():definitions[name]=sorted(set(libs[lib][key]));roles[name]='secondary_pathways';origins[name]=lib+'|'+key
panels={
 'inflammatory_secretory':['IL6','IL11','LIF','CXCL1','CXCL12'],
 'OS_iCAF_reference':['APOD','IGF1','ITM2A'],
 'contractile_markers':['ACTA2','TAGLN','MYL9','MYH11','DES'],
 'matrix_markers':['MMP11','POSTN','COL1A1','COL1A2','COL3A1','COL5A2','COL8A1','COL10A1','COL11A1','COL12A1','COMP','CTHRC1','INHBA','SULF1','VCAN','ASPN'],
 'osteogenic_effectors':['RUNX2','SP7','ALPL','IBSP','BGLAP'],
 'Liu_C2_OMD_OGN':['OMD','OGN'],
 'Liu_C5_ossification':['SMPD3','SOSTDC1','PHOSPHO1'],
 'LOX_crosslink':['LOX','LOXL1','LOXL2','LOXL3','LOXL4'],
 'mural_competitor':['RGS5','CSPG4','MCAM','NOTCH3','ABCC9','KCNJ8','MYH11'],
 'fibroblast_reference':['DCN','LUM','PDGFRA','COL14A1','DPT'],
 'LRRC15_activation':['LRRC15'],
}
refs={'inflammatory_secretory':['10.1084/jem.20162024','10.1158/2159-8290.CD-18-0710'],'OS_iCAF_reference':['10.3389/fimmu.2025.1680554'],'contractile_markers':['10.1038/s41467-020-20059-6','10.1158/2159-8290.CD-19-0094'],'matrix_markers':['10.1038/s41467-023-39762-1'],'osteogenic_effectors':['10.3389/fonc.2021.709210','10.1038/s41467-020-20059-6'],'Liu_C2_OMD_OGN':['10.3389/fonc.2021.709210'],'Liu_C5_ossification':['10.3389/fonc.2021.709210'],'LOX_crosslink':['10.7150/thno.73714'],'mural_competitor':['existing_project_lineage_panels_and_competitor_annotation'],'fibroblast_reference':['10.1038/s41467-020-20059-6','10.1158/2159-8290.CD-19-0094'],'LRRC15_activation':['10.1038/s41586-022-05272-1']}
for name,genes in panels.items():definitions[name]=genes;roles[name]='secondary_reference_panels';origins[name]='analyst_prespecified_reference_panel_not_validated_classifier'
for axis in primary:
 others=set().union(*(set(definitions[k]) for k in primary if k!=axis));name=axis+'_nonoverlap';definitions[name]=sorted(set(definitions[axis])-others);roles[name]='primary_nonoverlap_sensitivity';origins[name]='primary_axis_minus_other_three_axis_genes'
programs=list(definitions);genes=sorted(set().union(*map(set,definitions.values())));geneix={g:i for i,g in enumerate(genes)};members={p:np.array([geneix[g] for g in gs],int) for p,gs in definitions.items()}
protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Four overlapping functional axes in frozen original F1-F5; focus F2 contraction/matrix and F4 osteogenic source; no new clustering, no new CNV, no spatial eggshell or treatment claim','primary_definitions':primary,'definition_policy':'Primary full pathway sets copied verbatim from versioned phase58 canonical library; panels are secondary prespecified references, never a classifier; pre-existing phase58 findings known','normalization':'pseudobulk full-library CP10k; log2(CP10k+0.1); target minus matched reference, site-equal then patient-equal','N':'independent source-local donor_record, no cross-source pooling; descriptive only N<2; N<10 exploratory','program_coverage':'all assigned genes must be assayed across samples in pair; use common assayed genes; at least 80% of definition and at least min(3,definition_size), otherwise not evaluated; expression zero distinct from unmeasured','test':'two-sided patient signflip exact N<=18; otherwise MC99999 + observed, seed20261004; all-zero effects not tested','families':'dataset|matrix|analysis_context|primary_four_axes / secondary_pathways / secondary_reference_panels / primary_nonoverlap_sensitivity; all five class contrasts together for main; genes separately','cell_threshold':'primary at least 1 cell per paired arm; explicit sensitivity at least 10 cells per patient-site arm, not a new identity definition','identity_sensitivity':'review Fibro only; conservative review Fibro with mixed/doublet flags removed and high/medium lineage confidence (unassessed doublet not assumed negative); F2/F4 matched within-class lineages; no CNV-negativity identity claim','external_background':'raw cells outside original8174, same dataset/patient/site, expression review strata Fibro/Mural/Osteo; corrected sensitivity confined original8174','coactivity':'patient effects on all four axes retained as continuous vectors; axes not mutually exclusive; positive effect not absolute cell identity','gene_source_attribution':'within original class, raw/corrected selected-marker RNA contribution and cells by reviewed lineage, both pooled and patient-equal; expression lineage not proven origin or malignancy','known_ascertainment':'original8174 selected under historical CAF/CNV rules; selected CNV distribution cannot establish absence of malignant cells; review labels use RNA and are not independent functional validation'}
if (O/'resources/PROTOCOL.json').exists():
 old=json.loads((O/'resources/PROTOCOL.json').read_text()); assert old['primary_definitions']==json.loads(json.dumps(primary))
else:
 dump(O/'resources/PROTOCOL.json',protocol);dump(O/'resources/PREVIOUS_CURRENT.json',json.loads((P/'CAF_NEXT_STAGE_CURRENT.json').read_text()))
dump(O/'resources/GENE_SETS.json',{'definitions':definitions,'roles':roles,'origins':origins,'panel_references':refs});pd.DataFrame([{'program':p,'role':roles[p],'source':origins[p],'gene':g} for p,v in definitions.items() for g in v]).to_csv(O/'resources/gene_set_membership.tsv',sep='\t',index=False)
pd.DataFrame([{'axis1':i,'axis2':j,'shared':len(set(definitions[i])&set(definitions[j])),'genes':';'.join(sorted(set(definitions[i])&set(definitions[j])))} for i in primary for j in primary]).to_csv(O/'resources/primary_axis_overlap.tsv',sep='\t',index=False)
print('FROZEN',len(programs),'programs',len(genes),'genes',flush=True)
rawpath=A/'objects/all_candidate_full_genes.h5ad';corpath=B/'objects/expanded_corrected_full_genes.h5ad';annopath=P/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz';ob=read(B/'tables/frozen_expanded_cells.tsv.gz').set_index('cell_id');anno=read(annopath)
rv=anno.set_index('original_cell_id').loc[ob.index];assert len(ob)==8174 and rv.index.is_unique
for col in ['review_lineage_code','CNV_evidence_current','lineage_confidence','mixed_program_flag','doublet_assessed_any','doublet_positive_any','decontx_contamination','candidate_tier']:ob[col]=rv[col].to_numpy()
ob['lineage_stratum']=np.select([ob.review_lineage_code.eq('Fibro'),ob.review_lineage_code.isin(['Osteo','Osteochondro','Chondro']),ob.review_lineage_code.isin(['Mural','SMC'])],['Fibro','Osteo','Mural'],default='Other_unresolved')
ob['conservative_fibro']=ob.lineage_stratum.eq('Fibro')&~truth(ob.mixed_program_flag)&~truth(ob.doublet_positive_any)&ob.lineage_confidence.isin(['高','中'])
ob.to_csv(O/'tables/frozen_cells_with_identity.tsv.gz',sep='\t')
ob.groupby(['class_id','lineage_stratum','CNV_evidence_current'],observed=True).size().rename('cells').reset_index().to_csv(O/'tables/class_identity_CNV_counts.tsv',sep='\t',index=False)
mp=read(B/'resources/source_canonical_feature_map.tsv');sources=sorted(ob.dataset.unique());mapping={ds:mp[mp.dataset.eq(ds)].set_index('symbol').feature.to_dict() for ds in sources}
presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=read(A/'tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict()
a=ad.read_h5ad(rawpath);assert a.obs_names.tolist()==ob.index.tolist();features=a.var_names.astype(str);fidx={g:i for i,g in enumerate(features)}
assay=np.zeros((len(ob),len(genes)),bool);counts={};denoms={};coverage=[]
for mode,path in [('raw',rawpath),('corrected',corpath)]:
 x=a if mode=='raw' else ad.read_h5ad(path);assert x.obs_names.tolist()==ob.index.tolist();den=np.asarray(x.X.sum(1)).ravel();expected=ob.full_library_counts.to_numpy(float) if mode=='raw' else x.obs.corrected_full_library_counts.to_numpy(float);assert np.allclose(den,expected,rtol=1e-6,atol=1e-4);v=np.zeros((len(ob),len(genes)),float)
 for ds in sources:
  rows=np.flatnonzero(ob.dataset.eq(ds));ix=np.array([fidx.get(mapping[ds].get(g,''),-1) for g in genes]);ok=ix>=0;v[np.ix_(rows,np.flatnonzero(ok))]=x.X[rows][:,ix[ok]].toarray()
  if mode=='raw':
   for sample in ob.iloc[rows]['sample'].unique():
    ss=np.flatnonzero(ob['sample'].eq(sample));assay[np.ix_(ss,np.flatnonzero(ok))]=presence[sidx[sample],ix[ok]].toarray().ravel()>0
   for p,ii in members.items():coverage.append({'dataset':ds,'program':p,'definition_genes':len(ii),'assayed_genes_any':int(assay[rows][:,ii].any(0).sum()),'assayed_genes_all':int(assay[rows][:,ii].all(0).sum()),'missing_all_source':';'.join(genes[i] for i in ii if not assay[rows,i].any())})
 counts[mode]=v;denoms[mode]=den
 if mode=='corrected':del x
del a;gc.collect();pd.DataFrame(coverage).to_csv(O/'tables/program_coverage.tsv',sep='\t',index=False)
dump(O/'resources/INPUT_MANIFEST.json',[{'path':str(p),'bytes':p.stat().st_size,'SHA256':sha(p)} for p in [rawpath,corpath,annopath,B/'resources/libraries_canonical.json',B/'resources/source_canonical_feature_map.tsv']])
np.savez_compressed(O/'objects/selected_gene_counts.npz',raw=counts['raw'],corrected=counts['corrected'],raw_library=denoms['raw'],corrected_library=denoms['corrected'],assayed=assay,genes=np.array(genes),cells=ob.index.to_numpy(str))
signcache={}
def test(v):
 v=np.asarray(v);v=v[np.isfinite(v)];n=len(v)
 if n<2:return np.nan,'not_tested_N_lt2'
 if np.all(v==0):return np.nan,'not_tested_all_zero'
 if n not in signcache:
  if n<=18:signcache[n]=2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(float)-1
  else:signcache[n]=np.random.default_rng(20261004).choice([-1.,1.],size=(99999,n))
 null=np.abs(signcache[n]@v/n);hits=int((null>=abs(v.mean())-1e-12).sum());return (hits/len(null),'exact_patient_signflip') if n<=18 else ((hits+1)/(len(null)+1),'MC99999_plus_observed')
effects=[];geneeffects=[];roster=[];stats=[];genestats=[];trace=[];keymarkers=sorted(set().union(*map(set,panels.values())));keyix=[geneix[g] for g in keymarkers]
def compare(ds,mode,context,contrast,left,right,mincells=1,external=None):
 # left/right are index arrays of original cells, or external grouped pseudobulks.
 donor_rows={};gene_rows={};lframe=ob.iloc[left];rframe=ob.iloc[right] if external is None else None
 for (donor,site),lg in lframe.groupby(['donor_record','site'],observed=True):
  il=ob.index.get_indexer(lg.index)
  if len(il)<mincells:continue
  if external is None:
   rg=rframe[(rframe.donor_record==donor)&(rframe.site==site)];ir=ob.index.get_indexer(rg.index)
   if len(ir)<mincells:continue
   rc=counts[mode][ir].sum(0);rd=denoms[mode][ir].sum();ra=assay[ir].all(0);nr=len(ir)
  else:
   z=external.get((ds,str(donor),str(site)))
   if z is None or z['cells']<mincells:continue
   rc=z['counts'];rd=z['library'];ra=z['assayed'];nr=z['cells']
  ld=denoms[mode][il].sum()
  if ld<=0 or rd<=0:continue
  la=assay[il].all(0);valid=la&ra;dd=np.log2(counts[mode][il].sum(0)/ld*1e4+.1)-np.log2(rc/rd*1e4+.1);dd[~valid]=np.nan
  score=[]
  for p,ii in members.items():
   n=int(valid[ii].sum());score.append(float(np.nanmean(dd[ii])) if n>=max(min(3,len(ii)),int(np.ceil(.8*len(ii)))) else np.nan)
  donor_rows.setdefault(str(donor),[]).append(score);gene_rows.setdefault(str(donor),[]).append(dd[keyix]);roster.append({'dataset':ds,'matrix':mode,'context':context,'contrast':contrast,'donor_record':donor,'site':site,'target_cells':len(il),'reference_cells':nr,'minimum_cells':mincells})
  if ds=='SCPCP000017' and context=='main':trace.append({'matrix':mode,'contrast':contrast,'donor_record':donor,'site':site,'left_indices':il.tolist(),'right_indices':ir.tolist()})
 for j,p in enumerate(programs):
  vv=[]
  for donor,rows in donor_rows.items():
   value=float(validmean(np.array(rows)[:,j]));vv.append(value);effects.append({'dataset':ds,'matrix':mode,'context':context,'contrast':contrast,'program':p,'donor_record':donor,'effect':value,'paired_sites':len(rows)})
  v=np.array(vv);v=v[np.isfinite(v)];pv,method=test(v);n=len(v);stats.append({'dataset':ds,'matrix':mode,'context':context,'contrast':contrast,'program':p,'role':roles[p],'N':n,'mean_delta':float(v.mean()) if n else np.nan,'median_delta':float(np.median(v)) if n else np.nan,'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'P':pv,'test':method,'family':'|'.join([ds,mode,context,roles[p]]),'LOO_min':float(min((v.sum()-v)/(n-1))) if n>1 else np.nan,'LOO_max':float(max((v.sum()-v)/(n-1))) if n>1 else np.nan,'evidence':'exploratory_N_lt10' if n<10 else 'patient_level'})
 for j,gene in enumerate(keymarkers):
  vv=[]
  for donor,rows in gene_rows.items():
   val=float(validmean(np.array(rows)[:,j]));vv.append(val);geneeffects.append({'dataset':ds,'matrix':mode,'context':context,'contrast':contrast,'gene':gene,'donor_record':donor,'effect':val})
  v=np.array(vv);v=v[np.isfinite(v)];pv,method=test(v);genestats.append({'dataset':ds,'matrix':mode,'context':context,'contrast':contrast,'gene':gene,'N':len(v),'mean_delta':v.mean() if len(v) else np.nan,'positive_patients':int((v>0).sum()),'negative_patients':int((v<0).sum()),'P':pv,'test':method,'family':'|'.join([ds,mode,context,'marker_genes'])})
def ix(mask):return np.flatnonzero(np.asarray(mask))
for mode in counts:
 for ds in sources:
  dsok=ob.dataset.eq(ds)
  for cl in ['F1','F2','F3','F4','F5']:
   target=ob.class_id.eq(cl)
   for context,allowed,minimum in [('main',np.ones(len(ob),bool),1),('min10',np.ones(len(ob),bool),10),('Fibro_only',ob.lineage_stratum.eq('Fibro'),1),('conservative_Fibro',ob.conservative_fibro,1)]:compare(ds,mode,context,cl+'_vs_rest',ix(dsok&allowed&target),ix(dsok&allowed&~target),minimum)
  for cl in ['F2','F4']:
   clok=ob.class_id.eq(cl)&dsok
   for ls,rs in [('Fibro','Osteo'),('Fibro','Mural'),('Osteo','Other_unresolved')]:compare(ds,mode,'within_'+cl,ls+'_vs_'+rs,ix(clok&ob.lineage_stratum.eq(ls)),ix(clok&ob.lineage_stratum.eq(rs)))
  print('CANDIDATE_TESTS',mode,ds,flush=True)
def adjust(df):
 df=df.copy();df['FDR']=np.nan
 for fam,g in df.groupby('family',sort=False):
  # All prespecified rows retained in family; unevaluable rows use 1 only in BH denominator.
  q=multipletests(g.P.fillna(1).to_numpy(),method='fdr_bh')[1];df.loc[g.index,'FDR']=np.where(g.P.notna(),q,np.nan)
 df['nominal']=df.P.lt(.05);df['FDR_significant']=df.FDR.lt(.05);return df
def save_stats():
 adjust(pd.DataFrame(stats)).to_csv(O/'tables/axis_statistics.tsv',sep='\t',index=False);adjust(pd.DataFrame(genestats)).to_csv(O/'tables/marker_statistics.tsv.gz',sep='\t',index=False);pd.DataFrame(effects).to_csv(O/'tables/patient_axis_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(geneeffects).to_csv(O/'tables/patient_marker_effects.tsv.gz',sep='\t',index=False);pd.DataFrame(roster).to_csv(O/'tables/paired_patient_site_roster.tsv.gz',sep='\t',index=False)
save_stats();dump(O/'QA/main_pseudobulk_reconstruction_indices.json',trace)
# RNA contributions are descriptive. Sum only assayed genes; denominator within same patient-site/class/gene.
contrib=[];composition=[]
for (ds,cl,donor,site),g in ob.groupby(['dataset','class_id','donor_record','site'],observed=True):
 total=ob.index.get_indexer(g.index)
 for lineage in ['Fibro','Osteo','Mural','Other_unresolved']:
  h=g[g.lineage_stratum.eq(lineage)]
  ii=ob.index.get_indexer(h.index);composition.append({'dataset':ds,'class_id':cl,'donor_record':donor,'site':site,'lineage':lineage,'cells':len(ii),'class_cells':len(total),'cell_fraction':len(ii)/len(total)})
  if cl not in ['F2','F4']:continue
  for mode in counts:
   for gene in keymarkers:
    j=geneix[gene];ass=bool(assay[total,j].all());den=counts[mode][total,j].sum() if ass else np.nan;num=counts[mode][ii,j].sum() if ass else np.nan;contrib.append({'dataset':ds,'class_id':cl,'donor_record':donor,'site':site,'lineage':lineage,'matrix':mode,'gene':gene,'cells':len(ii),'class_cells':len(total),'RNA_counts':num,'total_class_RNA_counts':den,'RNA_fraction':num/den if den>0 else np.nan,'detected_cells':int((counts[mode][ii,j]>0).sum()) if ass else np.nan,'assayed':ass})
pd.DataFrame(contrib).to_csv(O/'tables/F2F4_lineage_RNA_contribution.tsv.gz',sep='\t',index=False);pd.DataFrame(composition).to_csv(O/'tables/patient_class_lineage_composition.tsv',sep='\t',index=False)
# Coexpression for description, not statistical replication: corrected abundance >0 and >=2 distinct panel genes.
co=[]
for mode in counts:
 for axis,panel in [('contraction','contractile_markers'),('matrix','matrix_markers'),('osteogenic','osteogenic_effectors'),('inflammation','inflammatory_secretory')]:
  ii=members[panel];co.append(pd.DataFrame({'cell_id':ob.index,'matrix':mode,'axis':axis,'detected_panel_genes':(counts[mode][:,ii]>=1).sum(1),'assayed_panel_genes':assay[:,ii].sum(1),'definition_genes':len(ii),'count_detection_threshold':1}))
pd.concat(co).to_csv(O/'tables/cell_axis_codetection.tsv.gz',sep='\t',index=False)
print('CANDIDATE_COMPLETE',flush=True);(O/'CANDIDATE_COMPLETE').write_text('SUCCESS\n')
# Build same-patient/site competitor references from outside frozen original cells, for F2/F4 source attribution.
backgrounds={g:{} for g in ['Fibro','Mural','Osteo']};ids=set(ob.index);focus_donors=set(ob.loc[ob.class_id.isin(['F2','F4']),'donor_record']);bgann=anno[anno.original_donor_record.isin(focus_donors)]
bgann=bgann.copy();bgann['stratum']=np.select([bgann.review_lineage_code.eq('Fibro'),bgann.review_lineage_code.isin(['Mural','SMC']),bgann.review_lineage_code.isin(['Osteo','Chondro','Osteochondro'])],['Fibro','Mural','Osteo'],default='Other');bgann=bgann[~bgann.original_cell_id.isin(ids)&bgann.stratum.ne('Other')]
bgledger=[]
for sn,(sample_id,sub) in enumerate(bgann.groupby('sample_id',sort=True),1):
 path=P/'phase53_de_novo_reannotation_20261002/blind_review/raw_counts'/f'{sample_id}.h5ad';x=ad.read_h5ad(path,backed='r');ds=sub.dataset.iloc[0];ixg=x.var_names.get_indexer([mapping[ds].get(g,'__absent__') for g in genes]);ok=ixg>=0
 for (donor,site,lineage),part in sub.groupby(['original_donor_record','site','stratum'],observed=True):
  rows=x.obs_names.get_indexer(part.cell_id);assert (rows>=0).all();v=np.zeros(len(genes)); mat=x.X[rows];v[ok]=np.asarray(mat[:,ixg[ok]].sum(0)).ravel();den=float(x.obs.iloc[rows].library_counts.sum());key=(ds,str(donor),str(site));z=backgrounds[lineage].get(key)
  if z is None:backgrounds[lineage][key]={'counts':v,'library':den,'cells':len(rows),'assayed':ok.copy()}
  else:z['counts']+=v;z['library']+=den;z['cells']+=len(rows);z['assayed']&=ok
  bgledger.append({'dataset':ds,'sample_id':sample_id,'donor_record':donor,'site':site,'lineage':lineage,'cells':len(rows),'library':den,'assayed_genes':int(ok.sum())})
 x.file.close();del x;gc.collect();print('BACKGROUND',sn,sample_id,flush=True)
np.savez_compressed(O/'objects/background_pseudobulks.npz',**{lineage:np.array([{'key':key,**v} for key,v in records.items()],dtype=object) for lineage,records in backgrounds.items()})
pd.DataFrame(bgledger).to_csv(O/'tables/background_sample_roster.tsv',sep='\t',index=False)
for ds in sources:
 for cl in ['F2','F4']:
  for stratum in ['All','Fibro','Osteo','Mural']:
   allowed=ob.dataset.eq(ds)&ob.class_id.eq(cl)
   if stratum!='All':allowed&=ob.lineage_stratum.eq(stratum)
   for bg in backgrounds:
    compare(ds,'raw','external_'+cl,cl+'_'+stratum+'_vs_BG_'+bg,ix(allowed),[],1,backgrounds[bg])
 print('EXTERNAL_TESTS',ds,flush=True)
save_stats();dump(O/'delivery/ANALYSIS_ACCEPTANCE.json',{'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cells':len(ob),'classes':ob.class_id.value_counts().to_dict(),'programs':len(programs),'primary_axes':list(primary),'source_count':len(sources),'program_tests':len(stats),'marker_tests':len(genestats),'new_CNV':False,'new_clustering':False,'external_background_cells':sum(x['cells'] for x in bgledger),'independent_QA':'pending'});(O/'ANALYSIS_COMPLETE').write_text('SUCCESS\n');print('ANALYSIS_COMPLETE',flush=True)
