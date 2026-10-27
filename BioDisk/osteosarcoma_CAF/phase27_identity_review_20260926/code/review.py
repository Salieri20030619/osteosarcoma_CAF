from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
import json,datetime,hashlib,traceback,fcntl,gc
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase26_client_reanalysis_20260922';O=P/'phase27_identity_review_20260926'
for d in ['code','tables','objects','logs','delivery','resources']:(O/d).mkdir(parents=True,exist_ok=True)
def now():return datetime.datetime.now().astimezone().isoformat()
def writej(p,v):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str));t.replace(p)
def status(s,**kw):writej(O/'STATUS.json',dict(status=s,at=now(),**kw));print(s,kw,flush=True)
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(d,n):d.to_csv(O/'tables'/n,sep='\t',index=False)
def clean(a):
 for frame in [a.obs,a.var]:
  for c in frame:
   if frame[c].dtype==object or isinstance(frame[c].dtype,pd.StringDtype):frame[c]=frame[c].fillna('').astype(str).astype(object)
  frame.index=pd.Index(frame.index.astype(str).to_numpy(dtype=object))
def digest(f):return hashlib.sha256(f.read_bytes()).hexdigest()
LINEAGE={
'Fibroblast':['DCN','LUM','DPT','PDGFRA','COL1A1','COL1A2','COL3A1','COL14A1','PI16','FBLN1','CFD','COL15A1'],
'Mural':['RGS5','CSPG4','NOTCH3','KCNJ8','ABCC9','MCAM','RBP1'],
'Osteogenic':['RUNX2','SP7','ALPL','IBSP','BGLAP','PHEX','DMP1','SOST'],
'T_NK':['CD3D','CD3E','TRAC','CD247','NKG7','GNLY'],
'B_Plasma':['CD79A','CD79B','MS4A1','CD37','MZB1','JCHAIN'],
'Myeloid':['PTPRC','LST1','TYROBP','FCER1G','AIF1','CTSS','CD68'],
'Endothelial':['PECAM1','VWF','CDH5','KDR','EMCN','RAMP2'],
'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT14'],
'Mast':['TPSAB1','TPSB2','CPA3','MS4A2','KIT'],
'Myogenic':['MYOD1','MYOG','ACTA1','DES','TNNC2','MYH3'],
'Osteoclast':['ACP5','CTSK','ATP6V0D2','OCSTAMP','DCSTAMP'],
'Chondrocyte':['ACAN','COL2A1','COL9A1','COL9A2','COMP'],
'Erythroid':['HBA1','HBA2','HBB','ALAS2','GYPA']}
PANELS=json.loads((P/'phase16_scpca_supplement_20260920/resources/customer_panels.json').read_text())
PANELS['mCAF_matrix']=['COL1A1','COL1A2','COL3A1','FN1','POSTN','TNC','COL11A1']
PANELS['stress_descriptive']=['FOS','JUN','JUNB','ATF3','HSPA1A','HSPA1B','DDIT3']
PANELS['cycle_descriptive']=['MKI67','TOP2A','UBE2C','CENPF','TYMS','PCNA']
SURFACE=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
TARGET=sorted(set(sum(LINEAGE.values(),[])+sum(PANELS.values(),[])+SURFACE+['S100A8','S100A9','MYH11','ACTA2','TAGLN']))
protocol={
'created':now(),'user_instruction':'Review existing evidence and advance without waiting for memory-failed CopyKAT.',
'input':'All 107 full provider-called inputs, no added UMI/gene/mitochondrial/marker-number/reference-count thresholds.',
'CNV':'Use completed native/provider calls; freeze first-pass/retry provenance; active uncapped S525 recovery is nonblocking and not silently incorporated.',
'lineage':'Competitive mean log1p(CP10k) expression of predefined biological lineage panels. Unique positive top panel gives provisional lineage; exact tie/zero unassigned. Scores and gene coverage retained. This is operational annotation, not a validated classifier or customer numerical cutoff.',
'CAF_pool':'Provisional fibroblast lineage plus at least one native/provider no-alteration call, no opposing alteration call, no existing explicit software/curated doublet label. All other cells retained in audit.',
'annotation_review':'Review aggregate markers and original/provider annotations before adopting CAF clusters; provider or old labels alone are not inclusion gates. Osteogenic/mural dominant cells are separate.',
'known_doublets':'Official scDblFinder doublet or historical doublet_exclude, not a new numerical doublet threshold; unassessed status preserved.',
'patient_N':'Sample records are not independent patients. Merge same dataset/participant; no unverified cross-study pooling.',
'next':'CAF clustering and customer marker typing followed by myCAF internal reclustering; clinical contrasts only where metadata permit.',
'original_results_immutable':True}
writej(O/'PROTOCOL.json',protocol);writej(O/'resources/lineage_panels.json',LINEAGE);writej(O/'resources/program_panels.json',PANELS)
lock=(O/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
if (O/'REVIEW_COMPLETE').exists():raise SystemExit('Already completed; keep frozen result')
rc=1
try:
 manifest=read(B/'tables/analysis_source_manifest.tsv').fillna('');assert len(manifest)==107 and manifest['sample'].is_unique
 clinical=read(P/'phase14_client_objectives_20260920/tables/01_clinical_sample_audit_updated.tsv').fillna('').set_index('sample')
 parts=[];summ=[];coverage=[];clinic=[];provenance=[];background=[];evpaths=[]
 status('REVIEWING',finished_samples=0)
 for i,row in enumerate(manifest.to_dict('records')):
  sm=row['sample'];src=B/'inputs'/sm/'full_called_counts.h5ad'
  ef=B/'retry_once_20260924/tables'/f'{sm}_identity_evidence.tsv.gz'
  if not ef.exists():ef=B/'tables'/f'{sm}_identity_evidence.tsv.gz'
  ev=read(ef).set_index('cell_id').fillna('');a=ad.read_h5ad(src)
  assert a.obs_names.is_unique and len(a)==len(ev) and set(a.obs_names)==set(ev.index)
  ev=ev.loc[a.obs_names].copy();x=a.X.tocsr();den=np.asarray(x.sum(1)).ravel();assert np.allclose(den,ev.n_counts)
  symbols=a.var.symbol.astype(str);gene_map={}
  for j,sy in enumerate(symbols):
   for z in sy.split(';'):
    if z in TARGET:gene_map.setdefault(z,[]).append(j)
  raw=np.zeros((len(a),len(TARGET)),dtype=np.float32);measured=[]
  for j,g in enumerate(TARGET):
   ix=gene_map.get(g,[]);measured.append(bool(ix))
   if ix:raw[:,j]=np.asarray(x[:,ix].sum(1)).ravel()
  log=np.log1p(raw/np.maximum(den,1)[:,None]*1e4);scores={}
  for kind,panels in [('lineage',LINEAGE),('program',PANELS)]:
   for name,gs in panels.items():
    ix=[TARGET.index(g) for g in gs if g in gene_map]
    scores[name]=log[:,ix].mean(1) if ix else np.full(len(a),np.nan)
    ev[kind+'_'+name]=scores[name]
    coverage.append(dict(sample=sm,kind=kind,panel=name,measured=len(ix),requested=len(gs),missing=';'.join(g for g in gs if g not in gene_map)))
  v=np.column_stack([scores[k] for k in LINEAGE]);mx=np.nanmax(v,axis=1);label=np.asarray(list(LINEAGE))[np.nanargmax(v,axis=1)]
  label[(mx<=0)|(np.sum(v==mx[:,None],axis=1)!=1)]='Unresolved'
  ev['provisional_lineage']=label
  db=np.zeros(len(a),bool);dbstatus=np.full(len(a),'unassessed',object)
  if 'scDblFinder_class' in a.obs:
   dbstatus=a.obs.scDblFinder_class.astype(str).to_numpy();db|=dbstatus=='doublet'
  for c in ['doublet_exclude','historical_doublet']:
   if c in a.obs:
    vals=a.obs[c].astype(str).str.lower();db|=vals.isin(['true','1','1.0']).to_numpy()
    dbstatus=np.where(vals.isin(['true','false','1','0','1.0','0.0']),np.where(db,'doublet','historical_not_excluded'),dbstatus)
  ev['known_doublet']=db;ev['doublet_evidence']=dbstatus
  low=ev.CNV_evidence.isin(['single_method_no_alteration','concordant_no_alteration'])
  fib=ev.provisional_lineage.eq('Fibroblast')
  ev['review_class']='other_lineage'
  ev.loc[ev.provisional_lineage.eq('Unresolved'),'review_class']='lineage_unresolved'
  for category in ['single_method_alteration','concordant_alteration','method_conflict','unassessed']:
   ev.loc[fib&ev.CNV_evidence.eq(category),'review_class']='fibroblast_like_'+category
  ev.loc[fib&low,'review_class']='CAF_candidate_'+ev.loc[fib&low,'CNV_evidence']
  ev.loc[db,'review_class']='known_doublet'
  keep=fib&low&~db;ev['CAF_candidate']=keep
  meta={'sample':sm,'dataset':row['dataset'],'patient_id':row['patient_id'],'site':row['site'],'treatment':row['treatment'],'clinical_group':'Unknown','group_evidence':'not_established','primary_or_metastasis':'not_reported'}
  if sm in clinical.index:
   cr=clinical.loc[sm];meta.update(patient_id=cr.patient_verified,site=cr.site_verified,treatment=cr.treatment_evidence,clinical_group=cr.updated_comparison_role,group_evidence=cr.source_evidence)
  elif row['kind']=='ScPCA':
   timing=str(a.obs.disease_timing.iloc[0]) if 'disease_timing' in a.obs else 'not_reported';meta['primary_or_metastasis']=timing
   if 'lung' in row['site'].lower() and 'metasta' in timing.lower():meta.update(clinical_group='B2',group_evidence='provider tissue_location and disease_timing=Metastasis')
  meta['donor_record']=row['dataset']+'|'+str(meta['patient_id'])
  tech=str(a.obs.suspension_type.iloc[0]) if 'suspension_type' in a.obs else str(a.obs.technology.iloc[0]) if 'technology' in a.obs else 'not_reported'
  meta['technology']=tech;meta['study_technology']=row['dataset']+'|'+tech;clinic.append(meta)
  for k,z in meta.items():ev[k]=z
  ev['full_library_counts']=den
  ev['provider_annotation']=ev.get('provider_consensus_celltype_annotation','');ev['historical_annotation']=ev.get('historical_broad_lineage','')
  # All cells remain available for lineage/background/surface evidence; gene absence is explicit.
  bgobs=ev[['sample','dataset','donor_record','CNV_evidence','provisional_lineage','review_class','known_doublet','CAF_candidate']].copy()
  bg=ad.AnnData(sparse.csr_matrix(raw),obs=bgobs,var=pd.DataFrame(index=TARGET));bg.uns['measured_genes']=[g for g,m in zip(TARGET,measured) if m];clean(bg)
  bg.write_h5ad(O/'objects'/f'{sm}_marker_counts.h5ad',compression='gzip')
  out=O/'tables'/f'{sm}_review.tsv.gz';ev.to_csv(out,sep='\t',index_label='cell_id');evpaths.append(str(out))
  if keep.any():
   ids=a.var_names.astype(str);valid=np.array([z.startswith('ENSG') for z in ids]);assert pd.Index(ids[valid]).is_unique
   part=ad.AnnData(x[keep.to_numpy()][:,valid],obs=ev.loc[keep].copy(),var=a.var.loc[valid,['symbol']].copy());clean(part);parts.append(part)
  counts=ev.review_class.value_counts().to_dict();summ.append(dict(sample=sm,dataset=row['dataset'],input_cells=len(a),CAF_candidates=int(keep.sum()),dual_supported=int((keep&ev.CNV_evidence.eq('concordant_no_alteration')).sum()),single_supported=int((keep&ev.CNV_evidence.eq('single_method_no_alteration')).sum()),**counts))
  provenance.append(dict(sample=sm,input=str(src),input_bytes=src.stat().st_size,input_mtime_ns=src.stat().st_mtime_ns,evidence=str(ef),evidence_sha256=digest(ef)))
  del a,x,raw,log,bg;gc.collect()
  status('REVIEWING',finished_samples=i+1,CAF_candidates=sum(s['CAF_candidates'] for s in summ))
 clinical2=pd.DataFrame(clinic)
 # Bone samples from a documented participant who also contributed lung metastasis: history supported; contemporaneous timing not asserted.
 pairrows=[]
 for dn,g in clinical2.groupby('donor_record'):
  haslung=g.clinical_group.eq('B2').any()
  bone=g.site.str.contains('bone|femur|tibia|humerus|ulna|calcaneus',case=False,regex=True)
  if haslung and bone.any():
   ids=g.index[bone];clinical2.loc[ids,'clinical_group']='B1_history_supported_timing_unresolved';clinical2.loc[ids,'group_evidence']='same provider participant has lung metastasis specimen; sampling chronology not established'
   pairrows.append(dict(donor_record=dn,dataset=g.dataset.iloc[0],bone_samples=';'.join(g.loc[bone,'sample']),lung_samples=';'.join(g.loc[g.clinical_group.eq('B2'),'sample']),pair_type='same_participant_bone_lung_not_confirmed_contemporaneous'))
 save(clinical2,'clinical_mapping.tsv');save(pd.DataFrame(pairrows,columns=['donor_record','dataset','bone_samples','lung_samples','pair_type']),'bone_lung_pairs.tsv')
 save(pd.DataFrame(summ).fillna(0),'sample_review_counts.tsv');save(pd.DataFrame(coverage),'gene_coverage.tsv');writej(O/'resources/SOURCE_SNAPSHOT.json',provenance);writej(O/'resources/review_tables.json',evpaths)
 common=sorted(set.intersection(*[set(a.var_names) for a in parts]));assert common
 obs=pd.concat([a.obs for a in parts]);assert obs.index.is_unique
 cm=clinical2.set_index('sample')
 for c in ['clinical_group','group_evidence']:obs[c]=obs['sample'].map(cm[c])
 v=parts[0].var.loc[common].copy()
 pooled=ad.AnnData(sparse.vstack([a[:,common].X for a in parts],format='csr'),obs=obs,var=v);clean(pooled)
 pooled.write_h5ad(O/'objects/CAF_candidates_counts.h5ad',compression='gzip')
 qa={'at':now(),'sample_records':len(summ),'all_input_cells':int(sum(z['input_cells'] for z in summ)),'CAF_candidates':pooled.n_obs,'candidate_samples':pooled.obs['sample'].nunique(),'candidate_datasets':pooled.obs.dataset.nunique(),'dual_supported_candidates':int(pooled.obs.CNV_evidence.eq('concordant_no_alteration').sum()),'single_supported_candidates':int(pooled.obs.CNV_evidence.eq('single_method_no_alteration').sum()),'shared_measured_genes':len(common),'CNV_conflict_in_main':int(pooled.obs.CNV_evidence.eq('method_conflict').sum()),'unassessed_in_main':int(pooled.obs.CNV_evidence.eq('unassessed').sum()),'known_doublets_in_main':int(pooled.obs.known_doublet.sum()),'source_preserved':True,'S525_uncapped_recovery':'not_waited_for; separate running branch; frozen evidence uses completed previous calls','cell_identity_level':'operational CAF candidates, cluster marker review required','candidate_bone_lung_donor_records':len(pairrows)}
 assert qa['all_input_cells']==876648 and qa['CNV_conflict_in_main']==qa['unassessed_in_main']==qa['known_doublets_in_main']==0
 writej(O/'delivery/IDENTITY_QA.json',qa);(O/'REVIEW_COMPLETE').write_text(now()+'\n');status('IDENTITY_EVIDENCE_REVIEW_COMPLETE_CAF_CLUSTERING_READY',**qa);rc=0
except Exception as e:status('FAILED',error=str(e));traceback.print_exc()
finally:(O/'exit_code').write_text(str(rc)+'\n')
raise SystemExit(rc)
