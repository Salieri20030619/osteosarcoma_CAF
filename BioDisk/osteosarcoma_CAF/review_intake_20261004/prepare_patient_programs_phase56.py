from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,hashlib,datetime,gzip,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');S=P/'phase53_de_novo_reannotation_20261002';A=P/'phase55_annotation_freeze_20261004';O=P/'phase56_patient_programs_20261004';OLD=P/'phase41_full_transcriptome_setup_20261001'
for s in ['code','tables','resources','objects','logs','delivery','figures']:(O/s).mkdir(parents=True,exist_ok=True)
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
assert (A/'freeze.exit_code').read_text().strip()=='0' and (A/'ANNOTATION_FREEZE_COMPLETE').is_file()
ann=read(A/'tables/all_cell_annotation_v2.tsv.gz');assert len(ann)==876648 and ann.cell_id.is_unique;scope=list(json.loads((A/'resources/ANNOTATION_AND_SCOPE_PROTOCOL_FROZEN.json').read_text())['scope_definitions']);groups=scope+['fixed_background','osteogenic_reference','all_input']
cw=read(P/'phase43_full_transcriptome_discovery_20261001/tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',keep_default_na=False).set_index('feature');hgnc=read(OLD/'resources/hgnc_complete_set.txt',dtype=str,keep_default_na=False);hgnc=hgnc[hgnc.status.eq('Approved')];approved=set(hgnc.symbol);aliases={}
for row in hgnc.itertuples(index=False):
 for field in ['prev_symbol','alias_symbol']:
  for v in getattr(row,field,'').split('|'):
   if v:aliases.setdefault(v,set()).add(row.symbol)
def canon(g):return g if g in approved else (next(iter(aliases[g])) if len(aliases.get(g,[]))==1 else g)
sets={};origin={}
for lib,name in [('hallmark','hallmark.gmt'),('matrisome','curated_cgp.gmt')]:
 for line in (OLD/'resources'/name).read_text().splitlines():
  z=line.split('\t')
  if lib=='matrisome' and not z[0].startswith('NABA_'):continue
  sets[z[0]]=list(dict.fromkeys(canon(g) for g in z[2:]));origin[z[0]]=lib
functional=json.loads((P/'phase46_spatial_functional_20261001/resources/FUNCTIONAL_SIGNATURES_FROZEN.json').read_text())['programs']
for name,genes in functional.items():
 if name not in sets:sets[name]=list(dict.fromkeys(canon(g) for g in genes));origin[name]='previously_frozen_functional'
custom=json.loads((P/'phase31_customer_completion_20260926/resources/customer_program_panels.json').read_text())
for name,genes in custom.items():sets['CUSTOM_'+name]=list(dict.fromkeys(canon(g) for g in genes));origin['CUSTOM_'+name]='previously_frozen_customer'
dump(O/'resources/PROGRAMS_FROZEN.json',{'programs':sets,'origin':origin})
proto={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'annotation_SHA256':sha(A/'tables/all_cell_annotation_v2.tsv.gz'),'scopes':scope,'groups':groups,'programs':len(sets),'primary_reference':'fixed_background','secondary_reference':'osteogenic_reference','patient_unit':'source-scoped donor; same-site pseudobulk before equal-site weighting; no source pooling','genes':'Deterministic HGNC-preferred feature per canonical symbol from each source common measured universe; no global intersection. Full-library denominator retains all raw/corrected features.','raw_and_DecontX_programs':True,'test':'Two-sided mean patient signflip; exact N<=15, 19999 random flips + observed otherwise; n<2 and all-zero effect not tested.','pseudocount_CP10k':.1,'FDR':'BH within source|scope|reference|normalization|programs; gene and ORA families separate. All predefined program slots retained; non-evaluable P recorded NA.','selection':'Raw P<0.05; FDR reported not mandatory for candidate retention. N<10 exploratory.','sensitivity':'Compare all four scopes and raw/corrected programs; common-patient comparisons exported. No choice of scope based on smallest P.','primary_display':'agreed11, plus explicit CNV flag strata; no_flag is sensitivity, not proven nonmalignant.','clinical_or_CAF_identity_claim':'Expression association/discovery only; same-data annotation selection not removed by FDR.','no_new_cell_QC_threshold':True}
dump(O/'resources/PROTOCOL_FROZEN.json',proto)
featuremaps={};coverage=[]
for source,g in ann.groupby('source_id'):
 gf=read(S/'blind_review'/source/'genes.tsv',keep_default_na=False);choices={}
 for f in gf.feature.astype(str):
  symbol=cw.at[f,'canonical_symbol'] if f in cw.index else ''
  if not symbol:continue
  priority=0 if f in cw.index and cw.at[f,'mapping_method']=='HGNC_Ensembl_approved' else 1
  k=(priority,f)
  if symbol not in choices or k<choices[symbol][0]:choices[symbol]=(k,f)
 genes=sorted(choices);features=[choices[z][1] for z in genes];featuremaps[source]=(genes,features)
 pd.DataFrame({'gene':genes,'feature_id':features}).to_csv(O/'resources'/f'{source}_gene_feature_mapping.tsv',sep='\t',index=False)
 for name,members in sets.items():coverage.append({'source_id':source,'program':name,'collection':origin[name],'defined_genes':len(members),'measured_genes':sum(g in choices for g in members),'missing_genes':';'.join(g for g in members if g not in choices)})
pd.DataFrame(coverage).to_csv(O/'tables/program_gene_coverage.tsv',sep='\t',index=False)
meta=[];inputqa=[]
for i,(sample,g) in enumerate(ann.groupby('sample',sort=True)):
 source=g.source_id.iloc[0];sid=g.sample_id.iloc[0];genes,features=featuremaps[source];f=S/'blind_review/raw_counts'/f'{sid}.h5ad';a=ad.read_h5ad(f);g=g.set_index('cell_id').loc[a.obs_names];x=a.X.tocsr().astype(float);rawlib=np.asarray(x.sum(1)).ravel();assert np.array_equal(rawlib,g.full_library.to_numpy())
 base=P/'phase26_client_reanalysis_20260922/inputs'/sample;cells=read(base/'cells.tsv');ig=read(base/'genes.tsv');assert cells.cell_id.tolist()==g.original_cell_id.tolist() and ig.ensembl.tolist()==a.var_names.tolist()
 de=P/'phase31_customer_completion_20260926/decontx'/sample;assert (de/'COMPLETE').is_file() and (de/'exit_code').read_text().strip()=='0';cc=read(de/'cell_contamination.tsv');assert cc.cell_id.tolist()==cells.cell_id.tolist()
 with gzip.open(de/'corrected_counts.mtx.gz','rb') as h:
  with threadpool_limits(limits=4):y=mmread(h).tocsr().T.tocsr()
 corlib=np.asarray(y.sum(1)).ravel();assert y.shape==x.shape and np.allclose(corlib,cc.corrected_UMI,rtol=1e-5,atol=1e-4)
 ix=a.var_names.get_indexer(features);assert (ix>=0).all();raw=[];corr=[];rl=[];cl=[];ns=[]
 for group in groups:
  if group=='all_input':mask=np.ones(len(g),dtype=bool)
  elif group in scope:mask=g['scope_'+group].to_numpy(bool)
  else:mask=g[group].to_numpy(bool)
  raw.append(np.asarray(x[mask][:,ix].sum(0)).ravel());corr.append(np.asarray(y[mask][:,ix].sum(0)).ravel());rl.append(rawlib[mask].sum());cl.append(corlib[mask].sum());ns.append(int(mask.sum()))
 np.savez_compressed(O/'objects'/f'{sid}_pseudobulk.npz',raw=np.vstack(raw),corrected=np.vstack(corr),raw_library=rl,corrected_library=cl,cells=ns,groups=np.asarray(groups),genes=np.asarray(genes),features=np.asarray(features))
 meta.append({'dataset':g.dataset.iloc[0],'source_id':source,'sample':sample,'sample_id':sid,'patient_id':g.patient_id.iloc[0],'donor_record':g.original_donor_record.iloc[0],'site':g.site.iloc[0],'object':str(O/'objects'/f'{sid}_pseudobulk.npz')})
 inputqa.append({'sample_id':sid,'sample':sample,'cells':len(g),'canonical_genes':len(genes),'cell_and_feature_order_verified':True,'raw_full_library_verified':True,'corrected_full_library_verified':True,'agreed11_cells':ns[0],'all15_cells':ns[2]})
 print(json.dumps({'pseudobulk_samples':i+1,'total':107,'sample':sid}),flush=True);dump(O/'PROGRESS.json',{'stage':'pseudobulk','completed_samples':i+1,'total_samples':107});del a,x,y;gc.collect()
pd.DataFrame(meta).to_csv(O/'tables/pseudobulk_sample_metadata.tsv',sep='\t',index=False);pd.DataFrame(inputqa).to_csv(O/'tables/pseudobulk_input_QA.tsv',sep='\t',index=False)
assert sum(q['cells'] for q in inputqa)==876648 and sum(q['agreed11_cells'] for q in inputqa)==45721 and sum(q['all15_cells'] for q in inputqa)==91012
dump(O/'delivery/PSEUDOBULK_ACCEPTANCE.json',{'samples':107,'all_cells':876648,'programs':len(sets),'scopes':scope,'gene_universe':'source-common measured canonical genes','raw_and_corrected_complete':True});(O/'PSEUDOBULK_COMPLETE').write_text('SUCCESS\n');print('PSEUDOBULK_COMPLETE',flush=True)