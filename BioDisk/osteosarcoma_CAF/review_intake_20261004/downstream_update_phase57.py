from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import sys,json,hashlib,datetime,importlib.util
import numpy as np,pandas as pd
from scipy.stats import hypergeom
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase56_patient_programs_20261004';O=P/'phase57_downstream_update_20261004';OLD=P/'phase41_full_transcriptome_setup_20261001'
for sub in ['code','tables/genes','tables/pathways','resources','delivery','logs','figures']:(O/sub).mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('patient_stats',B/'code/patient_program_statistics.py');st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
protocol={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'target_contrasts':'All four frozen scopes vs the same fixed noncandidate/nonmixed background, raw counts; no scope chosen from P values.','genes':'All source-common deterministic canonical genes used in phase56 pseudobulk; effects are conditional on selected expression-defined populations, not independent identity validation.','gene_statistic':'Patient mean log2 CP10k effect; exact/Monte Carlo signflip as phase56.','gene_FDR_family':'One source|scope|reference|raw gene family.','ORA_selection':'Raw gene P<0.05; sign of mean effect, up and down separately. No FDR entry gate.','ORA_background':'All validly tested canonical genes after the pre-test source mapping; no replacement feature chosen when a gene is untestable.','ORA_families':'source|scope|library|direction; Hallmark, NABA Matrisome, Reactome separately; all predefined terms retained.','ORA_P_type':'Overrepresentation model P; not additional independent patient-level validation.','spatial':'Existing fixed gene-set/spot correlations depend on unchanged spatial data and score definitions; retain numeric results, update RNA-source interpretation from new candidate shares. No spatial remapping or boundary inference.','historical_F3_F4':'Retain as operational-label historical analysis; do not transfer old labels, effects or signature validation claims to new fibroblast candidates.'}
dump(O/'resources/DOWNSTREAM_PROTOCOL_FROZEN.json',protocol)
assert (B/'PROGRAMS_COMPLETE').is_file() and (B/'programs.exit_code').read_text().strip()=='0'
registry=read(B/'tables/executed_contexts.tsv');registry=registry[registry.normalization.eq('raw')&registry.reference.eq('fixed_background')]
hg=pd.read_csv(OLD/'resources/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);hg=hg[hg.status.eq('Approved')];approved=set(hg.symbol);aliases={}
for row in hg.itertuples(index=False):
 for field in ['prev_symbol','alias_symbol']:
  for a in getattr(row,field,'').split('|'):
   if a:aliases.setdefault(a,set()).add(row.symbol)
def canon(g):return g if g in approved else (next(iter(aliases[g])) if len(aliases.get(g,[]))==1 else g)
libs={}
for lib,fn in [('hallmark','hallmark.gmt'),('matrisome','curated_cgp.gmt'),('reactome','reactome.gmt')]:
 d={}
 for line in (OLD/'resources'/fn).read_text().splitlines():
  z=line.split('\t')
  if lib=='matrisome' and not z[0].startswith('NABA_'):continue
  d[z[0]]=set(canon(g) for g in z[2:])
 libs[lib]=d
cache={};summary=[];orasummary=[];nominal=[]
for ctx in registry.itertuples(index=False):
 if not ctx.N:summary.append({'context_id':ctx.context_id,'source_id':ctx.source_id,'scope':ctx.scope,'N':0,'status':'no_matched_patients'});continue
 inp=np.load(B/'objects'/f'{ctx.context_id}_gene_effects.npz');d=inp['effects'];genes=inp['genes'].astype(str);n=d.shape[0];signature=(ctx.source_id,hashlib.sha256(d.tobytes()+inp['donors'].tobytes()).hexdigest())
 if signature in cache:pv,method=cache[signature]
 else:pv,method=st.tests(d);cache[signature]=(pv,method)
 mu=d.mean(0);q=st.bh(pv);nz=(d!=0).sum(0);floor=np.where(nz>0,2.**(1-nz),np.nan) if n<=15 else np.full(len(genes),1/20000)
 out=pd.DataFrame({'source_id':ctx.source_id,'dataset':ctx.dataset,'scope':ctx.scope,'reference':ctx.reference,'normalization':'raw','gene':genes,'N':n,'mean_effect':mu,'median_effect':np.median(d,axis=0),'positive_patients':(d>0).sum(0),'negative_patients':(d<0).sum(0),'P':pv,'FDR':q,'test':method,'nominal':pv<.05,'FDR_significant':q<.05,'family':ctx.context_id+'|genes','minimum_P_resolution':floor,'evidence_level':'exploratory_small_N' if n<10 else 'source_scoped_discovery'})
 out.to_csv(O/'tables/genes'/f'{ctx.context_id}.tsv.gz',sep='\t',index=False);hit=out[out.nominal];nominal.append(hit)
 summary.append({'context_id':ctx.context_id,'source_id':ctx.source_id,'dataset':ctx.dataset,'scope':ctx.scope,'N':n,'status':'inference_available' if n>=2 else 'direction_only','genes':len(out),'valid_tests':int(np.isfinite(pv).sum()),'nominal_up':int(((pv<.05)&(mu>0)).sum()),'nominal_down':int(((pv<.05)&(mu<0)).sum()),'FDR_genes':int((q<.05).sum())})
 background=set(out.loc[out.P.notna(),'gene']);M=len(background)
 for lib,sets in libs.items():
  for direction in ['up','down']:
   selected=set(hit.loc[hit.mean_effect.gt(0) if direction=='up' else hit.mean_effect.lt(0),'gene']);nn=len(selected);rr=[]
   for term,defined in sets.items():
    members=defined&background;K=len(members);overlap=selected&members;kk=len(overlap);p=float(hypergeom.sf(kk-1,M,K,nn)) if M>1 and nn>0 and K>0 and K<M else np.nan
    rr.append({'context_id':ctx.context_id,'source_id':ctx.source_id,'dataset':ctx.dataset,'scope':ctx.scope,'library':lib,'direction':direction,'term':term,'N_patients':n,'background_genes':M,'query_genes':nn,'set_background_genes':K,'defined_genes':len(defined),'overlap_genes':kk,'fold_enrichment':(kk/nn)/(K/M) if M and nn and K else np.nan,'P':p,'overlap_symbols':';'.join(sorted(overlap)),'family':ctx.context_id+'|ORA|'+lib+'|'+direction,'P_type':'gene_set_overrepresentation_model'})
   df=pd.DataFrame(rr);df['FDR']=st.bh(df.P);df['nominal']=df.P.lt(.05);df.to_csv(O/'tables/pathways'/f'{ctx.context_id}__{lib}__{direction}.tsv.gz',sep='\t',index=False);orasummary.append({'context_id':ctx.context_id,'source_id':ctx.source_id,'scope':ctx.scope,'library':lib,'direction':direction,'N':n,'background_genes':M,'query_genes':nn,'valid_tests':int(df.P.notna().sum()),'nominal_terms':int(df.nominal.sum()),'FDR_terms':int(df.FDR.lt(.05).sum())})
 print(json.dumps(summary[-1]),flush=True)
pd.DataFrame(summary).to_csv(O/'tables/gene_statistics_summary.tsv',sep='\t',index=False);pd.DataFrame(orasummary).to_csv(O/'tables/pathway_statistics_summary.tsv',sep='\t',index=False)
if nominal:pd.concat(nominal,ignore_index=True).to_csv(O/'tables/all_nominal_genes.tsv.gz',sep='\t',index=False)
# Retained spatial numeric results are copied as references with provenance; their phenotype interpretation is updated later in the report.
spatialfiles=[P/'phase46_spatial_functional_20261001/tables/spatial_patient_statistics.tsv',P/'phase50_F3_spatial_decomposition_20261002/tables/statistics.tsv']
links=[]
for f in spatialfiles:
 x=read(f);dest=O/'tables'/('retained_'+f.parent.parent.name+'_statistics.tsv');x.to_csv(dest,sep='\t',index=False);links.append({'original_path':str(f),'original_SHA256':hashlib.sha256(f.read_bytes()).hexdigest(),'copy':str(dest),'numeric_recomputation':False,'reason':'Fixed spatial gene sets and spatial inputs unchanged; annotation changes alter interpretation, not the existing correlation calculation.'})
dump(O/'resources/RETAINED_SPATIAL_PROVENANCE.json',links)
impact=[{'module':'annotation','action':'updated','current':str(P/'phase55_annotation_freeze_20261004/tables/all_cell_annotation_v2.tsv.gz'),'interpretation':'All cells preserved; expression lineage, mixed risk and CNV flags separate.'},{'module':'patient_distribution','action':'updated','current':str(P/'phase55_annotation_freeze_20261004/tables/patient_candidate_distribution.tsv'),'interpretation':'Source-scoped patient distributions; counts not independent N.'},{'module':'patient_programs','action':'recomputed_four_scopes_raw_and_corrected','current':str(B/'tables/all_patient_program_statistics.tsv.gz'),'interpretation':'Expression-defined candidate populations; nominal P and FDR retained.'},{'module':'gene_differences','action':'recomputed_new_populations','current':str(O/'tables/genes'),'interpretation':'New contrasts replace use of historical F3/F4 as confirmed CAF; original contrasts archived.'},{'module':'pathways','action':'recomputed_new_gene_sets','current':str(O/'tables/pathways'),'interpretation':'Nominal gene selection; ORA model evidence separate from patient program evidence.'},{'module':'spatial_fixed_programs','action':'numeric_results_retained_attribution_updated','current':str(O/'resources/RETAINED_SPATIAL_PROVENANCE.json'),'interpretation':'Matrix/RNA spatial association; no automatic attribution to confirmed CAF or boundary geometry.'},{'module':'old_F3_F4_TF_LR_clinical_external_surface','action':'historical_operational_label_results_not_transferred','current':str(P/'phase49_audit_revision_20261002/delivery/client/CURRENT_RESULT_INDEX.tsv'),'interpretation':'No new identity, clinical or external validation claim; optional specialized reanalysis is outside this request sequence.'}]
pd.DataFrame(impact).to_csv(O/'tables/DOWNSTREAM_RESULT_INDEX.tsv',sep='\t',index=False)
qa={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'gene_contexts':len(summary),'ORA_context_library_directions':len(orasummary),'all_nominal_genes_retained':True,'direction_uses_mean_effect':True,'mapping_frozen_before_tests':True,'spatial_numbers_preserved_with_provenance':True,'original_labels_and_baseline_not_overwritten':True};dump(O/'delivery/DOWNSTREAM_ACCEPTANCE.json',qa);(O/'DOWNSTREAM_COMPLETE').write_text('SUCCESS\n');print('DOWNSTREAM_COMPLETE',flush=True)