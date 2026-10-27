from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import pandas as pd,numpy as np,json
from scipy.stats import hypergeom,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');C=P/'phase43_full_transcriptome_discovery_20261001';R=P/'phase49_audit_revision_20261002'
def bh(p):
 v=np.asarray(p,float);q=false_discovery_control(np.where(np.isfinite(v),v,1));q[~np.isfinite(v)]=np.nan;return q
cross=pd.read_csv(C/'tables/HGNC_full_gene_symbol_crosswalk.tsv.gz',sep='	',keep_default_na=False);feature_priority=dict(zip(cross.feature,np.where(cross.mapping_method.eq('HGNC_Ensembl_approved'),0,1)))
reg=pd.read_csv(R/'tables/pathway_input_revision_registry.tsv',sep='	');mp=pd.read_csv(R/'resources/source_canonical_feature_map.tsv',sep='	');libraries=json.loads((R/'resources/libraries_canonical.json').read_text());checks=[];cached=None
for rr in reg.itertuples():
 if rr.N_matched_patients<2 or rr.new_background<2:continue
 if rr.gene_result_file!=cached:
  gene=pd.read_csv(C/rr.gene_result_file,sep='	',low_memory=False,float_precision='round_trip');cached=rr.gene_result_file
  expected_map={}
  for en,sy in gene[['feature','symbol']].drop_duplicates().itertuples(index=False,name=None):
   if pd.isna(sy) or not sy:continue
   candidate=(feature_priority.get(en,1),en)
   if sy not in expected_map or candidate<expected_map[sy]:expected_map[sy]=candidate
 g=gene[gene.contrast.eq(rr.contrast)];mapping=mp[mp.layer.eq(rr.layer)&mp.dataset.eq(rr.dataset)];assert not mapping.symbol.duplicated().any();assert dict(zip(mapping.symbol,mapping.feature))=={k:v[1] for k,v in expected_map.items()};chosen=g[g.feature.isin(mapping.feature)];chosen=chosen[chosen.P.notna()&chosen.N.ge(2)];bg=set(chosen.symbol);M=len(bg)
 rank=pd.read_csv(R/'resources'/f'{rr.context_id}__rank.tsv',sep='	',float_precision='round_trip');assert set(rank.symbol)==bg and M==rr.new_background
 # Independent membership and hypergeometric-tail reconstruction for every term.
 for lib,sets in libraries.items():
  f=R/'tables/pathways'/f'{rr.context_id}__{lib}__ORA.tsv.gz';t=pd.read_csv(f,sep='	',float_precision='round_trip');pv=[]
  for direction,z in t.groupby('direction',sort=False):
   q=set(chosen.loc[chosen.P.lt(.05)&(chosen.mean_delta.gt(0) if direction=='up' else chosen.mean_delta.lt(0)),'symbol']);assert len(q)==int(z.query_genes.iloc[0]);vals=[]
   for r in z.itertuples():
    members=set(sets[r.term])&bg;k=len(members&q);K=len(members);n=len(q);assert k==r.overlap_genes and K==r.set_background_genes and M==r.background
    vals.append(float(hypergeom.cdf(n-k,M,M-K,n)) if n and K else np.nan)
   vals=np.asarray(vals);assert np.array_equal(np.isnan(vals),z.P.isna());pe=float(np.nanmax(np.r_[0,np.abs(vals-z.P)]));qe=float(np.nanmax(np.r_[0,np.abs(bh(vals)-z.FDR)]));assert pe<1e-10 and qe<1e-10;checks.append(dict(context=rr.context_id,library=lib,kind='ORA_'+direction,rows=len(z),P_error=pe,FDR_error=qe,passed=True))
  gf=R/'tables/pathways'/f'{rr.context_id}__{lib}__GSEA.tsv.gz'
  if gf.exists():
   gt=pd.read_csv(gf,sep='	',float_precision='round_trip');err=float(np.nanmax(np.r_[0,np.abs(bh(gt.GSEA_model_P)-gt.model_BH_FDR)]));assert err<1e-10;checks.append(dict(context=rr.context_id,library=lib,kind='GSEA_model_BH',rows=len(gt),P_error=np.nan,FDR_error=err,passed=True))
pd.DataFrame(checks).to_csv(R/'tables/pathway_independent_QA.tsv',sep='	',index=False);out={'files_and_directions':len(checks),'rows':sum(x['rows'] for x in checks),'all_passed':True,'ORA_P_max_error':max(x['P_error'] for x in checks if np.isfinite(x['P_error'])),'FDR_max_error':max(x['FDR_error'] for x in checks),'GSEA_P_scope':'output range/native model P provenance plus independently recomputed BH; stochastic model P not independently replayed'};(R/'delivery/PATHWAY_INDEPENDENT_QA.json').write_text(json.dumps(out,indent=2));(R/'PATHWAY_INDEPENDENT_QA_COMPLETE').write_text('SUCCESS');print(json.dumps(out))
