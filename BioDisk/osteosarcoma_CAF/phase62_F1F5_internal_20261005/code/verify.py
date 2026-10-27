from pathlib import Path
import json,itertools,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase62_F1F5_internal_20261005';B=P/'phase58_F1F5_round1_20261004';A=P/'phase41_full_transcriptome_setup_20261001'
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
assert (O/'ANALYSIS_COMPLETE').exists()
ob=read(O/'tables/cells_original_with_review.tsv.gz').set_index('cell_id');orig=read(B/'tables/frozen_expanded_cells.tsv.gz').set_index('cell_id');labels=read(O/'tables/F1F5_internal_cell_labels.tsv.gz').set_index('cell_id');assert len(labels)==8174 and labels.index.is_unique and labels.index.tolist()==orig.index.tolist()
for col in ['round1_state','class_id']:assert labels[col].equals(orig[col])
for col in ['lineage_confidence','candidate_tier','CNV_evidence_current','mixed_program_flag','doublet_assessed_any','doublet_positive_any']:assert labels[col].fillna('NA').equals(ob[col].fillna('NA'))
assert labels.subgroup_id.nunique()==33 and labels.subgroup_id.notna().all();assert labels.expression_confidence.notna().all() and labels.cross_patient_support.notna().all()
cr=read(O/'tables/contrast_registry.tsv');mp=read(B/'resources/source_canonical_feature_map.tsv');probes=read(O/'tables/marker_raw_corrected_evidence.tsv');checked=0;maxerr=0
for mode,path in [('raw',A/'objects/all_candidate_full_genes.h5ad'),('corrected',B/'objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);den=np.asarray(a.X.sum(1)).ravel();assert a.obs_names.tolist()==ob.index.tolist()
 for st,pr in probes.groupby('subgroup_id'):
  ds=pr.primary_source.iloc[0];cl=st.split('_')[0];left=np.flatnonzero(ob.subgroup_id.eq(st)&ob.dataset.eq(ds));right=np.flatnonzero(ob.class_id.eq(cl)&ob.subgroup_id.ne(st)&ob.dataset.eq(ds));genes=pr.symbol.head(5).tolist();features=mp[mp.dataset.eq(ds)].set_index('symbol').loc[genes,'feature'];jj=a.var_names.get_indexer(features);lx=a.X[left][:,jj].toarray();rx=a.X[right][:,jj].toarray();calc=np.log2(lx.sum(0)/den[left].sum()*1e4+.1)-np.log2(rx.sum(0)/den[right].sum()*1e4+.1);expected=pr.set_index('symbol').loc[genes,mode+'_effect'].to_numpy();err=float(np.max(np.abs(calc-expected)));assert err<1e-6;maxerr=max(maxerr,err);checked+=len(genes)
 del a
stats=read(O/'tables/paired_gene_statistics.tsv.gz');families=0;bhmax=0
for fam,g in stats.groupby('family'):
 q=false_discovery_control(g.P.fillna(1).to_numpy(),method='bh');valid=g.P.notna().to_numpy();err=np.max(np.abs(q[valid]-g.loc[g.P.notna(),'FDR'])) if valid.any() else 0;assert err<1e-10;bhmax=max(bhmax,float(err));families+=1
assert stats.loc[stats.N.lt(2),'P'].isna().all();assert stats.nominal.equals(stats.P.lt(.05))
pchecks=0
for p in (O/'objects').glob('*__paired.npz'):
 z=np.load(p);n=len(z['effects']);ds,st=p.name.split('__')[:2]
 if n<2 or n>15:continue
 ix=np.linspace(0,len(z['symbols'])-1,7,dtype=int);d=z['effects'][:,ix];sgn=np.array(list(itertools.product([-1,1],repeat=n)));expected=(np.abs(sgn@d/n)>=np.abs(d.mean(0))[None,:]-1e-12).mean(0);f=stats[(stats.dataset==ds)&(stats.subgroup_id==st)].set_index('symbol');got=f.loc[z['symbols'][ix],'P'].to_numpy();active=np.any(d!=0,axis=0);assert np.allclose(expected[active],got[active],atol=1e-12);assert np.isnan(got[~active]).all();pchecks+=int(active.sum())
jobs=json.loads((O/'resources/GSEA_JOBS.json').read_text());gres=0
for j in jobs:
 g=read(j['output']);rank=read(j['rank']);bg=set(rank.symbol);assert len(rank)==len(bg) and len(g)>0;assert g.GSEA_model_P.between(0,1).all();q=false_discovery_control(g.GSEA_model_P.to_numpy(),method='bh');assert np.allclose(q,g.model_BH_FDR,atol=1e-12)
 for lead in g.Lead_genes.fillna(''):assert set(lead.split(';'))-{''}<=bg
 gres+=len(g)
g=read(O/'tables/all_GSEA_results.tsv.gz');assert len(g)==gres
report={'status':'PASS','cells':len(labels),'subgroups':labels.subgroup_id.nunique(),'original_labels_and_confidence_preserved':True,'gene_effect_probes':checked,'maximum_gene_effect_error':maxerr,'patient_P_probes':pchecks,'gene_BH_families':families,'maximum_BH_error':bhmax,'GSEA_jobs':len(jobs),'GSEA_rows':gres,'all_leading_edges_in_rank':True,'exit_codes':{stage:(O/f'{stage}.exit_code').read_text().strip() for stage in ['prepare','gsea']}}
assert set(report['exit_codes'].values())=={'0'}
(O/'QA/ACCEPTANCE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(O/'VERIFIED').write_text('SUCCESS\n');print(json.dumps(report,indent=2),flush=True)
