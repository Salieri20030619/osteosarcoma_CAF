from pathlib import Path
import os,sys,json,gc,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from current_metadata import load_current_cells
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase66_doublet_repair_DSC1_20261005';A=B/'phase41_full_transcriptome_setup_20261001';F=B/'phase58_F1F5_round1_20261004'
ob=load_current_cells().set_index('cell_id');genes=['DSC1','DSC2','DSC3','DSG1','DSG2','DSG3','DSG4','DSP','JUP','PKP1','PKP2','PKP3']
h=pd.read_csv(A/'resources/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
exact=h[h.apply(lambda r:any(v.upper()=='SR1004' for col in ['symbol','prev_symbol','alias_symbol'] for v in r[col].split('|')),axis=1)]
h[h.symbol.isin(genes)].to_csv(O/'resources/desmosome_HGNC_records.tsv',sep='\t',index=False)
ident={'query':'sr1004','query_kind':'antibody_clone','gene':'DSC1','protein':'Desmocollin-1','catalog':'NBP3-22640','official_product':'https://www.rndsystems.com/jp/products/desmocollin-1-antibody-sr1004_nbp3-22640','manufacturer_catalog':'https://www.novusbio.com/primary-antibodies/desmocollin-1','NCBI':'https://www.ncbi.nlm.nih.gov/gene/1823/','Entrez':'1823','HGNC':'HGNC:3035','Ensembl':'ENSG00000134765','UniProt':'Q08554','SR1004_HGNC_exact_matches':len(exact),'product_listed_application':'WB; other DSC1 clones on catalog are distinct reagents and their applications are not transferred to SR1004','interpretation':'Desmocollin member of desmosomal cadherins; clone code is not a transcript feature; do not automatically label DSC1 as CAF-specific or functional mediator.'}
(O/'resources/SR1004_identity.json').write_text(json.dumps(ident,indent=2))
mp=pd.read_csv(F/'resources/source_canonical_feature_map.tsv',sep='\t');presence=sparse.load_npz(A/'objects/sample_feature_presence.npz').tocsr();sidx=pd.read_csv(A/'tables/feature_presence_sample_index.tsv',sep='\t').set_index('sample').presence_row.to_dict()
frames=[];coverage=[]
for mode,path in [('raw',A/'objects/all_candidate_full_genes.h5ad'),('corrected',F/'objects/expanded_corrected_full_genes.h5ad')]:
 a=ad.read_h5ad(path);assert list(a.obs_names)==list(ob.index)
 x=a.X.tocsr();lib=np.asarray(x.sum(1)).ravel();allv=np.full((len(ob),len(genes)),np.nan)
 for ds,t in ob.groupby('dataset'):
  mapping=mp[mp.dataset==ds].set_index('symbol').feature.to_dict();ids=ob.index.get_indexer(t.index)
  for j,g in enumerate(genes):
   feature=mapping.get(g);idx=a.var_names.get_indexer([feature or '__absent__'])[0]
   for sm,smobs in t.groupby('sample'):
    ix=ob.index.get_indexer(smobs.index);assayed=idx>=0 and bool(presence[sidx[sm],idx])
    if assayed:allv[ix,j]=x[ix,idx].toarray().ravel()
    coverage.append({'matrix':mode,'dataset':ds,'sample':sm,'gene':g,'feature':feature or '', 'assayed':assayed,'cells':len(ix)})
 f=ob[['sample','dataset','donor_record','class_id','subgroup_id','lineage_stratum','doublet_positive_any','lineage_confidence']].copy();f['matrix']=mode;f['full_library']=lib
 for j,g in enumerate(genes):f[g]=allv[:,j]
 frames.append(f.reset_index());del a,x;gc.collect()
d=pd.concat(frames,ignore_index=True);d.to_csv(O/'tables/desmosome_cell_counts.tsv.gz',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/desmosome_feature_coverage.tsv',sep='\t',index=False)
summ=[]
for arm in ['keep_all','exclude_model_positive']:
 z=d if arm=='keep_all' else d[~d.doublet_positive_any]
 for level in ['class_id','subgroup_id','dataset']:
  for (mode,k),t in z.groupby(['matrix',level]):
   for g in genes:
    valid=t[g].notna();v=t.loc[valid,g]
    summ.append({'arm':arm,'matrix':mode,'level':level,'group':k,'gene':g,'cells':len(t),'assayed_cells':int(valid.sum()),'RNA_ge1_cells':int(sum(v>=1)),'RNA_gt0_cells':int(sum(v>0)),'detection_ge1_fraction':float((v>=1).mean()) if len(v) else np.nan,'positive_donor_records':t.loc[valid&(t[g]>=1),'donor_record'].nunique(),'total_counts':float(v.sum()) if len(v) else np.nan,'max_count':float(v.max()) if len(v) else np.nan,'pseudobulk_CP10k':float(v.sum()/t.loc[valid,'full_library'].sum()*1e4) if len(v) else np.nan})
su=pd.DataFrame(summ);su.to_csv(O/'tables/desmosome_expression_summary.tsv',sep='\t',index=False)
d[(d['DSC1']>=1)].to_csv(O/'tables/DSC1_RNA_detected_cells.tsv',sep='\t',index=False)
# Resolve absence of a published marker hit versus absence of assay coverage.
hits=[]
for folder in [B/'phase62_F1F5_internal_20261005/tables',B/'phase58_F1F5_round1_20261004/tables']:
 for p in folder.glob('*marker*'):
  if p.is_file() and p.stat().st_size<20000000:
   try:
    t=pd.read_csv(p,sep='\t',low_memory=False);cols=[c for c in t if c.lower() in ['gene','symbol','feature']]
    for c in cols:
     q=t[t[c].eq('DSC1')]
     if len(q):hits.append({'file':str(p),'rows':len(q)})
   except Exception:pass
(O/'resources/DSC1_PRIOR_MARKER_LIST_SEARCH.json').write_text(json.dumps({'scope':'Top-level marker-named tables in phase58/62, not exhaustive full-transcriptome DE','hits':hits},indent=2))
print(su[(su.gene=='DSC1')&(su.level=='class_id')].to_string(index=False),flush=True)
print('DSC1 total',d.groupby('matrix').DSC1.agg(['count','sum','max']).to_dict(),flush=True)
(O/'MARKER_COMPLETE.json').write_text(json.dumps({'success':True,'cells':8174,'genes':genes,'identity':ident,'at':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
