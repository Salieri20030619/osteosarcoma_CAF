import pathlib,json,gzip,hashlib
import numpy as np,pandas as pd
p=pathlib.Path(__file__).resolve().parents[1]
genes=json.loads((p/'resources/PROTOCOL_FROZEN.json').read_text())['genes']
aliases=json.loads((p/'resources/GENE_ALIAS_MAP.json').read_text())['aliases']
reverse={alias:g for g,aa in aliases.items() for alias in aa}
matched_labels={}
m=pd.read_csv(p/'inputs/lung_meta_verified.tsv',sep='\t',encoding='utf-8',dtype=str).set_index('cellId')
m['cell_type_base']=m.cell_type.str.replace(r'_P\d+$','',regex=True)
print(m.groupby(['cell_type_base','patient']).size().to_string(),flush=True)
sel=m[m.cell_type_base.str.contains('fibroblast',case=False,regex=False)&m.tissue.eq('lung')].copy()
assert not sel.index.duplicated().any()
sel.to_csv(p/'tables/lung_selected_cell_metadata.tsv',sep='\t',encoding='utf-8')
groups={}
for spec in ['all_fibroblasts','alveolar_only','adventitial_only']:
 q=sel if spec=='all_fibroblasts' else sel[sel.cell_type_base.eq('Alveolar Fibroblast' if spec=='alveolar_only' else 'Adventitial Fibroblast')]
 for donor,g in q.groupby('patient'):groups[(spec,donor)]=g.index
with gzip.open(p/'inputs/lung_counts_verified.tsv.gz','rt',encoding='utf-8') as f:
 header=f.readline().rstrip('\n\r').split('\t');cells=pd.Index(header[1:])
 assert cells.is_unique
 pos={k:cells.get_indexer(v) for k,v in groups.items()}
 assert all((v>=0).all() for v in pos.values())
 totals={k:0. for k in groups};counts={};nrows=0;seen=set()
 for line in f:
  label,values=line.rstrip('\n\r').split('\t',1)
  a=np.fromstring(values,sep='\t')
  assert len(a)==len(cells) and np.isfinite(a).all() and (a>=0).all() and (a==np.floor(a)).all()
  assert label not in seen;seen.add(label)
  sm={k:float(a[ix].sum()) for k,ix in pos.items()}
  for k,v in sm.items():totals[k]+=v
  if label in reverse:
   gene=reverse[label]
   assert gene not in counts,('ambiguous gene rows',gene,label)
   counts[gene]=sm;matched_labels[gene]=label
  nrows+=1
  if nrows%5000==0:print('genes parsed',nrows,flush=True)
libcheck=[]
for key,indices in groups.items():
 expected=float(sel.loc[indices,'nUMI'].astype(float).sum())
 libcheck.append(dict(specification=key[0],donor=key[1],matrix_UMIs=totals[key],metadata_UMIs=expected,difference=totals[key]-expected))
pd.DataFrame(libcheck).to_csv(p/'tables/lung_library_roundtrip.tsv',sep='\t',encoding='utf-8',index=False)
assert all(z['difference']==0 for z in libcheck),'lung matrix and original nUMI denominator disagree'
rows=[]
for (spec,donor),ids in groups.items():
 for gene in genes:
  raw=counts[gene][(spec,donor)] if gene in counts else np.nan
  rows.append(dict(reference='lung',specification=spec,donor=donor,gene=gene,raw_count=raw,library_UMIs=totals[(spec,donor)],cells=len(ids),measured=gene in counts,expression=np.log2(raw/totals[(spec,donor)]*1e4+.1) if gene in counts else np.nan))
pd.DataFrame([dict(gene=g,reference_row=matched_labels.get(g,''),measured=g in counts) for g in genes]).to_csv(p/'tables/lung_gene_coverage.tsv',sep='\t',encoding='utf-8',index=False)
pd.DataFrame(rows).to_csv(p/'tables/lung_reference_pseudobulk.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
(p/'resources/LUNG_EXTRACTION_QA.json').write_text(json.dumps(dict(matrix_genes=nrows,matrix_cells=len(cells),all_integer_nonnegative=True,missing_genes=[g for g in genes if g not in counts],selected_cells=len(sel),donors=int(sel.patient.nunique()),gene_matching='published symbol or documented stable-ID alias; no missing-as-zero'),indent=2),encoding='utf-8')
print('LUNG_DONE',flush=True)
