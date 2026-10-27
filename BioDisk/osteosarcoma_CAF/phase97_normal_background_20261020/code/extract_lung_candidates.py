import pathlib,json,gzip
import numpy as np,pandas as pd
p=pathlib.Path(__file__).resolve().parents[1]
old=json.loads((p/'resources/LUNG_EXTRACTION_QA.json').read_text())
(p/'resources/LUNG_MATRIX_VALIDATION.json').write_text(json.dumps(old,indent=2),encoding='utf-8')
gmap=json.loads((p/'resources/GENE_ALIAS_MAP.json').read_text());genes=list(gmap['aliases'])
reverse={a:g for g,aa in gmap['aliases'].items() for a in aa}
reverse.update({v:g for g,v in gmap['stable_ids'].items()})
sel=pd.read_csv(p/'tables/lung_selected_cell_metadata.tsv',sep='\t',encoding='utf-8',dtype=str).set_index('cellId')
libs=pd.read_csv(p/'tables/lung_library_roundtrip.tsv',sep='\t',encoding='utf-8',dtype={'donor':str})
assert libs.difference.eq(0).all()
groups={}
for spec in ['all_fibroblasts','alveolar_only','adventitial_only']:
 q=sel if spec=='all_fibroblasts' else sel[sel.cell_type_base.eq('Alveolar Fibroblast' if spec=='alveolar_only' else 'Adventitial Fibroblast')]
 for donor,z in q.groupby('patient'):groups[(spec,donor)]=z.index
totals={(r.specification,r.donor):r.matrix_UMIs for r in libs.itertuples()}
counts={};labels={};nrows=0
with gzip.open(p/'inputs/lung_counts_verified.tsv.gz','rt',encoding='utf-8') as f:
 cells=pd.Index(f.readline().rstrip('\r\n').split('\t')[1:]);pos={k:cells.get_indexer(v) for k,v in groups.items()}
 assert all((v>=0).all() for v in pos.values())
 for line in f:
  label,values=line.rstrip('\n\r').split('\t',1);nrows+=1
  matching={reverse[t] for t in label.split('|') if t in reverse}
  if not matching:continue
  assert len(matching)==1,(label,matching)
  gene=next(iter(matching));assert gene not in counts,('multiple rows for gene',gene,label)
  a=np.fromstring(values,sep='\t');assert len(a)==len(cells) and (a>=0).all() and np.isfinite(a).all() and (a==np.floor(a)).all()
  counts[gene]={k:float(a[ix].sum()) for k,ix in pos.items()};labels[gene]=label
assert nrows==old['matrix_genes']
rows=[]
for (spec,donor),ids in groups.items():
 for gene in genes:
  raw=counts[gene][(spec,donor)] if gene in counts else np.nan
  denom=totals[(spec,donor)];assert np.isnan(raw) or raw<=denom
  rows.append(dict(reference='lung',specification=spec,donor=donor,gene=gene,raw_count=raw,library_UMIs=denom,cells=len(ids),measured=gene in counts,expression=np.log2(raw/denom*1e4+.1) if gene in counts else np.nan))
pd.DataFrame(rows).to_csv(p/'tables/lung_reference_pseudobulk.tsv',sep='\t',encoding='utf-8',index=False,na_rep='NA')
pd.DataFrame([dict(gene=g,reference_row=labels.get(g,''),measured=g in counts) for g in genes]).to_csv(p/'tables/lung_gene_coverage.tsv',sep='\t',encoding='utf-8',index=False)
qa=dict(matrix_genes=nrows,matrix_cells=len(cells),all_integer_nonnegative=old['all_integer_nonnegative'],library_roundtrip=True,missing_genes=[g for g in genes if g not in counts],selected_cells=len(sel),donors=int(sel.patient.nunique()),gene_matching='UCSC ID|symbol parsed; exact stable ID or documented symbol/alias; full original labels retained')
(p/'resources/LUNG_EXTRACTION_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa,indent=2),flush=True)
