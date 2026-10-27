from pathlib import Path
import pandas as pd,json,datetime
B=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=B/'phase64_response_followup_20261005'
p=O/'code/analyze_response.py';s=p.read_text()
old="approved=set(h.symbol);aliases={}"
new="""approved=set(h.symbol);aliases={}
refmap={};entrezmap={}
for _,r in h.iterrows():
 for a in r.get('refseq_accession','').split('|'):
  if a:refmap.setdefault(a.split('.')[0],set()).add(r.symbol)
 for a in r.get('entrez_id','').split('|'):
  if a:entrezmap.setdefault(a,set()).add(r.symbol)
def map_probe(r):
 es=set()
 for a in re.findall(r'\\d+',r.get('Entrez_Gene_ID','')):es.update(entrezmap.get(a,set()))
 rs=set()
 for field in ['RefSeq_ID','Source_Reference_ID','GB_ACC']:
  for a in re.findall(r'[NX][MR]_\\d+',r.get(field,'')):rs.update(refmap.get(a,set()))
 if es and rs and es!=rs:return None,'identifier_conflict'
 if len(es)==1:return next(iter(es)),'Entrez'
 if len(es)>1:return None,'multi_Entrez'
 if len(rs)==1:return next(iter(rs)),'RefSeq'
 if len(rs)>1:return None,'multi_RefSeq'
 val=canon(r['Symbol'].strip())
 return val,'symbol_fallback' if val else 'unresolved_symbol'
"""
assert old in s;s=s.replace(old,new)
s=s.replace("ann['canonical_symbol']=ann.Symbol.map(canon)","ann[['canonical_symbol','mapping_method']]=ann.apply(map_probe,axis=1,result_type='expand')")
s=s.replace("mapping=ann[['ID','Symbol','canonical_symbol']].copy()","mapping=ann[['ID','Symbol','canonical_symbol','mapping_method']].copy()")
p.write_text(s)
d={'at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'Annotation QA revealed older GPL10295 PRG2 and AGC1 labels refer to distinct genes for different RefSeq accessions. Symbol-only mapping incorrectly assigns PXDN NM_024888 to current PRG2 and drops ambiguous AGC1.','rule':'HGNC Entrez and version-stripped RefSeq identifiers first; conflicting unique identifiers excluded; fallback to deterministic current/unique historical symbol. No expression/outcome used for mapping. Same rule applied to both assays, all probes.','status':'Supersedes exploratory interim21257 symbol-only results; final full analysis recalculated and independently verified.','primary_gene_set':'Unchanged fixed six genes; no outcome-based signature changes.'}
(O/'resources/ANNOTATION_QA_AMENDMENT.json').write_text(json.dumps(d,indent=2))
print('AMENDMENT_SAVED',flush=True)
