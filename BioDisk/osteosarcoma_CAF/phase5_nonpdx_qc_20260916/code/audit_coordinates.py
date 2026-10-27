"""SYNTAX_CHECK: python -m py_compile code/audit_coordinates.py
RUN_COMMAND: python code/audit_coordinates.py
EXPECTED_SUCCESS_MARKER: PHASE5_COORDINATE_AUDIT_COMPLETE
"""
import sys,json,gzip,hashlib
from pathlib import Path
import pandas as pd
r=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p5=r/'phase5_nonpdx_qc_20260916';p3=r/'phase3_identity_20260914'
sys.path.insert(0,str(p3/'code'))
from coordinate_audit import read_gtf
refs={v:read_gtf(p3/'resources'/name) for v,name in [('GRCh38_84','Homo_sapiens.GRCh38.84.gtf.gz'),('GRCh37_75','Homo_sapiens.GRCh37.75.gtf.gz')]}
refs['GRCh38_98']=read_gtf(p5/'resources/Homo_sapiens.GRCh38.98.gtf.gz')
records=json.loads((r/'phase4_download_20260914/matrix_input_resolution_20260916.json').read_text());rows=[];seen=set()
for x in records:
 if x['role']!='features':continue
 p=Path(x['selected']);sample=p.parent.parent.name+'__'+p.parent.name
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt') as f:v=pd.read_csv(f,sep='\t',header=None,dtype=str)
 ids=v[0].str.split('.').str[0];z={'sample':sample,'input_genes':len(ids)}
 for ref,g in refs.items():
  hits=ids.isin(g.ensembl_id);z[ref+'_matched']=int(hits.sum())
  mapping=pd.DataFrame({'ensembl_id':ids,'symbol_input':v[1]}).merge(g.rename(columns={'symbol':'symbol_reference'}),how='left',on='ensembl_id',validate='one_to_one')
  z[ref+'_exact_symbols']=int((mapping.symbol_input==mapping.symbol_reference).sum())
  if ref.startswith('GRCh38'):mapping.to_csv(p5/'resources'/(sample+'_'+ref+'_gene_mapping.tsv.gz'),sep='\t',index=False)
 z['in_GRCh38_not_GRCh37']=int((ids.isin(refs['GRCh38_84'].ensembl_id)&~ids.isin(refs['GRCh37_75'].ensembl_id)).sum())
 rows.append(z)
g=refs['GRCh38_98'];g=g[g.chr.isin([str(i) for i in range(1,23)])].copy();g['chrnum']=g.chr.astype(int);g=g.sort_values(['chrnum','start','end','ensembl_id'])
g[['ensembl_id','chr','start','end']].to_csv(p5/'resources/infercnv_gene_order_GRCh38_98.tsv',sep='\t',header=False,index=False)
pd.DataFrame(rows).to_csv(p5/'resources/coordinate_audit.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False));print('PHASE5_COORDINATE_AUDIT_COMPLETE')
