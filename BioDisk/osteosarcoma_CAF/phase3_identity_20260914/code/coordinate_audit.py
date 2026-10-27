"""SYNTAX_CHECK: python -m py_compile coordinate_audit.py
RUN_COMMAND: python coordinate_audit.py --root <phase2> --phase3 <phase3>
EXPECTED_SUCCESS_MARKER: COORDINATES_VERIFIED
"""
import argparse,gzip,re,json,hashlib
from pathlib import Path
import h5py,pandas as pd,numpy as np
from takeover_audit import frame
def read_gtf(p):
    rows=[]
    with gzip.open(p,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            a=line.rstrip('\n').split('\t')
            if a[2]!='gene':continue
            info=dict(re.findall(r'(\w+) "([^"]+)"',a[8]))
            rows.append({'ensembl_id':info['gene_id'].split('.')[0],'chr':a[0],'start':int(a[3]),'end':int(a[4]),'symbol':info.get('gene_name',''),'biotype':info.get('gene_biotype','')})
    return pd.DataFrame(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--phase3',type=Path,required=True);a=ap.parse_args()
    r=a.phase3/'resources'
    g38=read_gtf(r/'Homo_sapiens.GRCh38.84.gtf.gz');g37=read_gtf(r/'Homo_sapiens.GRCh37.75.gtf.gz')
    summary=[];first=None
    for p in sorted((a.root/'results').glob('GSM*/annotated.h5ad')):
        with h5py.File(p) as h:v=frame(h['var'])
        ids=v.index.astype(str);m38=v[['symbol']].rename_axis('ensembl_id').reset_index().merge(g38,on='ensembl_id',suffixes=('_input','_gtf'),how='left')
        n38=int(m38.chr.notna().sum());n37=int(ids.isin(g37.ensembl_id).sum())
        exact_symbols=int((m38.symbol_input==m38.symbol_gtf).sum())
        newer=ids[~ids.isin(g37.ensembl_id)&ids.isin(g38.ensembl_id)]
        summary.append({'sample':p.parent.name,'input_genes':len(v),'GRCh38_84_gene_matches':n38,'GRCh38_84_exact_symbols':exact_symbols,'GRCh37_75_gene_matches':n37,'genes_in_38_not_37':len(newer),'geo_reported_assembly':'hg19'})
        m38.to_csv(r/(p.parent.name+'_gene_coordinate_mapping.tsv.gz'),sep='\t',index=False)
        if first is None:first=m38
        else:assert m38.ensembl_id.tolist()==first.ensembl_id.tolist()
        assert n38/len(v)>.99,(p,n38)
        assert exact_symbols/len(v)>.99,(p,exact_symbols)
        assert len(newer)>=50,(p,len(newer))
    autos=[str(i) for i in range(1,23)]
    use=first[first.chr.isin(autos)].copy();assert not use.ensembl_id.duplicated().any()
    use['chrnum']=use.chr.astype(int);use=use.sort_values(['chrnum','start','end','ensembl_id'])
    use[['ensembl_id','chr','start','end']].to_csv(r/'infercnv_gene_order_GRCh38_84.tsv',sep='\t',header=False,index=False)
    pd.DataFrame(summary).to_csv(r/'coordinate_audit.tsv',sep='\t',index=False)
    (r/'coordinate_provenance.json').write_text(json.dumps({'inference':'RNA feature annotation matches Ensembl84 GRCh38; GEO hg19 preserved as conflicting submitted metadata. Original alignment command/reference not recovered. CNV gene positions mapped explicitly by stable Ensembl ID to GRCh38; no FASTQ realignment claimed.','coordinate_source':'https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz','comparison_source':'https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz','autosomal_genes':len(use),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in r.glob('Homo_sapiens*.gtf.gz')}},indent=2),encoding='utf-8')
    print(pd.DataFrame(summary).to_string(index=False));print('COORDINATES_VERIFIED')
if __name__=='__main__':main()
