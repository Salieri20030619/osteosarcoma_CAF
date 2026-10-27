"""SYNTAX_CHECK: python -m py_compile audit_cta_context.py
RUN_COMMAND: python audit_cta_context.py <phase2> <phase3>
EXPECTED_SUCCESS_MARKER: CTA_CONTEXT_AUDITED
"""
from pathlib import Path
import sys,pandas as pd
root=Path(sys.argv[1]);p=Path(sys.argv[2]);rows=[]
genes={'MAGEA3','MAGEA10','MAGEA12','SSX1','CTAG2','PAGE2B','PAGE5','CSAG1','FMR1NB'}
for f in (p/'identity').glob('*_identity.tsv.gz'):
    o=pd.read_csv(f,sep='\t');sample=o['sample'].iloc[0]
    m=pd.read_csv(root/'results'/sample/'cluster_markers.tsv',sep='\t')
    for c,d in m.groupby('group'):
        hits=genes&set(d.symbol)
        if len(hits)<2:continue
        z=o[o.cluster==c]
        rows.append({'sample':sample,'short':o.short.iloc[0],'cluster':c,'CTA_markers':';'.join(sorted(hits)),'cells':len(z),'Fib_candidate':int((z.final_identity=='Fib_candidate').sum()),'strict':int(z.caf_strict.sum()),'extended':int(z.caf_extended.sum())})
out=pd.DataFrame(rows);out.to_csv(p/'review/CTA_cluster_context_audit.tsv',sep='\t',index=False);print(out.to_string(index=False));print('CTA_CONTEXT_AUDITED')
