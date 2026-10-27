"""SYNTAX_CHECK: python -m py_compile quality_scan.py
RUN_COMMAND: python quality_scan.py <phase3>
EXPECTED_SUCCESS_MARKER: QUALITY_SCAN_COMPLETE
"""
from pathlib import Path
import sys,pandas as pd
p=Path(sys.argv[1]);rows=[]
for method in ['infercnv','copykat','copykat_balanced','scevan']:
    normal='low_signal' if method=='infercnv' else ('normal' if method=='scevan' else 'diploid')
    for f in (p/method).glob('*/cell_results.tsv'):
        d=pd.read_csv(f,sep='\t')
        for role in ['ref_immune','ref_endo','heldout','FibLike']:
            mask=d.cnv_role.str.startswith('holdout_') if role=='heldout' else ((d.lineage_phase3=='FibLike') if role=='FibLike' else d.cnv_role==role)
            q=d[mask&(d.call!='not_evaluable')]
            rows.append({'method':method,'sample':f.parent.name,'group':role,'evaluable':len(q),'normal':int((q.call==normal).sum()),'normal_fraction':round(float((q.call==normal).mean()),4)})
out=pd.DataFrame(rows);out.to_csv(p/'review/quality_scan.tsv',sep='\t',index=False);print(out.to_string(index=False));print('QUALITY_SCAN_COMPLETE')
