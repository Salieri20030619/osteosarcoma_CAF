"""SYNTAX_CHECK: python -m py_compile assess_balanced.py
RUN_COMMAND: python assess_balanced.py <phase3>
EXPECTED_SUCCESS_MARKER: BALANCED_REFERENCE_ASSESSED
"""
import sys,json
from pathlib import Path
import pandas as pd
p=Path(sys.argv[1]);out=[]
for s in ['S0058','SC072']:
    f=p/'copykat_balanced'/s/'cell_results.tsv'
    if not f.exists():continue
    d=pd.read_csv(f,sep='\t')
    out.append({'short':s,'by_lineage':pd.crosstab(d.lineage_phase3,d.call).to_dict(),'by_role':pd.crosstab(d.cnv_role,d.call).to_dict()})
print(json.dumps(out,indent=2));print('BALANCED_REFERENCE_ASSESSED')
