"""SYNTAX_CHECK: python -m py_compile summarize_cnv.py
RUN_COMMAND: python summarize_cnv.py <phase3>
EXPECTED_SUCCESS_MARKER: CNV_STATUS_SUMMARIZED
"""
import sys,json
from pathlib import Path
import pandas as pd
p=Path(sys.argv[1]);out=[]
for method in ['infercnv','copykat']:
    for f in sorted((p/method).glob('*/cell_results.tsv')):
        d=pd.read_csv(f,sep='\t');held=d.cnv_role.str.startswith('holdout_');refs=d.cnv_role.str.startswith('ref_')
        row={'method':method,'sample':f.parent.name,'calls':d.call.value_counts().to_dict(),'heldout_calls':d.loc[held,'call'].value_counts().to_dict(),'reference_calls':d.loc[refs,'call'].value_counts().to_dict(),'by_lineage':pd.crosstab(d.lineage_phase3,d.call).to_dict(),'median_signal':d.groupby('lineage_phase3').signal_strength.median().to_dict()}
        out.append(row)
cross=[]
for f in sorted((p/'infercnv').glob('*/cell_results.tsv')):
    g=p/'copykat'/f.parent.name/'cell_results.tsv'
    if g.exists():
        d=pd.read_csv(f,sep='\t');e=pd.read_csv(g,sep='\t');m=d.merge(e[['tool_id','call']],on='tool_id',suffixes=('_ic','_ck'))
        cross.append({'sample':f.parent.name,'cross':pd.crosstab(m.call_ic,m.call_ck).to_dict()})
result={'tool_results':out,'cross':cross}
(p/'review/cnv_status_summary.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2));print('CNV_STATUS_SUMMARIZED')
