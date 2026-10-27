"""Exclude five strong mixed-lineage cells from baseline/holdout without dropping them from CNV.
SYNTAX_CHECK: python -m py_compile refine_reference.py
RUN_COMMAND: python refine_reference.py <phase3>
EXPECTED_SUCCESS_MARKER: REFERENCE_REFINEMENT_COMPLETE
"""
import sys,json,hashlib
from pathlib import Path
import pandas as pd,numpy as np
p=Path(sys.argv[1]);records=[]
s=pd.read_csv(p/'review/input_summary.tsv',sep='\t')
for f in sorted((p/'review').glob('*_cell_evidence.tsv.gz')):
    o=pd.read_csv(f,sep='\t');short=o.short.iloc[0];mask=(o.cnv_role.str.startswith('ref_')|o.cnv_role.str.startswith('holdout_'))&o.mixed_lineage_flag
    if not mask.any():continue
    out=p/'inputs'/short
    before=out/'cells.tsv';backup=out/'cells_before_refinement.tsv'
    if not backup.exists():backup.write_bytes(before.read_bytes())
    for r in o.loc[mask].itertuples():records.append({'short':short,'cell_id':r.cell_id,'tool_id':r.tool_id,'old_role':r.cnv_role,'new_role':'target','reason':'strong_mixed_lineage_evidence','det_Osteo':r.det3_Osteo,'det_Endo':r.det3_Endo,'det_Myeloid':r.det3_Myeloid})
    o.loc[mask,'cnv_role']='target';o.loc[mask,'cnv_target_expanded']=True
    o.to_csv(f,sep='\t',index=False)
    sel=o.cnv_role!='not_selected';n=o.loc[sel].copy()
    old=pd.read_csv(before,sep='\t');assert n.tool_id.tolist()==old.tool_id.tolist()
    n.to_csv(before,sep='\t',index=False)
    anns=pd.DataFrame({'id':n.tool_id,'group':np.where(n.cnv_role.str.startswith('ref_'),n.cnv_role,'C'+n.cluster.astype(str)+'_'+n.lineage_phase3)})
    anns.to_csv(out/'annotations.tsv',sep='\t',header=False,index=False)
    ix=s.short==short
    s.loc[ix,'target_expanded_total']=o.cnv_target_expanded.sum();s.loc[ix,'target_run']=(o.cnv_role=='target').sum()
    s.loc[ix,'ref_immune']=(o.cnv_role=='ref_immune').sum();s.loc[ix,'ref_endo']=(o.cnv_role=='ref_endo').sum();s.loc[ix,'holdout']=o.cnv_role.str.startswith('holdout_').sum()
    meta=out/'count_export.json'
    if meta.exists():
        d=json.loads(meta.read_text());d['selection_sha256_before_reference_refinement']=d['selection_sha256'];d['selection_sha256']=hashlib.sha256(before.read_bytes()).hexdigest();d['count_matrix_unchanged']=True;meta.write_text(json.dumps(d,indent=2))
s.to_csv(p/'review/input_summary.tsv',sep='\t',index=False)
if records:pd.DataFrame(records).to_csv(p/'review/reference_exclusions.tsv',sep='\t',index=False)
print('REFERENCE_REFINEMENT_COMPLETE',len(records))
print(s.to_string(index=False))
