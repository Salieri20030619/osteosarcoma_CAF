"""Re-export counts only, using fixed, already-reviewed cell selection.
SYNTAX_CHECK: python -m py_compile build_count_inputs.py
RUN_COMMAND: python build_count_inputs.py --root <phase2> --phase3 <phase3>
EXPECTED_SUCCESS_MARKER: COUNTS_EXPORTED
"""
import argparse,gzip,json,hashlib
from pathlib import Path
import pandas as pd,numpy as np,h5py
from scipy.io import mmwrite
from takeover_audit import frame,sparse
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--phase3',type=Path,required=True);a=ap.parse_args()
for p in sorted((a.phase3/'inputs').iterdir()):
    if not p.is_dir():continue
    cells=pd.read_csv(p/'cells.tsv',sep='\t');gsm=cells['sample'].iloc[0]
    with h5py.File(a.root/'results'/gsm/'annotated.h5ad') as h:
        obs=frame(h['obs']);var=frame(h['var']);x=sparse(h['layers/counts'])
    ids=obs.index.get_indexer(cells.cell_id);assert (ids>=0).all();m=x[ids]
    assert np.allclose(np.asarray(m.sum(1)).ravel(),cells.n_counts)
    genes=pd.read_csv(p/'genes.tsv',sep='\t');assert genes.ensembl_id.tolist()==var.index.tolist()
    info={'shape_cells_by_genes':list(m.shape),'nnz':m.nnz,'sum_umi':int(m.sum()),'selected_cells':len(cells),'integer_counts':bool((m.data==np.floor(m.data)).all()),'selection_sha256':hashlib.sha256((p/'cells.tsv').read_bytes()).hexdigest()}
    if not (p/'counts.mtx.gz').exists():
        with gzip.open(p/'counts.mtx.gz','wb') as f:mmwrite(f,m.T,field='integer')
    (p/'count_export.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    print('COUNTS_READY',p.name,json.dumps(info),flush=True)
print('COUNTS_EXPORTED')
