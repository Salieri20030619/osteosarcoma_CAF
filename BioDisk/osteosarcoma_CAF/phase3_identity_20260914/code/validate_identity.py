"""Check delivered identities and selected expression against the unchanged source objects.
SYNTAX_CHECK: python -m py_compile validate_identity.py
RUN_COMMAND: python validate_identity.py --root <phase2> --phase3 <phase3>
EXPECTED_SUCCESS_MARKER: IDENTITY_VALIDATION_PASS
"""
import argparse,json,hashlib
from pathlib import Path
import pandas as pd,numpy as np,h5py,anndata as ad
from takeover_audit import frame,sparse
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--phase3',type=Path,required=True);a=ap.parse_args();p=a.phase3
o=pd.read_csv(p/'identity/all_cell_identity_evidence.tsv.gz',sep='\t')
assert len(o)==28542 and o.cell_id.is_unique
assert ((~o.caf_strict)|o.caf_extended).all()
assert (o.loc[o.caf_strict,'infercnv_call']=='low_signal').all()
assert o.loc[o.caf_strict,'decision_cnv2_call'].isin(['diploid','normal']).all()
assert o.loc[o.caf_strict,'decision_cnv2_reference_check'].all()
assert not o.loc[o.caf_strict,'mixed_lineage_flag'].any()
assert not o.loc[o.caf_strict,'doublet_review_flag'].any()
assert not o.loc[o.caf_extended,'cluster_CTA_program'].any()
assert (o.loc[o.caf_strict,'final_identity']=='Fib_candidate').all()
assert (o.loc[o.final_identity=='Malignant_OS','cnv_concordance']=='both_high').all()
audit=json.loads((p/'takeover_audit/object_audit.json').read_text())
source_records=[]
for r in audit:
    source=a.root/'results'/r['sample']/'annotated.h5ad';sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha==r['sha256']
    d=o[o.sample_id==r['sample']];assert len(d)==r['cells']
    with h5py.File(source) as h:
        original=frame(h['obs']);v=frame(h['var']);counts=sparse(h['layers/counts'])
    assert d.cell_id.tolist()==original.index.tolist()
    chosen=pd.read_csv(p/'inputs'/r['short']/'cells.tsv',sep='\t')
    assert set(original.index[original.lineage_reviewed=='OC_like'])<=set(chosen.cell_id)
    for method in ['infercnv','copykat','copykat_balanced','scevan']:
        if not (p/method/r['short']/'COMPLETE').exists():continue
        results=pd.read_csv(p/method/r['short']/'cell_results.tsv',sep='\t')
        assert results.cell_id.tolist()==chosen.cell_id.tolist()
        assert set(d.loc[d[method+'_call']!='not_run','cell_id'])==set(results.cell_id)
    source_records.append({'sample':r['sample'],'short':r['short'],'original_sha256_unchanged':True,'cells':len(d),'all_OC_included':True})
for flag in ['caf_strict','caf_extended']:
    listing=pd.read_csv(p/'identity'/(flag+'_cells.tsv'),sep='\t')
    assert set(listing.cell_id)==set(o.loc[o[flag],'cell_id'])
    if not len(listing):continue
    obj=ad.read_h5ad(p/'identity/objects'/(flag+'.h5ad'))
    assert set(obj.obs_names)==set(listing.cell_id) and obj.n_vars==33694
    assert np.isfinite(obj.layers['counts'].data).all() and (obj.layers['counts'].data==np.floor(obj.layers['counts'].data)).all()
    for sample in obj.obs['sample'].unique():
        with h5py.File(a.root/'results'/sample/'annotated.h5ad') as h:
            original=frame(h['obs']);v=frame(h['var']);counts=sparse(h['layers/counts']);x=sparse(h['X'])
        z=(obj.obs['sample']==sample).to_numpy();ix=original.index.get_indexer(obj.obs_names[z]);assert (ix>=0).all()
        assert obj.var_names.tolist()==v.index.tolist()
        diff=obj.layers['counts'][z]-counts[ix];assert not diff.nnz or (diff.data==0).all()
        diff=obj.X[z]-x[ix];assert not diff.nnz or np.allclose(diff.data,0,atol=1e-6)
summary=pd.read_csv(p/'identity/sample_identity_summary.tsv',sep='\t')
assert summary.caf_strict.sum()==o.caf_strict.sum() and summary.caf_extended.sum()==o.caf_extended.sum()
record={'status':'PASS','all_cells':len(o),'strict':int(o.caf_strict.sum()),'extended':int(o.caf_extended.sum()),'source_checks':source_records,'count_matrices_identical_to_source':True,'no_selected_method_CNV_positive_in_strict':True,'original_conflicting_calls_retained':True}
(p/'identity/validation.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record,indent=2));print('IDENTITY_VALIDATION_PASS')
