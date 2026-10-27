from pathlib import Path
import json,re,datetime,sys,importlib.metadata
import pandas as pd,anndata as ad
O=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/geo_identity_20260919')
summary=pd.read_csv(O/'tables/sample_identity_summary.tsv',sep='\t').set_index('sample')
checks=[]
for p in sorted((O/'objects').glob('*.h5ad')):
 a=ad.read_h5ad(p,backed='r')
 sm=p.name.replace('_supported_CAF_candidates.h5ad','')
 assert a.n_obs==int(summary.at[sm,'supported_candidate'])
 assert a.uns['annotation_status']=='joint_lineage_CNV_and_ambient_review_complete; operational_candidates'
 assert a.obs.caf_supported_candidate.all()
 checks.append({'sample':sm,'cells':a.n_obs,'status_metadata':'PASS'})
 a.file.close()
assert sum(x['cells'] for x in checks)==2065
pd.DataFrame(checks).to_csv(O/'tables/object_metadata_acceptance.tsv',sep='\t',index=False)
conv=[]
for p in sorted((O/'logs').glob('GSE*.log')):
 text=p.read_text()
 matches=re.findall(r'Completed iteration: (\d+) \| converge: ([0-9.eE+-]+)',text)
 assert matches and 'DECONTX_SAMPLE_COMPLETE' in text,p
 it,metric=matches[-1];assert int(it)<500 and float(metric)<.001,p
 conv.append({'sample':p.stem,'iterations':int(it),'final_convergence':float(metric),'pass':True})
assert len(conv)==12
pd.DataFrame(conv).to_csv(O/'tables/decontx_convergence.tsv',sep='\t',index=False)
versions={'python':sys.version,**{p:importlib.metadata.version(p) for p in ['anndata','pandas','numpy','scipy','matplotlib']}}
(O/'python_environment.json').write_text(json.dumps(versions,indent=2))
qa={'status':'PASS','reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'figures':['figures/identity_funnel.png','figures/supported_markers.png'],'visual_checks':['readable labels and complete legends','all FibLike cells conserved in stacked counts','lung samples identified','corrected detection-sensitive cells separately shown','fibroblast markers visible; residual ALPL signal not hidden','no inferred biological P values or independent-validation claims'],'input_integrity':'PASS','objects_metadata':'PASS','decontx_convergence':'12/12 PASS'}
(O/'VISUAL_QA_PASS.json').write_text(json.dumps(qa,indent=2))
(O/'README.md').write_text("""# CAF GEO identity review 2026-09-19

Primary report: CAF_GEO_identity_review_20260919.md.
Full all-cell evidence: tables/all_cell_identity_evidence.tsv.gz.
Primary candidates: tables/supported_CAF_candidates.tsv.
Detection-threshold-sensitive candidates: tables/sensitivity_only_CAF_candidates.tsv.
Candidate count overlap: tables/candidate_sensitivity_tiers.tsv.
Candidate expression objects: objects/ (one file per sample; original counts plus separate decontx_counts).

Recorded computation order:
1. code/run.sh -> code/prepare_and_run.py -> code/run_decontx.R
2. code/wait_integrate.py waits for the native completion event -> code/integrate.py
3. code/annotate_sensitivity.py preserves all decision flags and clarifies sensitivity labels
4. code/build_report.py; code/plot_review.py
5. code/acceptance.py and visual inspection; code/finalize.py

The original accepted GEO objects and CNV results remain on sgs10.
Full ambient correction inputs and corrected matrices remain under geo_identity_20260919/inputs and decontx.
The delivery archive excludes these large intermediate matrices.
Original stage paths are embedded in the scripts. Reproduction in a new version directory requires updating the output root.
SHA256SUMS.tsv covers delivered evidence, objects, code and logs.
No cross-cohort patient pooling, state clustering, clinical inference or CNV rerun is included.
""",encoding='utf-8')
print(json.dumps(qa),flush=True)
