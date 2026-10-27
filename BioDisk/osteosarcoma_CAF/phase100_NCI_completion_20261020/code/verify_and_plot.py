import pathlib,json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');p=root/'phase100_NCI_completion_20261020';q=root/'phase95_NCI_expanded_reanalysis_20261016'
build=json.loads((p/'BUILD_COMPLETE.json').read_text(encoding='utf-8'))
reg=pd.read_csv(p/'tables/NCI_expanded_candidate_registry.tsv.gz',sep='\t',encoding='utf-8')
summary=pd.read_csv(p/'tables/final_library_disposition.tsv',sep='\t',encoding='utf-8')
review=pd.read_csv(p/'tables/cluster_review_decisions.tsv',sep='\t',encoding='utf-8')
expr=pd.read_csv(p/'tables/NCI_14gene_library_descriptive.tsv',sep='\t',encoding='utf-8')
manifest=json.loads((p/'objects/COUNT_COLLECTION.json').read_text(encoding='utf-8'))
assert len(summary)==19 and summary['sample'].is_unique and len(reg)==build['expanded_candidates']
assert reg.cell_id.is_unique and reg.final_expanded_admission.all() and not reg.manual_cluster_review_pending.any()
assert reg.mechanical_expanded_candidate.all() and reg.reference_validation_pass.all()
assert reg.full_scDblFinder_class.eq('singlet').all() and reg.native_copykat_call.eq('diploid').all()
assert reg.cluster_review_decision.eq('retain_expanded_identity_unresolved').all()
assert not reg.confirmed_pure_CAF.any() and not reg.independent_patient_verified.any()
assert len(review)==295 and not review.duplicated(['sample','cluster']).any()
assert expr.inferential_P.isna().all() and summary.expanded_candidates.sum()==len(reg)
assert set(reg['sample'])==set(z['sample'] for z in manifest)
for z in manifest:
 f=pathlib.Path(z['path']);assert hashlib.sha256(f.read_bytes()).hexdigest()==z['sha256']
 a=ad.read_h5ad(f,backed='r');assert a.n_obs==z['cells'];assert a.obs.final_expanded_admission.all();a.file.close()
orig=pd.read_csv(p/'resources/EVIDENCE_INPUT_SHA256.tsv',sep='\t',encoding='utf-8')
for z in orig.itertuples():assert hashlib.sha256(pathlib.Path(z.path).read_bytes()).hexdigest()==z.sha256
run=json.loads((q/'COMPUTATION_FINISHED.json').read_text(encoding='utf-8'))
assert run['exit_code']==0 and len(run['samples'])==19
for r in run['samples']:
 assert r['ambient_exit']==0
 if r['status']=='native_completed':assert r['copykat_exit']==0 and r['copykat_marker']
 else:assert r['status']=='insufficient_references'
g=summary.sort_values('sample').copy();y=np.arange(len(g))
fig,ax=plt.subplots(figsize=(11,8.4))
other=g.expanded_candidates-g.author_tumor_conflict-g.author_CAF_agreement
ax.barh(y,g.author_CAF_agreement,color='#268d81',label='Author CAF agreement')
ax.barh(y,other,left=g.author_CAF_agreement,color='#91b7ce',label='Other author host labels')
ax.barh(y,g.author_tumor_conflict,left=g.author_CAF_agreement+other,color='#d98a61',label='Author tumor-label conflict')
for j,r in enumerate(g.itertuples()):
 label=str(r.expanded_candidates) if r.reference_pass else 'Reference not released'
 ax.text(r.expanded_candidates+45,j,label,va='center',fontsize=8)
ax.set_yticks(y,g['sample'].str.replace('NCI__','',regex=False),fontsize=9);ax.invert_yaxis();ax.set_xlim(0,6500);ax.set_xlabel('Expanded candidate nuclei (not independent patients)')
ax.set_title('NCI resource completion: retained candidates and identity conflicts',loc='left',fontsize=13,pad=14)
ax.legend(loc='lower right',frameon=False,fontsize=9);ax.spines[['top','right']].set_visible(False);fig.tight_layout()
fig.savefig(p/'figures/NCI_candidate_disposition.png',dpi=170);plt.close(fig)
genes=json.loads((root/'phase96_all14_literature_20261016/project_registry.json').read_text(encoding='utf-8'));order=[z['gene'] for z in genes]
samples=sorted(reg['sample'].unique())
fig,axes=plt.subplots(1,2,figsize=(13,7.8),sharey=True)
for ax,stratum,title in zip(axes,['expanded','author_CAF_agreement'],['Expanded candidates','Author-CAF-agreeing subset']):
 mat=expr[expr.stratum.eq(stratum)].pivot(index='sample',columns='gene',values='expression').reindex(index=samples,columns=order)
 im=ax.imshow(mat.to_numpy(),aspect='auto',cmap='viridis',vmin=-3.322,vmax=4.5)
 ax.set_xticks(range(len(order)),order,rotation=65,ha='right',fontsize=8)
 ax.set_yticks(range(len(samples)),[s.replace('NCI__','') for s in samples],fontsize=8);ax.set_title(title,fontsize=12)
fig.subplots_adjust(left=.17,right=.9,bottom=.2,top=.88,wspace=.1)
cax=fig.add_axes([.92,.22,.015,.62]);fig.colorbar(im,cax=cax,label='log2(CP10k + 0.1)')
fig.suptitle('Descriptive 14-gene expression by library; no patient-level significance test',fontsize=12,y=.96)
fig.savefig(p/'figures/NCI_14gene_descriptive.png',dpi=170);plt.close(fig)
qa=dict(status='PASS_PROCESSING_SCOPE',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),input_libraries=19,native_completed=18,insufficient_reference_libraries=1,reference_pass_libraries=int(summary.reference_pass.sum()),cluster_dispositions=len(review),expanded_candidates=len(reg),object_hashes_checked=len(manifest),roundtrip_reported_all=all(z['raw_counts_roundtrip'] for z in manifest),source_evidence_hashes_unchanged=19,canonical_admission_flags_reconciled=True,new_inferential_P=False,patient_identity_unverified=True)
(p/'resources/PROCESSING_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa,indent=2),flush=True)
