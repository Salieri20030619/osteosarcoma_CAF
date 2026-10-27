from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,datetime,fcntl,sys
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');O=P/'phase91_lung_patient_validation_20261015';sm='GSE234187__GSM7453693'
review_lock=(O/'evidence_review.lock').open('a');fcntl.flock(review_lock,fcntl.LOCK_EX)
if (O/'EVIDENCE_REVIEW.json').exists() and '--force' not in sys.argv:
 print((O/'EVIDENCE_REVIEW.json').read_text(encoding='utf-8'));raise SystemExit(0)
done=json.loads((O/'JOB_FINISHED.json').read_text(encoding='utf-8'));assert done['exit_code']==0
assert all(r['exit_code']==0 and r['completion_marker'] for r in done['jobs'])
assert json.loads((O/'resources/FINGERPRINT_EXIT.json').read_text(encoding='utf-8'))['exit_code']==0
fp=json.loads((O/'resources/FINGERPRINT_COMPLETE.json').read_text(encoding='utf-8'));assert fp['comparisons']==116
c=pd.read_csv(O/'inputs'/sm/'cells.tsv',sep='\t').set_index('cell_id')
cnv=pd.read_csv(O/'cnv/infercnv'/sm/'cell_results.tsv',sep='\t').set_index('cell_id')
assert cnv.index.is_unique and set(c.index)==set(cnv.index)
c['native_HMM_call']=cnv['call'];assert c.native_HMM_call.isin(['alteration_called','no_alteration_called','not_returned_by_software']).all()
review=pd.read_csv(O/'tables/cluster_adjudication_frozen.tsv',sep='\t').set_index('cluster')
c['cluster_decision']=c.cluster.map(review.decision);c['cluster_reason']=c.cluster.map(review.reason)
assert c.cluster_decision.notna().all()
cal=[]
for label,g in c[c.heldout].groupby('reference_role'):
 ev=g.native_HMM_call.isin(['no_alteration_called','alteration_called']);nn=int(ev.sum());normal=int(g.native_HMM_call.eq('no_alteration_called').sum())
 cal.append({'group':label,'heldout_cells':len(g),'evaluable':nn,'no_alteration':normal,'normal_fraction':normal/nn if nn else None,'pass':nn>=10 and normal/nn>=.9 if nn else False})
reference_pass=bool(cal) and all(r['pass'] for r in cal)
pd.DataFrame(cal).to_csv(O/'tables/native_reference_validation.tsv',sep='\t',index=False,encoding='utf-8')
a=ad.read_h5ad(O/'objects/full_called_counts.h5ad');assert a.obs_names.tolist()==c.index.tolist()
y=mmread(O/'ambient/corrected_counts.mtx.gz').T.tocsr();assert y.shape==a.shape and np.isfinite(y.data).all() and (y.data>=0).all()
amb=pd.read_csv(O/'ambient/contamination.tsv',sep='\t').set_index('cell_id');assert set(amb.index)==set(c.index)
c['decontX_contamination']=amb.contamination
panels=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text(encoding='utf-8'))
symbols=a.var.symbol.astype(str).to_numpy();allgenes=sorted(set(sum(panels.values(),[])))
gm={g:np.flatnonzero(symbols==g) for g in allgenes};raw=np.column_stack([np.asarray(y[:,gm[g]].sum(1)).ravel() for g in allgenes])
den=np.asarray(y.sum(1)).ravel();log=np.log1p(raw/np.maximum(den,1)[:,None]*1e4)
scores={k:log[:,[allgenes.index(g) for g in gs if len(gm[g])]].mean(1) for k,gs in panels.items()}
v=np.column_stack(list(scores.values()));mx=v.max(1);winner=np.array(list(panels),dtype=object)[v.argmax(1)];winner[(mx<=0)|(np.sum(v==mx[:,None],axis=1)!=1)]='Unresolved'
c['decontX_lineage_winner']=winner;c['decontX_detected_Fib_genes_ge1']=(raw[:,[allgenes.index(g) for g in panels['Fibroblast']]]>=1).sum(1)
for k,z in scores.items():c['decontX_phase27_'+k]=z
eligible=c.provisional_lineage.eq('Fibroblast')&c.full_scDblFinder_class.eq('singlet')&c.cluster_decision.eq('CAF_compatible_pending_native_validation')&c.native_HMM_call.eq('no_alteration_called')
c['native_reviewed_candidate']=eligible
c['reference_validated_candidate']=eligible&reference_pass
c['ambient_supported_candidate']=c.reference_validated_candidate&c.decontX_lineage_winner.eq('Fibroblast')&c.decontX_detected_Fib_genes_ge1.ge(2)
c['formal_admission']=False
c.to_csv(O/'tables/all_cell_validation.tsv.gz',sep='\t',index_label='cell_id',encoding='utf-8')
capacity=c.groupby('cluster').agg(cells=('sample','size'),software_doublets=('full_scDblFinder_class',lambda s:int(s.eq('doublet').sum())),raw_Fibroblast=('provisional_lineage',lambda s:int(s.eq('Fibroblast').sum())),native_no_alteration=('native_HMM_call',lambda s:int(s.eq('no_alteration_called').sum())),reviewed_native_candidates=('native_reviewed_candidate','sum'),reference_validated_candidates=('reference_validated_candidate','sum'),ambient_supported_candidates=('ambient_supported_candidate','sum'),median_contamination=('decontX_contamination','median'))
capacity['review_decision']=review.decision;capacity.to_csv(O/'tables/cluster_validation_summary.tsv',sep='\t',encoding='utf-8')
sel=c.reference_validated_candidate.to_numpy();candidate=a[sel].copy();candidate.obs=c.loc[sel].copy();candidate.uns['source_patient_count']=1;candidate.uns['candidate_patient_count']=int(candidate.n_obs>0);candidate.uns['scientific_admission']='review_pending';candidate.write_h5ad(O/'objects/reference_validated_candidate_counts.h5ad',compression='gzip')
check=ad.read_h5ad(O/'objects/reference_validated_candidate_counts.h5ad',backed='r');assert check.obs_names.tolist()==candidate.obs_names.tolist() and check.var_names.tolist()==a.var_names.tolist();check.file.close()
rawcounts=np.asarray(a.X.sum(1)).ravel();assert np.array_equal(rawcounts,c.n_counts.to_numpy())
assert not c.loc[sel,'full_scDblFinder_class'].eq('doublet').any()
assert c.loc[sel,'native_HMM_call'].eq('no_alteration_called').all()
steps=[('Full called',len(c)),('Software singlets',int(c.full_scDblFinder_class.eq('singlet').sum())),('Fibroblast winner',int((c.full_scDblFinder_class.eq('singlet')&c.provisional_lineage.eq('Fibroblast')).sum())),('Reviewed groups',int((c.full_scDblFinder_class.eq('singlet')&c.provisional_lineage.eq('Fibroblast')&c.cluster.isin([5,7,11])).sum())),('Native no-alteration',int(eligible.sum())),('Reference validated',int(sel.sum())),('Ambient supported',int(c.ambient_supported_candidate.sum()))]
assert all(steps[i][1]>=steps[i+1][1] for i in range(len(steps)-1))
pd.DataFrame(steps,columns=['stage','cells']).to_csv(O/'tables/candidate_waterfall.tsv',sep='\t',index=False,encoding='utf-8')
(O/'figures').mkdir(exist_ok=True)
fig,axes=plt.subplots(1,3,figsize=(17,5.6),layout='constrained')
axes[0].barh([s[0] for s in steps][::-1],[s[1] for s in steps][::-1],color='#2a788e')
for i,(_,n) in enumerate(steps[::-1]):axes[0].text(n+max(len(c)*.015,1),i,str(n),va='center',fontsize=9)
axes[0].set_xlim(0,len(c)*1.25);axes[0].set_xlabel('Cells');axes[0].set_title('GSE234187: one lung-metastasis patient')
axes[1].scatter(c.UMAP1,c.UMAP2,s=2,c='#d5d5d5',rasterized=True)
for cl,color in zip([5,7,11],['#21918c','#440154','#f7b32b']):
 q=c[c.reference_validated_candidate&c.cluster.eq(cl)];axes[1].scatter(q.UMAP1,q.UMAP2,s=6,c=color,label=f'Cluster {cl}: {len(q)}',rasterized=True)
axes[1].set_title('Native/reference-validated candidates');axes[1].set_xlabel('UMAP1');axes[1].set_ylabel('UMAP2');axes[1].legend(frameon=False,fontsize=9)
sub=c[c.native_reviewed_candidate];dots=axes[2].scatter(sub.phase27_Fibroblast,sub.decontX_phase27_Fibroblast,c=sub.decontX_contamination,cmap='viridis',s=8,vmin=0,vmax=1,rasterized=True)
fig.colorbar(dots,ax=axes[2],label='decontX contamination',shrink=.7)
mx=max(sub.phase27_Fibroblast.max(),sub.decontX_phase27_Fibroblast.max(),1) if len(sub) else 1
axes[2].plot([0,mx],[0,mx],color='grey',ls='--',lw=1);axes[2].set_xlabel('Raw fibroblast mean log1p(CP10k)');axes[2].set_ylabel('decontX fibroblast mean log1p(CP10k)');axes[2].set_title('Ambient RNA sensitivity (descriptive)')
fig.savefig(O/'figures/patient_validation.png',dpi=180);plt.close(fig)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'evidence_review_ready_not_admitted','source':'GSE234187','patient':'GSE234187_patient1','anatomical_site':'Lung','source_patient_N':1,'cells':len(c),'native_HMM_calls':c.native_HMM_call.value_counts().to_dict(),'reference_validation':cal,'reference_pass':reference_pass,'reviewed_native_candidates':int(eligible.sum()),'reference_validated_candidates':int(sel.sum()),'ambient_supported_candidates':int(c.ambient_supported_candidate.sum()),'fingerprint_comparisons':fp['comparisons'],'exact_shared_count_vectors':fp['matched_vectors'],'raw_counts_unchanged':True,'candidate_object_readback_passed':True,'new_inferential_tests':False,'new_formal_patients':0,'formal_admission_pending':True}
(O/'EVIDENCE_REVIEW.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2));print(capacity.to_string())
