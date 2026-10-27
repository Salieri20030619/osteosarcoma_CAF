from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[k]='2'
import json, hashlib, datetime, gc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF')
Q=P/'phase88_expansion_identity_20261014'
O=P/'audit_expanded_admission_20261015'
def dump(name,obj):
    (O/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(f):
    h=hashlib.sha256()
    with f.open('rb') as hnd:
        for b in iter(lambda:hnd.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def truth(s):
    return s.fillna('').astype(str).str.lower().isin(['true','1','1.0'])
def save(df,name):df.to_csv(O/'tables'/name,sep='\t',index=False,encoding='utf-8')
panels=json.loads((P/'phase27_identity_review_20260926/resources/lineage_panels.json').read_text(encoding='utf-8'))
keys=list(panels)
genes=sorted(set(sum(panels.values(),[])))
manifest=pd.read_csv(Q/'preprocess/manifest.tsv',sep='\t',keep_default_na=False)
baseline=[P/'CAF_FUNCTIONAL_BONE_LUNG_CURRENT.json',P/'CAF_ALL14_RESEARCH_CURRENT.json',P/'CAF_EXPANSION_IDENTITY_CURRENT.json']
tracked=baseline+[Q/'PROTOCOL.json',Q/'code/finalize.py',Q/'code/run_cnv.R',Q/'tables/all_cell_identity_evidence.tsv.gz']
before={str(f):sha(f) for f in tracked}
dump('resources/REVIEW_PROTOCOL.json',{
 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'scope':'Read-only expansion feasibility review; temporary outputs isolated; no baseline admission or new inferential tests',
 'current_baseline':'Phase27 competitive lineage/native CNV -> Phase28 cluster review -> Phase30 expanded8174, with Phase66 doublet flags retained',
 'lineage':'Exact Phase27 panels; per-gene log1p(full-library CP10k), mean over measured panel genes; unique positive largest panel',
 'native_entry_audit':'Fibroblast winner AND native CopyKAT diploid AND no inferCNV elevated signal AND no explicit any-method doublet positive',
 'cnv_limit':'Phase88 HMM=False; low_signal is a custom descriptive call, not native no-alteration. Candidate tally is a compatible starting subset, not exact Phase30 recreation.',
 'reference':'Reference-failed native-entry cells are retained in audit counts but separated from calibration-supported candidates; no call changed or failure erased',
 'ambient':'Corrected Fib support recorded as evidence, not a new mandatory expanded-layer entry filter',
 'upstream_limit':'Phase87/88 retained cells only; earlier full-provider-called entry is not reconstructed; cluster review remains required',
 'sample_size':'Within-source donor labels; BS4/BS14 one patient, EXT_OS_1_L1/L2 one donor; no cross-source pooled independent N',
 'old_labels':'Original sample/cluster IDs and all input evidence unchanged'})
allrows=[];summaries=[];clusters=[];dots=[];coverage=[];source_stat=[]
for row in manifest.to_dict('records'):
    sm=row['sample']
    f=Q/'preprocess/accepted'/sm/'analysis_ready.h5ad'
    source_stat.append({'path':str(f),'size':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns})
    a=ad.read_h5ad(f)
    x=a.layers['counts'].tocsr()
    ev=pd.read_csv(Q/'tables'/f'{sm}_identity.tsv.gz',sep='\t',low_memory=False,dtype={'cluster':str}).set_index('cell_id',drop=False)
    assert ev.index.is_unique and ev.index.tolist()==a.obs_names.tolist()
    den=np.asarray(x.sum(1)).ravel()
    assert np.allclose(den,ev.n_counts)
    syms=a.var.symbol.astype(str).to_numpy()
    mapping={g:[] for g in genes}
    for j,sy in enumerate(syms):
        for g in sy.split(';'):
            if g in mapping:mapping[g].append(j)
    raw=np.column_stack([np.asarray(x[:,mapping[g]].sum(1)).ravel() if mapping[g] else np.zeros(len(a)) for g in genes]).astype(np.float32)
    log=np.log1p(raw/np.maximum(den,1)[:,None]*1e4)
    score=[]
    for k,gs in panels.items():
        ix=[genes.index(g) for g in gs if mapping[g]]
        score.append(log[:,ix].mean(1) if ix else np.full(len(a),np.nan))
        coverage.append({'sample':sm,'panel':k,'measured':len(ix),'requested':len(gs),'missing':';'.join(g for g in gs if not mapping[g])})
    v=np.column_stack(score);mx=np.nanmax(v,axis=1)
    label=np.array(keys,dtype=object)[np.nanargmax(v,axis=1)]
    label[(mx<=0)|(np.sum(v==mx[:,None],axis=1)!=1)]='Unresolved'
    ev['phase27_panel_winner']=label
    for j,k in enumerate(keys):ev['phase27_'+k]=v[:,j]
    db=truth(ev.doublet_positive_any)
    for k in ['doublet_exclude','author_doublet_positive']:
        if k in ev:db |= truth(ev[k])
    if 'scDblFinder_class' in ev:db|=ev.scDblFinder_class.eq('doublet')
    if 'scrublet' in ev:db|=ev.scrublet.eq('doublet')
    ev['audit_doublet_positive']=db
    fib=ev.phase27_panel_winner.eq('Fibroblast')
    low=(ev.infer_call.eq('low_signal')|ev.copy_call.eq('diploid'))
    high=(ev.infer_call.eq('elevated_signal')|ev.copy_call.eq('aneuploid'))
    ev['compatible_native_entry']=fib & ev.copy_call.eq('diploid') & ~ev.infer_call.eq('elevated_signal') & ~db
    ev['signal_only_entry']=fib & ev.infer_call.eq('low_signal') & ev.copy_call.eq('not_evaluable') & ~db
    ev['native_entry_calibrated']=ev.compatible_native_entry & truth(ev.reference_QC_pass)
    ev['native_entry_reference_failed']=ev.compatible_native_entry & ~truth(ev.reference_QC_pass)
    ev['fib_CNV_conflict']=fib & low & high & ~db
    ev['fib_CNV_high_no_low']=fib & high & ~low & ~db
    ev['fib_CNV_unassessed']=fib & ~high & ~low & ~db
    ev['not_ready_for_formal_admission']=True
    ev['clinical_Bone_Lung_group']='unverified'
    ev['review_scope']='preprocessed_subset_not_full_called'
    nums={k:int(truth(ev[k]).sum()) for k in ['supported_CAF_candidate','expanded_Fib_expression_candidate','compatible_native_entry','signal_only_entry','native_entry_calibrated','native_entry_reference_failed','fib_CNV_conflict','fib_CNV_high_no_low','fib_CNV_unassessed']}
    z={k:row[k] for k in ['sample','dataset','patient_id','site','technology','treatment_raw','disease_stage_raw']}
    z.update(input_cells=len(ev),phase27_fibro_winner=int(fib.sum()),**nums)
    m=ev.compatible_native_entry
    z.update(native_corrected_Fib_ge2=int((m & ev.corrected_det_Fib.ge(2)).sum()),native_corrected_lineage_supported=int((m & truth(ev.corrected_lineage_supported)).sum()),native_mixed_flag=int((m&truth(ev.mixed_lineage_flag)).sum()),native_CTA_context=int((m&truth(ev.CTA_context)).sum()),native_osteo_program=int((m&ev.raw_det_Osteo.ge(3)&ev.raw_pct_Osteo.gt(.2)).sum()),native_mural_program=int((m&ev.raw_det_Mural.ge(3)&ev.raw_pct_Mural.gt(.2)).sum()))
    summaries.append(z)
    for cl,g in ev.groupby('cluster',sort=False):
        if not (g.compatible_native_entry.any() or g.signal_only_entry.any()):continue
        gg=g[g.compatible_native_entry]
        zz={'sample':sm,'cluster':cl,'original_lineage':g.lineage_phase3.iloc[0],'cluster_cells':len(g),'native_entry':int(g.compatible_native_entry.sum()),'signal_only':int(g.signal_only_entry.sum()),'reference_pass':bool(truth(g.reference_QC_pass).all()),'any_high_fraction':float(g.cluster_any_high_fraction.iloc[0]) if pd.notna(g.cluster_any_high_fraction.iloc[0]) else None,'CTA_context':bool(truth(g.CTA_context).any())}
        for k in ['Fib','Osteo','Mural','CTA']:
            zz['native_raw_det_'+k+'_median']=float(gg['raw_det_'+k].median()) if len(gg) else None
        for k in keys:zz['native_mean_'+k]=float(gg['phase27_'+k].mean()) if len(gg) else None
        clusters.append(zz)
    for scope,mask in [('native',ev.compatible_native_entry.to_numpy()),('calibrated',ev.native_entry_calibrated.to_numpy())]:
        if not mask.any():continue
        for j,g in enumerate(genes):dots.append({'sample':sm,'scope':scope,'gene':g,'cells':int(mask.sum()),'pct_detected':float((raw[mask,j]>=1).mean()*100),'mean_logCP10k':float(log[mask,j].mean()),'measured':bool(mapping[g])})
    take=['cell_id','sample','dataset','patient_id','site','technology','cluster','lineage_phase3','phase27_panel_winner','n_genes','n_counts','infer_call','copy_call','reference_QC_pass','audit_doublet_positive','supported_CAF_candidate','expanded_Fib_expression_candidate','compatible_native_entry','native_entry_calibrated','native_entry_reference_failed','signal_only_entry','fib_CNV_conflict','fib_CNV_high_no_low','fib_CNV_unassessed','corrected_lineage_supported','corrected_det_Fib','decontx_contamination','CTA_context','mixed_lineage_flag','raw_det_Fib','raw_det_Osteo','raw_det_Mural','raw_pct_Osteo','raw_pct_Mural','not_ready_for_formal_admission','clinical_Bone_Lung_group','review_scope']+['phase27_'+k for k in keys]
    allrows.append(ev[take].reset_index(drop=True))
    print('AUDIT_SAMPLE',json.dumps(z,ensure_ascii=True),flush=True)
    del a,x,raw,log,ev;gc.collect()
cells=pd.concat(allrows,ignore_index=True);s=pd.DataFrame(summaries)
assert cells.cell_id.is_unique and len(cells)==73645
save(cells,'all_cells_review.tsv.gz');save(cells[cells.compatible_native_entry],'native_entry_review_candidates.tsv.gz')
save(s,'sample_expansion_review.tsv');save(pd.DataFrame(clusters),'candidate_cluster_review.tsv');save(pd.DataFrame(dots),'candidate_marker_expression.tsv');save(pd.DataFrame(coverage),'panel_coverage.tsv')
patients=s.groupby(['dataset','patient_id'],sort=False).agg(libraries=('sample','size'),cells=('input_cells','sum'),native_entry=('compatible_native_entry','sum'),calibrated=('native_entry_calibrated','sum'),signal_only=('signal_only_entry','sum')).reset_index()
patients['independent_across_studies']='unresolved';patients['Bone_Lung_admission']='site_unverified'
save(patients,'patient_capacity_review.tsv')
for z in source_stat:
    f=Path(z['path']);assert f.stat().st_size==z['size'] and f.stat().st_mtime_ns==z['mtime_ns']
after={str(f):sha(f) for f in tracked};assert before==after
dump('resources/INPUT_IMMUTABILITY.json',{'sha256_before':before,'sha256_after':after,'count_objects_stat':source_stat,'passed':True})
fig,ax=plt.subplots(figsize=(10,6));yy=np.arange(len(s))
ax.barh(yy,s.native_entry_calibrated,label='Native entry / reference pass',color='#277f8e')
ax.barh(yy,s.native_entry_reference_failed,left=s.native_entry_calibrated,label='Native entry / reference failed',color='#e0a13b')
ax.barh(yy,s.signal_only_entry,left=s.compatible_native_entry,label='Signal-only / native CNV unresolved',color='#a6aab2')
ax.set_yticks(yy,s['sample']);ax.invert_yaxis();ax.set_xlabel('Cells (audit candidates, not admitted patients)');ax.set_title('Expanded-layer feasibility on Phase88 retained cells')
ax.legend(loc='lower right',fontsize=8);fig.tight_layout();fig.savefig(O/'figures/expansion_capacity.png',dpi=160);plt.close(fig)
summary={'status':'REVIEW_COMPUTATION_COMPLETE_MANUAL_INTERPRETATION_PENDING','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_cells':len(cells),'libraries':len(s),'within_source_donor_labels':len(patients),'native_entry_cells':int(cells.compatible_native_entry.sum()),'native_reference_pass_cells':int(cells.native_entry_calibrated.sum()),'native_reference_failed_cells':int(cells.native_entry_reference_failed.sum()),'signal_only_cells':int(cells.signal_only_entry.sum()),'native_entry_donor_labels':int(patients.native_entry.gt(0).sum()),'calibrated_native_donor_labels':int(patients.calibrated.gt(0).sum()),'new_formal_patients':0,'new_Bone_Lung_patients':0,'new_inferential_tests':False,'input_immutability_pass':True,'limits':['Existing preprocessed cells only','No native inferCNV HMM in Phase88','No new cluster integration or fixed-label assignment','Anatomical sites and cross-study deduplication unresolved']}
dump('COMPUTATION_COMPLETE.json',summary)
print(json.dumps(summary),flush=True)
