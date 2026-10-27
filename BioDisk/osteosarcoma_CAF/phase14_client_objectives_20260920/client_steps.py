"""SYNTAX_CHECK: python -m py_compile client_steps.py
RUN_COMMAND: python -u client_steps.py
EXPECTED_SUCCESS_MARKER: ANALYSIS_COMPLETE and exit_code=0
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import json,itertools,hashlib
import numpy as np,pandas as pd,anndata as ad
P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase10_batch_gates_20260919';O=P/'phase14_client_objectives_20260920'
for d in ['tables','resources','logs']: (O/d).mkdir(parents=True,exist_ok=True)
def save(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
read=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
gate=read(B/'tables/sample_gate_changes.tsv');geo=read(P/'geo_preprocess_20260918/manifest.tsv').set_index('sample');old=read(P/'phase7_review_states_20260917/review/clinical_enriched_20260917.tsv')
early=read(P/'phase3_identity_20260914/delivery/source_sample_manifest.tsv').set_index('sample')
clinical=[]
early_subject={'GSM8337654':'NCHS-009','GSM8337655':'NCHS-009','GSM8337656':'NCHS-015','GSM8337657':'NCHS-015','GSM8337658':'SC069','GSM8337659':'SC072','GSM8337660':'SC073'}
for r in gate.to_dict('records'):
    sm=r['sample'];ds=r['dataset'];patient='';site='';treatment='not_reported';stage='';source='';group='Unknown';reason='lung_metastasis_history_missing';history='not_reported';status='unresolved'
    if ds in ['GSE152048','GSE299023']:
        short=sm.split('__',1)[1];c=old[(old.dataset==ds)&(old['sample']==short)].iloc[0];patient=str(c.patient_id);site=str(c.sampling_site_verified);treatment=str(c.preoperative_treatment_raw);stage=str(c.disease_stage_raw);source=str(c.source)
        if str(c.clinical_group_verified)=='B2':group='B2';history='yes';reason='explicit_lung_metastasis_specimen';status='verified'
        if short=='BC17':treatment+='; anti-PD-L1 per Methods'
    elif ds=='GSE270231':
        patient=early_subject[sm];site='Lung';group='B2';history='yes';reason='explicit_lung_metastasis_specimen';status='verified';stage='Lung metastasis';source='GSE270231; source_sample_manifest; BioSample mapping verified 20260914'
    else:
        c=geo.loc[sm];patient=str(c.patient_id);site=str(c.site);stage=str(c.disease_stage_raw);treatment=str(c.treatment_raw);source=str(c.source_evidence)
        if ds=='GSE162454':treatment='no chemotherapy or radiotherapy before resection (original Methods)';source+='; Front Oncol 2021 doi:10.3389/fonc.2021.709210'
        if ds=='GSE318840' and site.lower()=='lung':group='B2';history='yes';reason='explicit_lung_metastasis_specimen';status='verified'
        if ds=='GSE318840' and site.lower()=='soft_tissue':group='Other_recurrent_soft_tissue';reason='not_bone_primary_or_lung';status='verified_site'
        if ds=='GSE250015' and sm.endswith('GSM7969979'):group='A_reported_nonmetastasis';history='no_reported_metastasis';reason='GEO primary bone Non-metastasis; timing/followup unavailable';status='reported_status_only';source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7969979'
        if ds=='GSE250015' and sm.endswith('GSM7969982'):group='Metastatic_primary_lung_unconfirmed';history='metastasis_yes_lung_not_confirmed';reason='GEO primary bone Metastasis; metastatic organ not established';status='not_B1_verified';source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7969982'
    clinical.append(dict(**r,patient_verified=patient,donor_record_verified=ds+'|'+patient,site_verified=site,stage_raw=stage,treatment_evidence=treatment,clinical_group=group,lung_history=history,group_status=status,group_reason=reason,source_evidence=source,sampling_date='not_established',prior_surgery_vs_sampling='not_separated_in_current_sources',cross_study_independence='not_fully_verified',bone_lung_pair='not_established'))
c=pd.DataFrame(clinical);assert len(c)==40 and c['sample'].is_unique and c.unified_strict.sum()==3976
save(c,'01_clinical_sample_audit.tsv')
summary=c.groupby(['dataset','clinical_group']).agg(sample_records=('sample','size'),strict_CAF=('unified_strict','sum'),records_with_CAF=('unified_strict',lambda x:int((x>0).sum()))).reset_index();save(summary,'01_group_availability.tsv')
pair=c.groupby('donor_record_verified').agg(samples=('sample',lambda x:';'.join(x)),sites=('site_verified',lambda x:';'.join(sorted(set(x)))),records=('sample','size')).reset_index();save(pair,'01_patient_record_links.tsv')
design=pd.DataFrame([
 {'comparison':'A_vs_B1','status':'not_estimable','reason':'reported A sample has zero strict CAF; no verified lung-metastatic bone-primary B1'},
 {'comparison':'A_vs_B2','status':'not_estimable','reason':'no reported-A strict CAF available'},
 {'comparison':'paired_B1_B2','status':'not_estimable','reason':'no verified same-patient bone-lung pair'},
 {'comparison':'GSE152048_lung_vs_bone_primary','status':'descriptive_site_only','reason':'BC10 versus BC2 BC5 BC21 BC22; all preoperative chemotherapy; bone lung history unknown; one lung donor'},
 {'comparison':'GSE318840_lung_vs_treated_bone','status':'descriptive_site_only','reason':'GSM9504304 versus GSM9504300; one donor each, no inference; lung history of bone donor unknown'}])
save(design,'01_comparison_design.tsv')
(O/'STEP1_COMPLETE').write_text('CLINICAL_AUDIT_COMPLETE_A_B1_COMPARISONS_UNAVAILABLE\n')
print(summary.to_string(index=False),flush=True)

panels={'myCAF_customer':['ACTA2','TAGLN','MYH11','MYL9','TPM2','CNN1'],'iCAF_customer':['IL6','CXCL12','PDGFRA','DPT','HAS1'],'apCAF_customer':['CD74','HLA-DRA','HLA-DRB1'],
 'Osteogenic_customer':['RUNX2','SP7','DLX5','MSX2','ALPL','BGLAP','IBSP','SPP1','SPARC','COL1A1','SOX9','BMP2','BMP4','PHEX','DMP1','SOST'],
 'Mineral_customer':['MGP','ENPP1','STC1','ANXA2','ANXA5','CALM1','CALM2','CACNA1C','SLC8A1','PTH1R','VDR','CASR','TRPM7','TNFRSF11B','TNFSF11']}
surfaces=json.loads((P/'phase8_original_plan_20260919/resources/surface_candidates_frozen.json').read_text())['surface_candidates']
genes=sorted(set(sum(panels.values(),[])+surfaces+['PTPRC','PECAM1','EPCAM','S100A8','S100A9','RGS5','CSPG4','NOTCH3']))
protocol={'clinical':'preserve unknown lung history; A reported-status only; no fake A/B1 or paired tests','panel_source':'client docx exact myCAF/iCAF/apCAF/osteogenic/mineral lists',
 'operational_support':'within frozen strict CAF, >=2 panel genes counts>=1 in BOTH raw and DecontX; >=80% measured; labels nonexclusive; >=3 gene sensitivity',
 'MYH11':'included in myCAF; full mural reference evaluated separately; not automatically promoted to CAF',
 'group_differences':'only verified tissue site descriptive comparisons; patient pseudobulk log2(CP10k+0.1); no inference with one lung donor per study',
 'surface_screen':'16 existing surface candidates and all pairs; unaltered source raw counts >=1 or>=2; PTPRC PECAM1 EPCAM all<1 RNA exclusion; all QC cells as denominator',
 'target':'customer-panel-supported myCAF within strict CAF, not a new or metastasis-validated subtype','rank':'donor record summaries; RNA known-target fraction and coverage, CNV-risk and mural fractions retained; no FACS purity claim',
 'completion':'clinical absence is a data dependency; marker screening is preparatory until metastatic target exists'}
(O/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2));(O/'resources/customer_panels.json').write_text(json.dumps(panels,indent=2))
manifest=read(B/'resources/object_manifest.tsv');cellrows=[];panelrows=[];generows=[];targets={}
for row in manifest.to_dict('records'):
    b=ad.read_h5ad(row['path']);ob=b.obs.copy();sy=b.var.symbol.astype(str).to_numpy();raw=np.full((len(b),len(genes)),np.nan);cor=raw.copy()
    for j,g in enumerate(genes):
        ix=np.flatnonzero(sy==g)
        if len(ix):raw[:,j]=np.asarray(b.layers['counts'][:,ix].sum(1)).ravel();cor[:,j]=np.asarray(b.layers['decontx_counts'][:,ix].sum(1)).ravel()
    masks={};out=ob[['dataset','sample','unified_strict','mural_reference','osteo_sensitive','competitor_sensitive']].copy()
    for name,gs in panels.items():
        ix=[genes.index(g) for g in gs];coverage=float(np.isfinite(raw[:,ix]).all(0).mean());nr=(raw[:,ix]>=1).sum(1);nc=(cor[:,ix]>=1).sum(1);supported=(nr>=2)&(nc>=2)&(coverage>=.8)
        out[name+'_raw_genes']=nr;out[name+'_corrected_genes']=nc;out[name+'_support']=supported;out[name+'_coverage']=coverage;masks[name]=supported
        for label,mask in [('strict',ob.unified_strict.to_numpy()),('mural_reference',ob.mural_reference.to_numpy()),('osteo_sensitive',ob.osteo_sensitive.to_numpy())]:
            panelrows.append(dict(sample=row['sample'],dataset=row['dataset'],set=label,panel=name,cells=int(mask.sum()),coverage=coverage,supported2=int((supported&mask).sum()),supported3=int(((nr>=3)&(nc>=3)&mask&(coverage>=.8)).sum())))
    my=ob.unified_strict.to_numpy()&masks['myCAF_customer'];targets[row['sample']]=set(ob.index[my])
    out['customer_myCAF_target']=my;out['MYH11_raw']=raw[:,genes.index('MYH11')];out['MYH11_corrected']=cor[:,genes.index('MYH11')]
    cellrows.append(out.reset_index(names='cell_id'))
    sm=row['sample'];cm=c.set_index('sample').loc[sm]
    for label,mask in [('strict_all',ob.unified_strict.to_numpy()),('customer_myCAF',my),('strict_other',ob.unified_strict.to_numpy()&~my)]:
        for kind,mat,den_col in [('raw',raw,'full_library_counts'),('corrected',cor,'corrected_library_counts')]:
            lib=float(ob.loc[mask,den_col].sum())
            for j,g in enumerate(genes):
                measured=bool(np.isfinite(mat[:,j]).all());val=float(mat[mask,j].sum()) if measured else np.nan
                generows.append(dict(sample=sm,dataset=row['dataset'],donor_record=cm.donor_record_verified,clinical_group=cm.clinical_group,site=cm.site_verified,set=label,expression=kind,gene=g,cells=int(mask.sum()),measured=measured,library_total=lib,gene_counts=val,CP10k=val/lib*1e4 if lib>0 and measured else np.nan,detection_fraction=float((mat[mask,j]>=1).mean()) if mask.any() and measured else np.nan))
    print('PANELS',sm,int(my.sum()),flush=True)
cells=pd.concat(cellrows);pr=pd.DataFrame(panelrows);gr=pd.DataFrame(generows)
save(cells,'02_customer_panel_cell_evidence.tsv.gz');save(pr,'02_panel_sample_counts.tsv');save(gr,'02_patient_gene_expression.tsv.gz')
assert cells.loc[cells.unified_strict].shape[0]==3976
# Descriptive group contrasts use patient/sample pseudobulk, not pooled cell tests.
contrasts=[]
for ds,lung,bones in [('GSE152048',['GSE152048__BC10'],['GSE152048__BC2','GSE152048__BC5','GSE152048__BC21','GSE152048__BC22']),('GSE318840',['GSE318840__GSM9504304'],['GSE318840__GSM9504300'])]:
    for label in ['strict_all','customer_myCAF']:
        for kind in ['raw','corrected']:
            sub=gr[(gr.dataset==ds)&(gr['set']==label)&(gr.expression==kind)&(gr.cells>=20)&gr.measured]
            for gene in genes:
                q=sub[sub.gene==gene];l=q[q['sample'].isin(lung)];b=q[q['sample'].isin(bones)]
                if not len(l) or not len(b):continue
                contrasts.append(dict(dataset=ds,set=label,expression=kind,gene=gene,lung_records=len(l),bone_primary_records=len(b),lung_samples=';'.join(l['sample']),bone_samples=';'.join(b['sample']),lung_median_log2CP10k=float(np.log2(l.CP10k+.1).median()),bone_median_log2CP10k=float(np.log2(b.CP10k+.1).median()),delta_log2CP10k=float(np.log2(l.CP10k+.1).median()-np.log2(b.CP10k+.1).median()),interpretation='tissue_site_descriptive_not_A_vs_B_or_metastatic_risk'))
save(pd.DataFrame(contrasts),'02_site_descriptive_gene_effects.tsv')
(O/'STEP2_COMPLETE').write_text('CUSTOMER_PANELS_REVIEWED_FORMAL_METASTASIS_COMPARISON_UNAVAILABLE\n')
# Full-cell surface gates are evaluated after target definition; no gate selected to maximize metastasis results.
combos=[(g,) for g in surfaces]+list(itertools.combinations(surfaces,2));screen=[];sm_manifest=read(B/'resources/source_manifest.tsv')
for row in sm_manifest.to_dict('records'):
    sm=row['sample'];a=ad.read_h5ad(row['source_file']);ev=read(B/'tables'/(sm+'_unified_evidence.tsv.gz')).set_index('cell_id').loc[a.obs_names];cm=c.set_index('sample').loc[sm]
    assert len(ev)==len(a);sy=a.var.symbol.astype(str).to_numpy();gs=surfaces+['PTPRC','PECAM1','EPCAM'];x=np.full((len(a),len(gs)),np.nan)
    for j,g in enumerate(gs):
        ix=np.flatnonzero(sy==g)
        if len(ix):x[:,j]=np.asarray(a.layers['counts'][:,ix].sum(1)).ravel()
    dump=np.isfinite(x[:,-3:]).all(1)&(x[:,-3:]<1).all(1);target=np.array([v in targets.get(sm,set()) for v in a.obs_names]);strict=ev.unified_strict.to_numpy();mural=ev.mural_reference.to_numpy();risk=(ev.infer_call.eq('elevated_signal')|ev.cnv2_call.eq('aneuploid')).to_numpy();unknown=~(strict|mural|risk)
    for threshold in [1,2]:
        for pair in combos:
            ix=[gs.index(g) for g in pair]
            if not np.isfinite(x[:,ix]).all():continue
            gm=dump&(x[:,ix]>=threshold).all(1)
            screen.append(dict(sample=sm,dataset=cm.dataset,donor_record=cm.donor_record_verified,clinical_group=cm.clinical_group,markers='+'.join(pair),count_threshold=threshold,all_QC_cells=len(a),target_cells=int(target.sum()),gate_cells=int(gm.sum()),target_in_gate=int((gm&target).sum()),strict_CAF_in_gate=int((gm&strict).sum()),mural_in_gate=int((gm&mural).sum()),CNV_risk_in_gate=int((gm&risk).sum()),other_or_unresolved_in_gate=int((gm&unknown).sum())))
    print('SCREEN',sm,len(a),flush=True)
scn=pd.DataFrame(screen);save(scn,'03_surface_screen_sample_counts.tsv.gz')
cols=['all_QC_cells','target_cells','gate_cells','target_in_gate','strict_CAF_in_gate','mural_in_gate','CNV_risk_in_gate','other_or_unresolved_in_gate']
don=scn.groupby(['dataset','donor_record','clinical_group','markers','count_threshold'])[cols].sum().reset_index()
for name,num,den in [('target_coverage','target_in_gate','target_cells'),('known_target_fraction','target_in_gate','gate_cells'),('known_CAF_fraction','strict_CAF_in_gate','gate_cells'),('CNV_risk_fraction','CNV_risk_in_gate','gate_cells'),('mural_fraction','mural_in_gate','gate_cells')]:don[name]=don[num]/don[den].replace(0,np.nan)
save(don,'03_surface_screen_donor.tsv.gz')
rank=[]
for scope in ['all_target_records','B2_target_records']:
    sub=don[(don.target_cells>=20)&don.dataset.ne('GSE250015')]
    if scope=='B2_target_records':sub=sub[sub.clinical_group=='B2']
    for (markers,t),g in sub.groupby(['markers','count_threshold']):
        rank.append(dict(scope=scope,markers=markers,count_threshold=t,eligible_donor_records=len(g),studies=g.dataset.nunique(),records_with_gate=int((g.gate_cells>0).sum()),median_target_coverage=g.target_coverage.median(),median_known_target_fraction=g.known_target_fraction.median(),median_known_CAF_fraction=g.known_CAF_fraction.median(),worst_CNV_risk_fraction=g.CNV_risk_fraction.max(),worst_mural_fraction=g.mural_fraction.max(),score=float((g.target_coverage*g.known_target_fraction.fillna(0)).median()),interpretation='RNA preparatory myCAF screen; not metastasis-linked subtype or protein sorting purity'))
rr=pd.DataFrame(rank).sort_values(['scope','count_threshold','score'],ascending=[True,True,False]);save(rr,'03_surface_candidate_ranking.tsv')
(O/'STEP3_SCREEN_COMPLETE').write_text('PREPARATORY_SURFACE_SCREEN_COMPLETE_TARGET_METASTASIS_LINK_UNRESOLVED\n')
status={'status':'CLIENT_THREE_STEPS_COMPUTED_REVIEW_PENDING','samples':40,'strict_CAF':3976,'customer_myCAF_supported':int(cells.customer_myCAF_target.sum()),'A_strict_CAF':int(c.loc[c.clinical_group=='A_reported_nonmetastasis','unified_strict'].sum()),'verified_B1_samples':0,'verified_bone_lung_pairs':0,'clinical_metastasis_group_DE_performed':False,'surface_screen_QC_cells':int(gate.retained_cells.sum()),'final_metastasis_specific_panel_delivered':False}
(O/'STATUS.json').write_text(json.dumps(status,indent=2));(O/'ANALYSIS_COMPLETE').write_text('ANALYSIS_COMPLETE\n');print(json.dumps(status),flush=True)
