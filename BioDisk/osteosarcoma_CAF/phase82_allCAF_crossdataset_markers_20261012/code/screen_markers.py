from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,traceback
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse

P=Path('/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF');B=P/'phase81_allCAF_bone_lung_20261012';O=P/'phase82_allCAF_crossdataset_markers_20261012'
for d in ['code','tables','resources','QA','figures','delivery','logs']:(O/d).mkdir(parents=True,exist_ok=True)
def read(p,**kw):return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def save(d,name):d.to_csv(O/'tables'/name,sep='\t',index=False,na_rep='NA')
def dump(v,p):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def truth(s):return s.astype(str).str.lower().isin(['true','1','1.0'])
def now():return datetime.now(timezone.utc).isoformat()

def main():
    assert not (O/'SCREEN_COMPLETE.json').exists()
    protocol=dict(frozen_at=now(),user_request='在这个层面重新寻找跨数据集支持的显著marker',
      level='All CAF candidate cells pooled within each original layer; anatomical Bone vs Lung',
      primary_rule='Same candidate layer, same direction, valid independent patient-group raw P<0.05 in at least two datasets; both increased and decreased candidates retained',
      layers='Supported and expanded overlap; never counted as separate datasets or validation cohorts',
      secondary_rule='Cross-layer same-direction nominal replication reported separately; insufficient-replication sources provide direction only',
      family='Preserve complete Phase81 original P/FDR/family; no new combined or intersection-selected P value',
      canonical='Use frozen per-source canonical feature-to-symbol registry; do not choose transcript by minimum P',
      sensitivity='Descriptive leave-one-patient-out effect, paired patient direction, existing doublet-excluded P/FDR, fibro-lineage-restricted patient effects and detection; no new screening cutoffs',
      flow='Membrane annotation appended after gene screen, no antibody or membrane filter',
      claim='Cross-dataset nominally replicated candidate expression markers, not independently validated cell-surface assay markers or CAF-specific identity markers',
      independence='Source-scoped donor IDs and sample IDs checked; absence of repeated identifiers alone does not establish independent recruitment')
    dump(protocol,O/'resources/PROTOCOL_FROZEN.json')
    paths=dict(stats=B/'tables/all_genes_statistics.tsv.gz',canonical=P/'phase58_F1F5_round1_20261004/resources/source_canonical_feature_map.tsv',clinical=B/'tables/clinical_sample_manifest.tsv',roster=B/'tables/patient_expression_roster.tsv',counts=P/'phase41_full_transcriptome_setup_20261001/objects/all_candidate_full_genes.h5ad',surface=P/'phase44_patient_mechanism_20261001/resources/intercell/receptor_roles.tsv')
    inputs={k:dict(path=str(f),sha256=sha(f),bytes=f.stat().st_size) for k,f in paths.items()};dump(inputs,O/'resources/INPUT_MANIFEST.json')
    allstats=read(paths['stats']);can=read(paths['canonical'],keep_default_na=False)
    assert not can.duplicated(['dataset','symbol']).any()
    mapping=can[['dataset','feature','symbol']].drop_duplicates();assert not mapping.duplicated(['dataset','feature']).any()
    stats=allstats.drop(columns='symbol').merge(mapping,on=['dataset','feature'],how='inner',validate='many_to_one')
    assert not stats.duplicated(['sensitivity','layer','dataset','design','symbol']).any()
    main=stats[(stats.sensitivity=='all_cells')&(stats.design=='independent')]
    nominal=main[main.P.lt(.05)].copy();nominal['direction']=np.sign(nominal.mean_lung_minus_bone).astype(int)
    registry=[]
    for (layer,symbol,direction),z in nominal.groupby(['layer','symbol','direction']):
        if z.dataset.nunique()<2:continue
        registry.append(dict(layer=layer,symbol=symbol,direction='lung_higher' if direction>0 else 'bone_higher',datasets=';'.join(sorted(z.dataset.unique())),nominal_datasets=z.dataset.nunique(),FDR_supported_datasets=int(z.FDR.lt(.05).sum()),minimum_abs_effect=float(z.mean_lung_minus_bone.abs().min()),maximum_P=float(z.P.max()),feature_ids=';'.join(sorted(z.feature.unique()))))
    primary=pd.DataFrame(registry).sort_values(['layer','direction','symbol']);save(primary,'primary_same_layer_crossdataset_markers.tsv')
    cross=[];conflicts=[]
    for symbol,z in nominal.groupby('symbol'):
        for i,a in enumerate(z.itertuples()):
            for b in list(z.itertuples())[i+1:]:
                if a.dataset==b.dataset:continue
                row=dict(symbol=symbol,dataset1=a.dataset,layer1=a.layer,effect1=a.mean_lung_minus_bone,P1=a.P,FDR1=a.FDR,dataset2=b.dataset,layer2=b.layer,effect2=b.mean_lung_minus_bone,P2=b.P,FDR2=b.FDR,same_direction=a.mean_lung_minus_bone*b.mean_lung_minus_bone>0)
                if not row['same_direction']:conflicts.append(row)
                elif a.layer!=b.layer:cross.append(row)
    cross=pd.DataFrame(cross);conflicts=pd.DataFrame(conflicts);save(cross,'secondary_cross_layer_replication.tsv');save(conflicts,'crossdataset_nominal_direction_conflicts.tsv')
    candidate_symbols=sorted(set(primary.symbol)|set(cross.symbol)|set(conflicts.symbol))
    evidence=stats[stats.symbol.isin(candidate_symbols)].copy();save(evidence,'all_candidate_source_layer_design_statistics.tsv')
    cm=read(paths['clinical'],keep_default_na=False)
    sources=sorted(set(';'.join(primary.datasets).split(';')))
    independence=[]
    for i,ds1 in enumerate(sources):
        for ds2 in sources[i+1:]:
            a=cm[cm.dataset.eq(ds1)];b=cm[cm.dataset.eq(ds2)]
            independence.append(dict(dataset1=ds1,dataset2=ds2,shared_sample_ids=';'.join(sorted(set(a['sample'])&set(b['sample']))),shared_recorded_patient_ids=';'.join(sorted(set(a.patient_id)-{'','Unknown','NA'} & (set(b.patient_id)-{'','Unknown','NA'}))),claim='Two source datasets; unlinked identifiers alone do not prove recruitment independence'))
    save(pd.DataFrame(independence),'dataset_identity_overlap_check.tsv')
    donors=read(paths['roster']);cached={};patient=[];sensitivity=[]
    for (layer,ds,design),e in evidence[evidence.sensitivity.eq('all_cells')].groupby(['layer','dataset','design']):
        key=f'all_cells__{layer}__{ds}__All_CAF__{design}'
        with np.load(B/'objects'/f'{key}.npz',allow_pickle=True) as z:
            features=z['features'].astype(str);x=z['bone_expression'];y=z['lung_expression'];b=z['bone_patients'].astype(str);l=z['lung_patients'].astype(str)
        cached[key]=(features,x,y,b,l)
        for r in e.itertuples():
            j=int(np.flatnonzero(features==r.feature)[0]);bx=x[:,j];ly=y[:,j];effect=r.mean_lung_minus_bone
            if design=='independent':
                vals=[np.nanmean(ly)-np.nanmean(np.delete(bx,i)) for i in range(len(bx)) if np.isfinite(np.delete(bx,i)).sum()>0]
                vals += [np.nanmean(np.delete(ly,i))-np.nanmean(bx) for i in range(len(ly)) if np.isfinite(np.delete(ly,i)).sum()>0]
                agree=sum(float(a)>float(v) for a in ly[np.isfinite(ly)] for v in bx[np.isfinite(bx)]);ties=sum(float(a)==float(v) for a in ly[np.isfinite(ly)] for v in bx[np.isfinite(bx)])
                denom=np.isfinite(ly).sum()*np.isfinite(bx).sum();auc=(agree+.5*ties)/denom if denom else np.nan
                paired_same=paired_opposite=0
            else:
                diff=ly-bx;vals=[np.nanmean(np.delete(diff,i)) for i in range(len(diff)) if np.isfinite(np.delete(diff,i)).sum()>0];auc=np.nan;paired_same=int((diff>0).sum());paired_opposite=int((diff<0).sum())
            vals=np.array(vals,float);vals=vals[np.isfinite(vals)]
            sensitivity.append(dict(layer=layer,dataset=ds,design=design,symbol=r.symbol,feature=r.feature,bone_N=r.bone_N,lung_N=r.lung_N,effect=effect,P=r.P,FDR=r.FDR,LOPO_effect_min=float(vals.min()) if len(vals) else np.nan,LOPO_effect_max=float(vals.max()) if len(vals) else np.nan,LOPO_direction_consistent=bool(np.all(vals*effect>0)) if len(vals) else None,patient_rank_AUC_lung_higher=auc,paired_lung_higher_patients=paired_same,paired_bone_higher_patients=paired_opposite))
            for site,ids,values in [('Bone',b,bx),('Lung',l,ly)]:
                for donor,val in zip(ids,values):
                    nr=donors[(donors.sensitivity=='all_cells')&(donors.layer==layer)&(donors.dataset==ds)&(donors.donor_record==donor)&(donors.site==site)]
                    assert len(nr)==1
                    patient.append(dict(layer=layer,dataset=ds,design=design,symbol=r.symbol,feature=r.feature,donor_record=donor,site=site,pseudobulk_log2_CP10k_plus01=val,cells=int(nr.cells.iloc[0]),library_UMIs=float(nr.library_UMIs.iloc[0])))
    sensitivity=pd.DataFrame(sensitivity);patient=pd.DataFrame(patient);save(sensitivity,'candidate_patient_direction_sensitivity.tsv');save(patient,'candidate_patient_expression.tsv')
    # Per-cell detection and Fibro-restricted descriptive effects, never cell-level inferential P.
    a=ad.read_h5ad(paths['counts']);raw=a.X.tocsr();pos={f:i for i,f in enumerate(a.var_names)}
    presence=sparse.load_npz(P/'phase41_full_transcriptome_setup_20261001/objects/sample_feature_presence.npz').tocsr();sr=read(P/'phase41_full_transcriptome_setup_20261001/tables/feature_presence_sample_index.tsv').set_index('sample').presence_row.to_dict()
    cell_layers={layer:read(P/'phase73_F134_bone_lung_20261009/strict_v2/tables'/f'{layer}_fixed_cells_with_sites.tsv.gz',keep_default_na=False) for layer in ['supported','expanded']}
    detection=[]
    for (layer,ds,site,donor),z in patient.groupby(['layer','dataset','site','donor_record']):
        ob=cell_layers[layer];ob=ob[(ob.dataset==ds)&(ob.strict_site==site)&(ob.donor_record==donor)]
        for subset,qq in [('all_candidates',ob),('Fibro_only',ob[ob.review_lineage_code.eq('Fibro')])]:
            ix=a.obs_names.get_indexer(qq.cell_id);assert (ix>=0).all()
            library=float(raw[ix].sum()) if len(ix) else 0
            for r in z[['symbol','feature']].drop_duplicates().itertuples():
                j=pos[r.feature];ss=qq['sample'].unique();measured=bool(len(ss) and presence[[sr[s] for s in ss],j].sum()==len(ss))
                v=raw[ix,j].toarray().ravel() if len(ix) else np.array([])
                detection.append(dict(layer=layer,dataset=ds,site=site,donor_record=donor,subset=subset,symbol=r.symbol,feature=r.feature,cells=len(ix),measured=measured,detected_cells=int((v>0).sum()) if measured else None,detection_fraction=float(np.mean(v>0)) if measured else np.nan,library_UMIs=library,raw_gene_UMIs=float(v.sum()) if measured else np.nan,pseudobulk_log2_CP10k_plus01=float(np.log2(v.sum()/library*1e4+.1)) if measured and library else np.nan))
    detection=pd.DataFrame(detection);save(detection,'candidate_detection_and_Fibro_patient_expression.tsv')
    fibro=[]
    for r in sensitivity[sensitivity.design.eq('independent')].itertuples():
        pp=patient[(patient.layer==r.layer)&(patient.dataset==r.dataset)&(patient.design=='independent')&(patient.symbol==r.symbol)]
        zz=detection[(detection.layer==r.layer)&(detection.dataset==r.dataset)&(detection.symbol==r.symbol)&(detection.subset=='Fibro_only')]
        vals={}
        for side in ['Bone','Lung']:
            q=zz[(zz.site==side)&zz.donor_record.isin(pp[pp.site.eq(side)].donor_record)];vals[side]=q.pseudobulk_log2_CP10k_plus01.dropna()
        effect=vals['Lung'].mean()-vals['Bone'].mean() if len(vals['Bone']) and len(vals['Lung']) else np.nan
        fibro.append(dict(layer=r.layer,dataset=r.dataset,symbol=r.symbol,Fibro_bone_N=len(vals['Bone']),Fibro_lung_N=len(vals['Lung']),Fibro_effect=effect,all_candidate_effect=r.effect,same_direction=bool(effect*r.effect>0) if np.isfinite(effect) else None,evidence='Descriptive lineage-restricted effect, not a new P or additional screening gate'))
    fibro=pd.DataFrame(fibro);save(fibro,'Fibro_restricted_direction_summary.tsv')
    surface=read(paths['surface'],keep_default_na=False);surf=surface[surface.genesymbol.isin(candidate_symbols)];save(surf,'candidate_existing_membrane_annotation.tsv')
    membrane_symbols=set(surf.loc[truth(surf.plasma_membrane_transmembrane)|truth(surf.plasma_membrane_peripheral),'genesymbol'])
    # One row per strongest same-layer candidate; exact per-source statistics retained.
    rows=[]
    for r in primary.itertuples():
        row=r._asdict();row.pop('Index',None)
        row['existing_membrane_annotation']=r.symbol in membrane_symbols
        for ds in sources:
            g=evidence[(evidence.layer==r.layer)&(evidence.dataset==ds)&(evidence.design=='independent')&(evidence.sensitivity=='all_cells')&(evidence.symbol==r.symbol)].iloc[0]
            sen=sensitivity[(sensitivity.layer==r.layer)&(sensitivity.dataset==ds)&(sensitivity.design=='independent')&(sensitivity.symbol==r.symbol)].iloc[0]
            dd=evidence[(evidence.layer==r.layer)&(evidence.dataset==ds)&(evidence.design=='independent')&(evidence.sensitivity=='exclude_doublets')&(evidence.symbol==r.symbol)].iloc[0]
            ff=fibro[(fibro.layer==r.layer)&(fibro.dataset==ds)&(fibro.symbol==r.symbol)].iloc[0]
            for k in ['bone_N','lung_N','mean_lung_minus_bone','P','FDR','family']:row[ds+'__'+k]=g[k]
            for k in ['LOPO_effect_min','LOPO_effect_max','LOPO_direction_consistent','patient_rank_AUC_lung_higher']:row[ds+'__'+k]=sen[k]
            row[ds+'__doublet_excluded_P']=dd.P;row[ds+'__doublet_excluded_effect']=dd.mean_lung_minus_bone
            for k in ['Fibro_bone_N','Fibro_lung_N','Fibro_effect','same_direction']:row[ds+'__'+k]=ff[k]
        supported=evidence[(evidence.layer=='supported')&(evidence.dataset=='SCPCP000017')&(evidence.design=='independent')&(evidence.sensitivity=='all_cells')&(evidence.symbol==r.symbol)]
        if len(supported):
            g=supported.iloc[0];row['supported17_effect']=g.mean_lung_minus_bone;row['supported17_P']=g.P;row['supported17_FDR']=g.FDR;row['supported17_nominal_same_direction']=bool(g.P<.05 and g.mean_lung_minus_bone*(1 if r.direction=='lung_higher' else -1)>0)
        pp=sensitivity[(sensitivity.layer==r.layer)&(sensitivity.dataset=='SCPCP000017')&(sensitivity.design=='paired')&(sensitivity.symbol==r.symbol)]
        if len(pp):
            g=pp.iloc[0];row['paired17_N']=g.bone_N;row['paired17_effect']=g.effect;row['paired17_lung_higher_patients']=g.paired_lung_higher_patients;row['paired17_bone_higher_patients']=g.paired_bone_higher_patients
        rows.append(row)
    registry=pd.DataFrame(rows);save(registry,'MARKER_EVIDENCE_REGISTRY.tsv')
    # Independent integrity checks reproduce primary membership directly from original P.
    expected=set()
    original_main=allstats[(allstats.sensitivity=='all_cells')&(allstats.design=='independent')&(allstats.layer=='expanded')]
    s17=original_main[(original_main.dataset=='SCPCP000017')&(original_main.P<.05)]
    s23=original_main[(original_main.dataset=='SCPCP000023')&(original_main.P<.05)]
    overlap=s17.merge(s23,on='feature',suffixes=('_17','_23'));overlap=overlap[overlap.mean_lung_minus_bone_17*overlap.mean_lung_minus_bone_23>0]
    expected=set(overlap.feature);observed=set(';'.join(primary[primary.layer.eq('expanded')].feature_ids).split(';'));assert observed==expected
    aliases=overlap[['feature','symbol_17','symbol_23']].merge(mapping[mapping.dataset.eq('SCPCP000017')][['feature','symbol']],on='feature',validate='one_to_one')
    save(aliases,'stable_feature_and_symbol_crosswalk.tsv')
    aa=detection[detection.subset.eq('all_candidates')]
    check=patient.merge(aa[['layer','dataset','site','donor_record','symbol','pseudobulk_log2_CP10k_plus01']],on=['layer','dataset','site','donor_record','symbol'],suffixes=('_saved','_reconstructed'),validate='many_to_one')
    error=(check.pseudobulk_log2_CP10k_plus01_saved-check.pseudobulk_log2_CP10k_plus01_reconstructed).abs().max();assert error<1e-9
    assert all(sha(paths[k])==info['sha256'] for k,info in inputs.items())
    receipt=dict(status='screen_complete_review_pending',completed_at=now(),primary_same_layer_markers=len(primary),primary_lung_higher=int(primary.direction.eq('lung_higher').sum()),primary_bone_higher=int(primary.direction.eq('bone_higher').sum()),cross_layer_support_records=len(cross),cross_layer_direction_conflicts=len(conflicts),total_candidate_symbols=len(candidate_symbols),primary_FDR_supported_datasets=int(primary.FDR_supported_datasets.sum()),source_statistics_unchanged=True,canonical_mapping_unique=True,primary_membership_reproduced=True,raw_patient_expression_checks=len(check),max_raw_patient_expression_error=float(error),new_inferential_tests=0)
    dump(receipt,O/'SCREEN_COMPLETE.json');dump(receipt,O/'QA/SCREEN_ACCEPTANCE.json')
    print(json.dumps(receipt,ensure_ascii=False),flush=True)
    cols=['symbol','direction','supported17_nominal_same_direction','paired17_effect','paired17_lung_higher_patients','paired17_bone_higher_patients','SCPCP000017__Fibro_effect','SCPCP000023__Fibro_effect','existing_membrane_annotation']
    print(registry[cols].to_string(index=False),flush=True)
    print('SECONDARY\n'+cross.to_string(index=False),flush=True)
    print('MEMBRANE_COLUMNS',surf.columns.tolist(),flush=True)
    print('LINEAGE_CODES',{k:v.review_lineage_code.value_counts().to_dict() for k,v in cell_layers.items()},flush=True)

if __name__=='__main__':
    try:main()
    except BaseException:
        (O/'logs/screen.exit_code').write_text('1\n');traceback.print_exc();raise
    else:(O/'logs/screen.exit_code').write_text('0\n')
